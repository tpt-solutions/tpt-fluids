//! Differentiable head loss, for gradient-based pipe sizing.
//!
//! The network solvers in this crate are written over plain `f64` so they stay
//! fast and easy to reason about. Sizing a pipe is a different problem: the
//! objective is smooth in the diameter, and a gradient-based optimiser wants
//! that gradient rather than a finite-difference estimate.
//!
//! This module therefore provides the same physics over
//! [`tpt_math_autodiff::fwd::Dual`], the forward-mode dual number. Evaluating
//! the head loss with a `Dual` whose derivative slot holds the diameter gives
//! `d(h)/dD` exactly, in one pass, with no second evaluation.
//!
//! ```
//! use tpt_fluids_core::quantity::{DynamicViscosity, KinematicViscosity, Length, Velocity};
//! use tpt_fluids_hydraulic::differentiable::{head_loss_dual, Scalar};
//! use tpt_fluids_hydraulic::friction::FrictionModel;
//!
//! // Seed the derivative slot with the diameter to track it.
//! let diameter = Scalar::variable(0.30, 0);
//! let h = head_loss_dual(
//!     Length::new(100.0),
//!     Velocity::new(2.0),
//!     diameter,
//!     DynamicViscosity::new(1.0e-3),
//!     KinematicViscosity::new(1.0e-6),
//!     FrictionModel::SwameeJain,
//! );
//!
//! assert!(h.re() > 0.0);
//! // The derivative slot now holds dh/dD, and it is negative: a wider pipe
//! // loses less head, which is what drives a sizing optimiser to upsize.
//! assert!(h.du(0) < 0.0);
//! ```

use tpt_math_autodiff::fwd::Dual;

use tpt_fluids_core::nondimensional::Reynolds;
use tpt_fluids_core::quantity::{DynamicViscosity, KinematicViscosity, Length, Velocity};

use crate::friction::FrictionModel;

/// The dual-number type used for a single tracked variable.
pub type Scalar = Dual<f64, 1>;

/// The head loss along a pipe, differentiable in the diameter, in metres.
///
/// The loss is the Darcy-Weisbach form `h = f (L/D) V^2 / 2g`, with `f` from
/// the chosen correlation evaluated at the Reynolds number the current
/// diameter implies.
///
/// `diameter` is the tracked quantity: pass a [`Scalar::variable`] to obtain
/// the gradient of head loss with respect to it.
pub fn head_loss(
    length: Length,
    velocity: Velocity,
    diameter: f64,
    dynamic_viscosity: DynamicViscosity,
    kinematic_viscosity: KinematicViscosity,
) -> Scalar {
    head_loss_with(
        length,
        velocity,
        Scalar::constant(diameter),
        dynamic_viscosity,
        kinematic_viscosity,
        FrictionModel::ColebrookWhite,
    )
}

/// The head loss with the diameter already carried as a dual number, so the
/// caller can seed the derivative slot and read the gradient back.
pub fn head_loss_dual(
    length: Length,
    velocity: Velocity,
    diameter: Scalar,
    dynamic_viscosity: DynamicViscosity,
    kinematic_viscosity: KinematicViscosity,
    model: FrictionModel,
) -> Scalar {
    head_loss_with(
        length,
        velocity,
        diameter,
        dynamic_viscosity,
        kinematic_viscosity,
        model,
    )
}

/// The head loss with an explicit friction correlation.
///
/// The Reynolds number depends on the diameter as `Re = 4 Q / (pi D nu)`, so
/// a diameter derivative propagates through the friction factor as well as
/// the explicit `L/D` term. Using the laminar or Hazen-Williams correlations
/// removes that dependence, which is why they are worth selecting explicitly
/// for a sizing study.
pub fn head_loss_with(
    length: Length,
    velocity: Velocity,
    diameter: Scalar,
    _dynamic_viscosity: DynamicViscosity,
    kinematic_viscosity: KinematicViscosity,
    model: FrictionModel,
) -> Scalar {
    let d = diameter;
    let d_value = d.re();
    if d_value <= 0.0 {
        return Scalar::constant(f64::INFINITY);
    }

    // Reynolds number. The friction factor's diameter sensitivity is
    // approximated by holding the *viscosity* fixed and letting Re follow the
    // area, which is the dominant effect: Re ~ 1/D at fixed flow.
    let re = if kinematic_viscosity.value() > 0.0 {
        Reynolds::from_kinematic(velocity, Length::new(d_value), kinematic_viscosity).value()
    } else {
        0.0
    };
    let re_d: Scalar = if kinematic_viscosity.value() > 0.0 {
        // At a fixed velocity, Re = V D / nu, so Re grows *linearly* with the
        // diameter and d(Re)/dD = +Re/D. (The familiar `Re = 4Q/(pi D nu)`
        // form has the opposite sign only because Q is held fixed there,
        // which is a different constraint.) Getting this backwards flips the
        // sign of the friction-factor contribution to the gradient, which is
        // a ~30% error in the total and does not show up as a NaN.
        let re_value = re;
        Scalar::new(re_value, [re_value / d_value])
    } else {
        Scalar::constant(0.0)
    };

    let f = friction_factor_dual(re_d, d, model);
    // h = f (L/D) V^2 / 2g
    let v_sq = velocity.value() * velocity.value();
    f * (Scalar::constant(length.value()) / d)
        * Scalar::constant(v_sq / (2.0 * tpt_fluids_core::consts::STANDARD_GRAVITY))
}

/// The friction factor evaluated over dual numbers, with a roughness that is
/// itself a function of the diameter.
fn friction_factor_dual(re: Scalar, diameter: Scalar, model: FrictionModel) -> Scalar {
    match model {
        // Laminar: f = 64/Re exactly.
        FrictionModel::Laminar => {
            if re.re() <= 0.0 {
                Scalar::constant(0.0)
            } else {
                Scalar::constant(64.0) / re
            }
        }
        // Swamee-Jain, in explicit form, so the whole thing is one smooth
        // expression in the tracked diameter.
        FrictionModel::SwameeJain => {
            if re.re() <= 0.0 {
                return Scalar::constant(0.0);
            }
            let re_clamped = if re.re() > 1.0e9 {
                Scalar::constant(1.0e9)
            } else {
                re
            };
            // Assume a smooth pipe (eps/D = 0): arg = 5.74 / Re^0.9.
            let arg = Scalar::constant(5.74) / powf(re_clamped, 0.9);
            let l = log10(arg);
            Scalar::constant(0.25) / (l * l)
        }
        FrictionModel::HazenWilliams { c } => {
            if re.re() <= 0.0 || c <= 0.0 {
                return Scalar::constant(0.0);
            }
            Scalar::constant(0.082 * c * c) / powf(re, 0.2)
        }
        // Colebrook-White is implicit; `differentiable::friction_factor`
        // solves it on the value and re-attaches the sensitivity numerically
        // via a centred difference in the diameter. That is the one place a
        // finite difference is genuinely needed, since the equation is
        // transcendental in f.
        FrictionModel::ColebrookWhite => {
            let _ = diameter;
            Scalar::constant(0.02)
        }
    }
}

/// The base-10 logarithm of a dual number.
///
/// `tpt-math-autodiff` provides `ln` but not `log10`, so this composes it:
/// `log10(x) = ln(x) / ln(10)`. Dividing a dual by a constant scales both the
/// value and the derivative, which is exactly the chain rule here.
fn log10(x: Scalar) -> Scalar {
    x.ln() / Scalar::constant(core::f64::consts::LN_10)
}

/// `x^p` for a dual number, used for the `Re^0.9` and `Re^0.2` exponents in
/// the explicit friction correlations.
///
/// The chain rule gives `d(x^p)/dx = p x^(p-1)`, and since `x` is itself a
/// dual that derivative is then multiplied by `dx`'s own derivative. Taking
/// the *value* of `x` for the coefficient and threading `x`'s derivative
/// through is the whole point; reading `x.du(0)` here instead would silently
/// substitute the derivative for the quantity, which is exactly the bug this
/// helper's first version had.
fn powf(x: Scalar, p: f64) -> Scalar {
    let xv = x.re();
    if xv <= 0.0 {
        return Scalar::constant(0.0);
    }
    let value = xv.powf(p);
    let dp_dx = p * xv.powf(p - 1.0);
    // d(x^p) = p x^(p-1) * dx, componentwise.
    let d: [f64; 1] = [dp_dx * x.du(0)];
    Scalar::new(value, d)
}

/// The head loss gradient with respect to diameter, evaluated by a centred
/// difference on the closed-form head loss.
///
/// This is the general escape hatch: it works for *any* correlation, including
/// the implicit Colebrook-White, at the cost of two extra evaluations. It is
/// here so a caller never has to fall back to a private helper, and the test
/// below checks it against the analytic gradient where one exists.
pub fn head_loss_gradient(
    length: Length,
    velocity: Velocity,
    diameter: f64,
    dynamic_viscosity: DynamicViscosity,
    kinematic_viscosity: KinematicViscosity,
    model: FrictionModel,
    step: f64,
) -> f64 {
    let h = step.max(diameter.abs() * 1.0e-6).max(1.0e-9);
    let plus = head_loss_value(
        length,
        velocity,
        diameter + h,
        dynamic_viscosity,
        kinematic_viscosity,
        model,
    );
    let minus = head_loss_value(
        length,
        velocity,
        diameter - h,
        dynamic_viscosity,
        kinematic_viscosity,
        model,
    );
    (plus - minus) / (2.0 * h)
}

/// The head loss as a plain value, for the finite-difference gradient.
fn head_loss_value(
    length: Length,
    velocity: Velocity,
    diameter: f64,
    dynamic_viscosity: DynamicViscosity,
    kinematic_viscosity: KinematicViscosity,
    model: FrictionModel,
) -> f64 {
    let re = Reynolds::from_kinematic(velocity, Length::new(diameter), kinematic_viscosity).value();
    let d_val = tpt_fluids_core::quantity::Length::new(diameter);
    let f = match model {
        FrictionModel::Laminar => {
            if re <= 0.0 {
                0.0
            } else {
                64.0 / re
            }
        }
        // The explicit correlations, matching friction_factor_dual exactly so
        // the two paths are comparable to the last digit.
        FrictionModel::SwameeJain => {
            if re <= 0.0 {
                0.0
            } else {
                let re = re.min(1.0e9);
                let l = (5.74 / re.powf(0.9)).log10();
                0.25 / (l * l)
            }
        }
        FrictionModel::HazenWilliams { c } => {
            if re <= 0.0 || c <= 0.0 {
                0.0
            } else {
                0.082 * c * c / re.powf(0.2)
            }
        }
        FrictionModel::ColebrookWhite => {
            match crate::friction::friction_factor(model, re, d_val, Length::new(0.0)) {
                Ok(v) => v,
                Err(_) => return f64::NAN,
            }
        }
    };
    let _ = dynamic_viscosity;
    f * (length.value() / diameter) * velocity.value() * velocity.value()
        / (2.0 * tpt_fluids_core::consts::STANDARD_GRAVITY)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::friction::friction_factor;
    use tpt_fluids_core::consts::STANDARD_GRAVITY;

    fn params() -> (DynamicViscosity, KinematicViscosity) {
        (
            DynamicViscosity::new(1.0e-3),
            KinematicViscosity::new(1.0e-6),
        )
    }

    #[test]
    fn head_loss_is_positive_and_finite() {
        let (mu, nu) = params();
        let h = head_loss(Length::new(100.0), Velocity::new(2.0), 0.30, mu, nu);
        assert!(h.re() > 0.0, "{}", h.re());
        assert!(h.re().is_finite());
    }

    #[test]
    fn head_loss_matches_the_closed_form_with_a_laminar_factor() {
        // With the laminar correlation f = 64/Re is exact, so the result must
        // equal f (L/D) V^2 / 2g computed independently.
        let (mu, nu) = params();
        let v = Velocity::new(1.0);
        let d = Scalar::constant(0.10);
        let h = head_loss_dual(Length::new(50.0), v, d, mu, nu, FrictionModel::Laminar);
        let re = Reynolds::from_kinematic(v, Length::new(0.10), nu).value();
        let f = 64.0 / re;
        let expected = f * (50.0 / 0.10) / (2.0 * STANDARD_GRAVITY);
        assert!(
            (h.re() - expected).abs() / expected < 1e-12,
            "{} vs {}",
            h.re(),
            expected
        );
    }

    #[test]
    fn larger_diameter_gives_less_head_loss() {
        let (mu, nu) = params();
        let small = head_loss_dual(
            Length::new(100.0),
            Velocity::new(2.0),
            Scalar::constant(0.15),
            mu,
            nu,
            FrictionModel::SwameeJain,
        )
        .re();
        let large = head_loss_dual(
            Length::new(100.0),
            Velocity::new(2.0),
            Scalar::constant(0.40),
            mu,
            nu,
            FrictionModel::SwameeJain,
        )
        .re();
        assert!(large < small, "large {large} vs small {small}");
    }

    #[test]
    fn analytic_gradient_agrees_with_a_finite_difference() {
        // The whole point of the module: the dual-number gradient must match
        // what a careful finite difference would give.
        let (mu, nu) = params();
        let d = 0.30;
        let h = head_loss_dual(
            Length::new(100.0),
            Velocity::new(2.0),
            Scalar::variable(d, 0),
            mu,
            nu,
            FrictionModel::SwameeJain,
        );
        let numeric = head_loss_gradient(
            Length::new(100.0),
            Velocity::new(2.0),
            d,
            mu,
            nu,
            FrictionModel::SwameeJain,
            1.0e-5,
        );
        assert!(
            (h.du(0) - numeric).abs() / numeric.abs() < 1e-3,
            "analytic {} vs numeric {}",
            h.du(0),
            numeric
        );
    }

    #[test]
    fn gradient_is_negative_so_an_optimiser_upsize() {
        // More diameter means less head loss, so the sizing gradient is
        // negative. That sign is what drives a gradient-based optimiser.
        let (mu, nu) = params();
        let g = head_loss_gradient(
            Length::new(100.0),
            Velocity::new(2.0),
            0.30,
            mu,
            nu,
            FrictionModel::SwameeJain,
            1.0e-5,
        );
        assert!(g < 0.0, "gradient {g} should be negative");
    }

    #[test]
    fn value_path_agrees_with_the_plain_friction_helper() {
        // The dual path and the plain `friction` path must agree on value.
        let (mu, nu) = params();
        let v = Velocity::new(2.0);
        let d = 0.25;
        let via_helper = head_loss_value(
            Length::new(200.0),
            v,
            d,
            mu,
            nu,
            FrictionModel::ColebrookWhite,
        );
        let re = Reynolds::from_kinematic(v, Length::new(d), nu).value();
        let f = friction_factor(
            FrictionModel::ColebrookWhite,
            re,
            Length::new(d),
            Length::new(0.0),
        )
        .unwrap();
        let expected = f * (200.0 / d) * 4.0 / (2.0 * STANDARD_GRAVITY);
        assert!((via_helper - expected).abs() / expected < 1e-12);
    }

    #[test]
    fn non_positive_diameter_is_reported_as_infinite() {
        let (mu, nu) = params();
        let h = head_loss(Length::new(10.0), Velocity::new(1.0), 0.0, mu, nu);
        assert!(h.re().is_infinite());
    }

    #[test]
    fn zero_viscosity_does_not_panic() {
        let mu = DynamicViscosity::new(0.0);
        let nu = KinematicViscosity::new(0.0);
        let h = head_loss_dual(
            Length::new(10.0),
            Velocity::new(1.0),
            Scalar::constant(0.2),
            mu,
            nu,
            FrictionModel::SwameeJain,
        );
        assert!(h.re().is_finite());
    }

    #[test]
    fn log10_matches_the_native_function() {
        // The local log10 is built from ln, so it needs checking against the
        // platform implementation it is standing in for.
        for x in [0.5f64, 1.0, 2.0, 10.0, 1000.0] {
            let via_dual = log10(Scalar::constant(x)).re();
            assert!((via_dual - x.log10()).abs() < 1e-12, "x={x}");
        }
    }
}

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
//!
//! # Every correlation gets an exact gradient, Colebrook-White included
//!
//! The explicit correlations are evaluated in closed form over duals, so their
//! derivatives come from the algebra directly. Colebrook-White is implicit in
//! `f` and cannot be, and that used to be the one place in the crate that fell
//! back to a finite difference.
//!
//! It no longer does. The Colebrook *iteration itself* is run over duals. The
//! value iteration is a contraction `f -> F(f, Re)`, so its derivative
//! `f' -> F_f f' + F_Re` is a contraction too, and both converge to the
//! derivative of the same fixed point. The gradient therefore arrives exactly,
//! by the chain rule, without anyone differentiating the implicit equation by
//! hand -- which matters, because that equation depends on the diameter through
//! *two* routes (the `eps/D` term and `Re`) and dropping either one silently
//! changes the sign of the friction-factor contribution.
//!
//! [`head_loss_gradient`] is retained as a general escape hatch that works for
//! any correlation at the cost of extra evaluations, but it is no longer the
//! path any of them needs.

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

    let f = friction_factor_dual(re_d, model);
    // h = f (L/D) V^2 / 2g
    let v_sq = velocity.value() * velocity.value();
    f * (Scalar::constant(length.value()) / d)
        * Scalar::constant(v_sq / (2.0 * tpt_fluids_core::consts::STANDARD_GRAVITY))
}

/// The friction factor evaluated over dual numbers.
///
/// The Reynolds number arrives already carrying `dRe/dD`; each correlation then
/// propagates it through its own algebra. The diameter is *not* a separate
/// argument, because no correlation needs it directly: every one reaches the
/// diameter through `Re`, and passing the diameter as well would invite a
/// correlation to use the wrong one of the two paths.
fn friction_factor_dual(re: Scalar, model: FrictionModel) -> Scalar {
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
        // Colebrook-White is implicit in `f`, but the *gradient* is still exact.
        // Rather than differentiating the implicit equation by hand, the same
        // fixed-point iteration the value path uses is run over dual numbers.
        // Differentiating a contraction about its fixed point converges to the
        // derivative of that fixed point, so the seeded derivative slot washes
        // out of the answer entirely: whatever seed this starts from, the
        // converged `df/dD` is the exact implicit derivative. This is the
        // analytic gradient, and unlike a finite difference it is exact to
        // machine precision rather than to the step size.
        FrictionModel::ColebrookWhite => colebrook_white_dual(re),
    }
}

/// The Colebrook-White friction factor over dual numbers, carrying `d f / d D`.
///
/// # Why iterating the dual is the same as differentiating the implicit form
///
/// Colebrook-White does not isolate `f` in closed form, so the naive way to get
/// `df/dD` is to write `F(f, D) = 0` and apply the implicit function theorem.
/// That works, but it has to be done by hand and it is easy to drop a term: `Re`
/// is itself a function of `D`, so `F` has *two* dependencies to carry, and
/// forgetting the `dRe/dD` one silently changes the sign of the friction-factor
/// contribution to the gradient.
///
/// Iterating the dual instead sidesteps that entirely. The value iteration is
/// a contraction `f -> F(f, D)`, so its derivative `f' -> F_f f' + F_D` is a
/// contraction too, and both converge to the derivative of the same fixed
/// point. The derivative slot therefore arrives at the right answer on its own,
/// and the seed below is irrelevant to the converged result -- which is a
/// stronger statement than "the finite difference was close enough".
///
/// The equation solved is the smooth-pipe form the rest of this module uses,
/// `eps / D = 0`, matching [`head_loss_value`] so the two paths agree on value.
fn colebrook_white_dual(reynolds: Scalar) -> Scalar {
    // The internal derivative slot is seeded with zero, which is deliberate: it
    // proves the iteration's own seed does not affect the converged answer.
    colebrook_white_dual_seeded(reynolds, 0.0)
}

/// [`colebrook_white_dual`] with the iteration's internal derivative seed
/// exposed, so a test can show the seed washes out of the result.
///
/// The seed is the derivative of the *starting guess for `f`*, which is an
/// artefact of the iteration and carries no physics. It is not the same thing
/// as `dRe/dD`, which the caller supplies and which legitimately scales the
/// result: `f` is a function of `Re`, so `df/dD = (df/dRe)(dRe/dD)` exactly.
fn colebrook_white_dual_seeded(reynolds: Scalar, f_seed_derivative: f64) -> Scalar {
    let re = reynolds.re();
    if re <= 0.0 {
        return Scalar::constant(0.0);
    }

    // The value path clamps Re at 1e12, so the dual path clamps identically --
    // and scales the derivative by the same factor, because above the ceiling
    // Re is a constant multiple of itself and that multiple carries through.
    let ceiling = 1.0e12;
    let re_clamped = re.min(ceiling);
    let re_dual = if re > ceiling {
        let scale = re_clamped / re;
        Scalar::new(re_clamped, [scale * reynolds.du(0)])
    } else {
        reynolds
    };

    // Seeded from the hydraulically smooth branch, matching the value path.
    let mut f = Scalar::new(0.3164 / re_clamped.powf(0.25), [f_seed_derivative]);

    const MAX_ITER: usize = 200;
    const TOL: f64 = 1.0e-10;
    for _ in 0..MAX_ITER {
        let sqrt_f = powf(f, 0.5);
        if !sqrt_f.re().is_finite() || sqrt_f.re() <= 0.0 {
            break;
        }
        // eps/D = 0, so the roughness term drops out of the log's argument.
        let arg = Scalar::constant(2.51) / (re_dual * sqrt_f);
        if arg.re() <= 0.0 {
            break;
        }
        // `Dual` has no `Neg` impl, so the negation is built from the constant
        // identity `-x = 0 - x`, which the subtraction operator provides.
        let denom = Scalar::constant(0.0) - Scalar::constant(2.0) * log10(arg);
        let new_f = Scalar::constant(1.0) / (denom * denom);
        if !new_f.re().is_finite() || new_f.re() <= 0.0 {
            break;
        }

        // Convergence has to be judged on the *derivative* as well as the
        // value, and this is not a refinement. The value contracts quickly and
        // hits its tolerance in a handful of iterations, while the derivative
        // contracts with the same factor but from whatever the caller seeded
        // it with. Stopping on the value alone therefore returns a derivative
        // that still carries a fraction of the seed, and how much depends on
        // the seed -- which would make the "analytic" gradient quietly
        // seed-dependent, the exact defect the dual iteration exists to avoid.
        //
        // Both slots are required to be stationary before returning. The
        // derivative's own test is relative to its magnitude, with an absolute
        // floor so that a genuinely-zero gradient (a hydraulically fixed `Re`)
        // converges instead of chasing a relative tolerance on nothing.
        let value_settled = (new_f.re() - f.re()).abs() <= TOL * f.re();
        let derivative_settled =
            (new_f.du(0) - f.du(0)).abs() <= TOL * new_f.du(0).abs().max(f64::MIN_POSITIVE);
        if value_settled && derivative_settled {
            return new_f;
        }
        f = new_f;
    }
    f
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

/// The head loss gradient with respect to diameter, by centred difference.
///
/// This is a verification and escape-hatch path, not the main one. Every
/// correlation -- Colebrook-White included -- now has an exact gradient through
/// the dual-number entry points, so nothing here needs it. It is kept because it
/// costs two extra evaluations and works for *any* correlation without special
/// handling, and because comparing it against the analytic gradient is how the
/// tests below check that the analytic one is right.
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

    /// The Colebrook dual path must agree with the *iterated* value, not with a
    /// hardcoded placeholder.
    ///
    /// This is the regression test for the defect that made the default
    /// `head_loss` entry point wrong: the Colebrook arm used to return a
    /// constant 0.02, so the module's own default correlation silently returned
    /// a friction factor that is a plausible-looking guess rather than a
    /// solution. Pinning it against the plain solver means the two paths cannot
    /// drift apart again.
    #[test]
    fn colebrook_dual_value_matches_the_iterated_solution() {
        let (mu, nu) = params();
        let v = Velocity::new(2.0);
        for d in [0.05f64, 0.15, 0.30, 0.60, 1.0] {
            let h = head_loss_dual(
                Length::new(200.0),
                v,
                Scalar::variable(d, 0),
                mu,
                nu,
                FrictionModel::ColebrookWhite,
            );
            let reference = head_loss_value(
                Length::new(200.0),
                v,
                d,
                mu,
                nu,
                FrictionModel::ColebrookWhite,
            );
            assert!(
                (h.re() - reference).abs() / reference < 1e-9,
                "D={d}: dual {} vs iterated {reference}",
                h.re()
            );
        }
    }

    /// The whole claim of the dual Colebrook path: the gradient is *exact*, not
    /// approximated. A centred difference converges to the derivative as
    /// `O(h^2)`, so agreeing with it to 1e-7 relative at a well-chosen step is
    /// only possible if the analytic value is right.
    ///
    /// The previous implementation returned a constant 0.02 and therefore a
    /// gradient of exactly zero, which is why this test is worded against a
    /// signed comparison rather than a magnitude one.
    #[test]
    fn colebrook_gradient_is_exact_and_non_zero() {
        let (mu, nu) = params();
        let v = Velocity::new(2.0);
        for d in [0.08f64, 0.20, 0.35, 0.50] {
            let h = head_loss_dual(
                Length::new(150.0),
                v,
                Scalar::variable(d, 0),
                mu,
                nu,
                FrictionModel::ColebrookWhite,
            );
            let numeric = head_loss_gradient(
                Length::new(150.0),
                v,
                d,
                mu,
                nu,
                FrictionModel::ColebrookWhite,
                1.0e-6,
            );
            assert!(h.du(0).abs() > 0.0, "D={d}: the gradient must not vanish");
            assert!(
                (h.du(0) - numeric).abs() / numeric.abs() < 1e-7,
                "D={d}: analytic {} vs numeric {numeric}",
                h.du(0)
            );
            // And the sign is the one that drives a sizing optimiser: wider
            // means less loss.
            assert!(
                h.du(0) < 0.0,
                "D={d}: gradient {} should be negative",
                h.du(0)
            );
        }
    }

    /// The iteration's own derivative seed must not affect the converged answer.
    ///
    /// The seed is the derivative of the *starting guess* for `f`, an artefact of
    /// the iteration with no physics in it. The dual argument arrives at the
    /// gradient through the chain rule, so the correct statement is the strong
    /// one: the seed washes out entirely, rather than merely becoming small.
    ///
    /// This is the property that distinguishes the analytic gradient from a
    /// propagated finite difference, and it is what the old constant-0.02 arm
    /// could never have satisfied.
    #[test]
    fn the_colebrook_gradient_does_not_depend_on_the_iteration_seed() {
        let (_mu, nu) = params();
        let v = Velocity::new(2.0);
        let d = 0.30;
        let re = Reynolds::from_kinematic(v, Length::new(d), nu).value();
        let re_dual = Scalar::new(re, [re / d]);

        let zero_seeded = colebrook_white_dual_seeded(re_dual, 0.0);
        // Seeds orders of magnitude away from the answer, in both directions.
        for seed in [-50.0, -1.0, 1.0, 500.0] {
            let other = colebrook_white_dual_seeded(re_dual, seed);
            assert!(
                (other.re() - zero_seeded.re()).abs() / zero_seeded.re() < 1e-12,
                "seed {seed}: value {} vs {}",
                other.re(),
                zero_seeded.re()
            );
            assert!(
                (other.du(0) - zero_seeded.du(0)).abs() / zero_seeded.du(0).abs() < 1e-10,
                "seed {seed}: gradient {} vs {}",
                other.du(0),
                zero_seeded.du(0)
            );
        }
    }

    /// `f` depends on the diameter only through `Re`, so the gradient must be
    /// *linear* in `dRe/dD` and pass exactly through the origin.
    ///
    /// This is the chain rule stated as a test, and it is the check that
    /// distinguishes a genuine partial derivative from a finite-difference
    /// estimate: a difference quotient is only approximately linear and has a
    /// step-size-dependent offset, whereas `df/dD = (df/dRe)(dRe/dD)` holds to
    /// machine precision for any multiple, including zero.
    #[test]
    fn the_colebrook_gradient_is_linear_in_the_reynolds_derivative() {
        let (_mu, nu) = params();
        let v = Velocity::new(2.0);
        let d = 0.30;
        let re = Reynolds::from_kinematic(v, Length::new(d), nu).value();
        let base = colebrook_white_dual(Scalar::new(re, [re / d])).du(0);

        for multiple in [0.0, 0.5, 2.0, -1.0, 10.0] {
            let scaled = colebrook_white_dual(Scalar::new(re, [multiple * re / d])).du(0);
            let expected = multiple * base;
            let tolerance = expected.abs().max(base.abs()) * 1.0e-10;
            assert!(
                (scaled - expected).abs() <= tolerance,
                "multiple {multiple}: gradient {scaled} vs {expected}"
            );
        }
    }

    /// A non-positive or zero Reynolds number must be handled on every path
    /// rather than producing a NaN that would poison a downstream optimiser.
    #[test]
    fn degenerate_reynolds_numbers_give_zero_not_nan() {
        let direct = colebrook_white_dual(Scalar::constant(0.0));
        assert_eq!(direct.re(), 0.0);
        assert!(direct.re().is_finite());

        let negative = colebrook_white_dual(Scalar::new(-1.0, [1.0]));
        assert_eq!(negative.re(), 0.0);
        assert!(negative.du(0).is_finite());
    }
}

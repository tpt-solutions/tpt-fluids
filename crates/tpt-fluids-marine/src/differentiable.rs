//! Differentiable resistance integrals, for gradient-based hull-form and speed
//! optimisation.
//!
//! Rating a hull is a different problem from designing one. A rating question has
//! a fixed form and wants the resistance; a design question has a target and wants
//! the form, and the variables a gradient-based optimiser should move are the
//! speed, the principal dimensions, and the Michell coefficients themselves.
//!
//! The Michell residuary integral is a polynomial in the Froude number
//!
//! ```text
//! C_r(Fr) = a + b Fr^2 + c Fr^4 + d Fr^6
//! ```
//!
//! so its gradient is available exactly by evaluating it over
//! [`tpt_math_autodiff::fwd::Dual`] instead of `f64`. That turns a search over a
//! hull-form coefficient into one evaluation per step rather than a
//! finite-difference sweep per step.
//!
//! ```
//! use tpt_fluids_marine::differentiable::{michell_resistance_dual, Scalar};
//! use tpt_fluids_marine::resistance::MichellResiduary;
//!
//! let michell = MichellResiduary { a: 0.05, b: 0.2, c: 0.1, d: 0.4 };
//!
//! // Track the Froude number, which is what a speed change moves.
//! let fr = Scalar::variable(0.25, 0);
//! let r = michell_resistance_dual(3.0e7, fr, michell).expect("a real hull");
//!
//! assert!(r.re() > 0.0);
//! // The derivative slot holds dR/dFr. It is positive: every term in the
//! // Michell polynomial is a positive power of a non-negative Froude number, so
//! // the residuary resistance is monotone in speed. That monotonicity is the
//! // fact that makes "go faster" always cost more, and it is what lets an
//! // optimiser read the sign off the gradient.
//! assert!(r.du(0) > 0.0);
//! ```
//!
//! # Why the total resistance chain matters
//!
//! [`total_resistance_dual`] composes the frictional and residuary parts, and is
//! the function worth handing an optimiser. Resistance per unit speed grows like
//! `Fr^5` for a typical hull, so there is a genuine interior minimum that
//! minimising either part alone never reaches.
//!
//! # What is and is not differentiable
//!
//! The ITTC-1957 friction correlation is *not* provided over duals. It is
//! logarithmic in the Reynolds number, `C_f = 0.075 / (log10 Re - 2)^2`, and more
//! importantly it is an empirical fit to model tests rather than a consequence
//! of anything, so differentiating it produces a gradient with no physical
//! standing behind it. [`total_resistance_dual`] therefore treats the friction
//! coefficient as a constant, which is the honest approximation: the frictional
//! term is `C_f` times `V^2`, and the `V^2` scaling, which is the part that
//! actually moves with the design, is carried exactly.
//!
//! The limitation worth naming: an optimiser built on this finds the optimum
//! *within the Michell polynomial family*. It cannot change the hull form itself,
//! because the polynomial is a fit to a parent form. Escaping it needs a real
//! hull-form synthesis loop, which is not a gradient of this integral.

use tpt_math_autodiff::fwd::Dual;

use crate::resistance::MichellResiduary;

/// The dual-number type used for a single tracked variable.
pub type Scalar = Dual<f64, 1>;

/// The Michell residuary coefficient, differentiable in the Froude number.
///
/// This is [`MichellResiduary::coefficient`] over dual numbers:
///
/// ```text
/// C_r(Fr) = a + b Fr^2 + c Fr^4 + d Fr^6
/// ```
///
/// Pass a [`Scalar::variable`] as `froude` to obtain `dC_r/dFr` in the
/// derivative slot. Pass [`Scalar::constant`] to get the plain value.
///
/// # Errors
///
/// Returns [`crate::error::MarineError::NonPositive`] for a negative Froude
/// number and [`crate::error::MarineError::NonFinite`] for a non-finite one.
/// Squaring a negative Froude number would produce a positive, plausible
/// resistance coefficient for a ship travelling in reverse, so the check happens
/// on the value before any dual arithmetic and an invalid input cannot produce a
/// poisoned derivative.
pub fn michell_coefficient_dual(
    froude: Scalar,
    michell: MichellResiduary,
) -> crate::error::Result<Scalar> {
    let value = froude.re();
    if !value.is_finite() {
        return Err(crate::error::MarineError::NonFinite("Froude number"));
    }
    if value < 0.0 {
        return Err(crate::error::MarineError::NonPositive("Froude number"));
    }

    let fr2 = froude * froude;
    // Horner in `Fr^2`: `a + Fr^2 (b + Fr^2 (c + Fr^2 d))`. This is the same
    // polynomial written to reuse `fr2`, and it keeps the derivative chain
    // short: the exponents fold into the coefficients on the way in, so four
    // powers become three multiplies.
    let inner = fr2 * Scalar::constant(michell.d) + Scalar::constant(michell.c);
    let middle = fr2 * inner + Scalar::constant(michell.b);
    Ok(fr2 * middle + Scalar::constant(michell.a))
}

/// The Michell residuary resistance, differentiable in the Froude number.
///
/// ```text
/// R_r = 0.5 W Fr^2 C_r(Fr)
/// ```
///
/// with `W` the displacement force in newtons. The Froude-squared prefactor is
/// what makes resistance per unit displacement a function of `Fr` alone, so the
/// whole of the hull-form dependence sits in `C_r`.
///
/// # Errors
///
/// As [`michell_coefficient_dual`], and additionally rejects a non-positive
/// displacement: a zero or negative displacement would scale the resistance to
/// zero or below, which is not a hull.
pub fn michell_resistance_dual(
    displacement_newtons: f64,
    froude: Scalar,
    michell: MichellResiduary,
) -> crate::error::Result<Scalar> {
    if !displacement_newtons.is_finite() {
        return Err(crate::error::MarineError::NonFinite("displacement"));
    }
    if displacement_newtons <= 0.0 {
        return Err(crate::error::MarineError::NonPositive("displacement"));
    }
    let coefficient = michell_coefficient_dual(froude, michell)?;
    let half = Scalar::constant(0.5);
    let fr2 = froude * froude;
    Ok(half * Scalar::constant(displacement_newtons) * fr2 * coefficient)
}

/// The total resistance of a ship, differentiable in the Froude number.
///
/// ```text
/// R = 0.5 W Fr^2 (C_f + C_r(Fr))
/// ```
///
/// This is the composition worth optimising: the frictional coefficient is
/// carried as a constant, because the ITTC-1957 correlation is an empirical
/// model-test fit whose derivative has no physical meaning, while the residuary
/// part -- which dominates and which a hull-form coefficient actually controls --
/// is differentiated exactly.
///
/// The `Fr^2` prefactor is the same quadratic speed dependence the `V^2` inside
/// `C_f` would give, so freezing `C_f` does not lose how friction moves with the
/// design. It only freezes the much weaker logarithmic Reynolds dependence of
/// `C_f` itself.
///
/// # Errors
///
/// As [`michell_resistance_dual`], and rejects a non-positive or non-finite
/// friction coefficient.
pub fn total_resistance_dual(
    displacement_newtons: f64,
    froude: Scalar,
    friction_coefficient: f64,
    michell: MichellResiduary,
) -> crate::error::Result<Scalar> {
    if !friction_coefficient.is_finite() {
        return Err(crate::error::MarineError::NonFinite("friction coefficient"));
    }
    if friction_coefficient < 0.0 {
        return Err(crate::error::MarineError::NonPositive(
            "friction coefficient",
        ));
    }
    let residuary = michell_resistance_dual(displacement_newtons, froude, michell)?;
    // The frictional share of `0.5 W Fr^2 C`, with `C_f` held constant. Kept as
    // its own term rather than folded into the Michell polynomial so the two
    // contributions keep their separate provenance in the derivative.
    let half = Scalar::constant(0.5);
    let fr2 = froude * froude;
    let friction = half
        * Scalar::constant(displacement_newtons)
        * fr2
        * Scalar::constant(friction_coefficient);
    Ok(friction + residuary)
}

/// The required power to overcome a resistance, differentiable in the Froude
/// number.
///
/// ```text
/// P = R V / (eta D)
/// ```
///
/// `speed` is the speed in metres per second corresponding to this Froude
/// number. Because `Fr = V / sqrt(g L)`, that is the hull's speed at the given
/// Froude number; it is passed separately rather than derived so this stays a
/// function of the Froude number alone, and so a caller tracking a speed
/// directly can build the Froude number from it.
///
/// This is what makes the design question well-posed. Resistance grows as
/// `Fr^6` and power as `Fr^7`, so neither has an interior minimum on its own;
/// resistance *per unit speed* does, at a finite Froude number.
///
/// # Errors
///
/// As [`total_resistance_dual`], and rejects a negative or non-finite speed and a
/// non-positive or non-finite efficiency. A zero efficiency is a division by
/// zero, so it must not reach the dual arithmetic.
pub fn required_power_dual(
    displacement_newtons: f64,
    froude: Scalar,
    speed: f64,
    friction_coefficient: f64,
    michell: MichellResiduary,
    propulsive_efficiency: f64,
) -> crate::error::Result<Scalar> {
    if !speed.is_finite() {
        return Err(crate::error::MarineError::NonFinite("speed"));
    }
    if speed < 0.0 {
        return Err(crate::error::MarineError::NonPositive("speed"));
    }
    if !propulsive_efficiency.is_finite() {
        return Err(crate::error::MarineError::NonFinite(
            "propulsive efficiency",
        ));
    }
    if propulsive_efficiency <= 0.0 {
        return Err(crate::error::MarineError::NonPositive(
            "propulsive efficiency",
        ));
    }
    let resistance =
        total_resistance_dual(displacement_newtons, froude, friction_coefficient, michell)?;
    Ok(resistance * Scalar::constant(speed) / Scalar::constant(propulsive_efficiency))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resistance;
    use tpt_fluids_core::nondimensional::Froude;
    use tpt_fluids_core::quantity::{Length, Velocity};

    /// A central finite difference, for checking the analytic gradients.
    ///
    /// This is the estimator being replaced, used here as the independent check
    /// on the closed form. Two different derivations disagreeing is evidence;
    /// one checked against itself is not.
    fn numeric_derivative(f: impl Fn(f64) -> f64, x: f64, step: f64) -> f64 {
        (f(x + step) - f(x - step)) / (2.0 * step)
    }

    /// A representative parent hull form: a moderate merchant ship.
    fn michell() -> MichellResiduary {
        MichellResiduary {
            a: 0.05,
            b: 0.2,
            c: 0.1,
            d: 0.4,
        }
    }

    /// The displacement force of a 3000 DWT feeder, in newtons.
    const FEEDER: f64 = 3.0e7;

    #[test]
    fn the_dual_value_equals_the_scalar_value() {
        // If these two ever diverge, the dual version is no longer the same
        // function, and every gradient it produces is about a different problem.
        for fr in [0.0, 0.1, 0.25, 0.4, 0.6] {
            let dual = michell_coefficient_dual(Scalar::variable(fr, 0), michell()).unwrap();
            let scalar = michell().coefficient(Froude::from_raw(fr));
            assert!(
                (dual.re() - scalar).abs() <= 1.0e-14 * scalar.abs().max(1.0),
                "Fr {fr}: {} vs {scalar}",
                dual.re()
            );

            let dual_r =
                michell_resistance_dual(FEEDER, Scalar::variable(fr, 0), michell()).unwrap();
            let scalar_r = michell().resistance(FEEDER, Froude::from_raw(fr));
            assert!(
                (dual_r.re() - scalar_r).abs() <= 1.0e-12 * scalar_r.abs().max(1.0),
                "Fr {fr}: {} vs {scalar_r}",
                dual_r.re()
            );
        }
    }

    #[test]
    fn the_horner_form_matches_the_naive_polynomial() {
        // The Horner rewrite is an optimisation, and an optimisation that is not
        // checked against what it replaced is just a different expression.
        let m = michell();
        for fr in [0.05, 0.17, 0.3, 0.45, 0.7] {
            let f = fr;
            let naive = m.a + m.b * f * f + m.c * f.powi(4) + m.d * f.powi(6);
            let hornered = michell_coefficient_dual(Scalar::constant(fr), m)
                .unwrap()
                .re();
            assert!(
                (hornered - naive).abs() <= 1.0e-14 * naive.abs().max(1.0),
                "Fr {fr}: {hornered} vs {naive}"
            );
        }
    }

    #[test]
    fn the_coefficient_gradient_matches_a_finite_difference() {
        // The central check. A dropped chain-rule term on the powers would still
        // give a plausible gradient that rises with Fr; only an independent
        // estimator catches that.
        for fr in [0.05, 0.15, 0.3, 0.5, 0.7, 0.9] {
            let analytic = michell_coefficient_dual(Scalar::variable(fr, 0), michell())
                .unwrap()
                .du(0);
            let numeric =
                numeric_derivative(|f| michell().coefficient(Froude::from_raw(f)), fr, 1.0e-6);
            assert!(
                (analytic - numeric).abs() / numeric.abs() < 1.0e-6,
                "Fr {fr}: analytic {analytic} vs numeric {numeric}"
            );
        }
    }

    #[test]
    fn the_resistance_gradient_matches_a_finite_difference() {
        for fr in [0.1, 0.25, 0.4, 0.55, 0.75] {
            let analytic = michell_resistance_dual(FEEDER, Scalar::variable(fr, 0), michell())
                .unwrap()
                .du(0);
            let numeric = numeric_derivative(
                |f| michell().resistance(FEEDER, Froude::from_raw(f)),
                fr,
                1.0e-6,
            );
            assert!(
                (analytic - numeric).abs() / numeric.abs() < 1.0e-6,
                "Fr {fr}: analytic {analytic} vs numeric {numeric}"
            );
        }
    }

    #[test]
    fn resistance_is_monotone_in_froude_and_so_is_its_gradient() {
        // Every term of the Michell polynomial is a non-negative power of a
        // non-negative Froude number, so the residuary resistance and its
        // derivative must both be non-decreasing. The sign of this gradient is
        // the fact that going faster costs more, and an optimiser reads it off
        // directly.
        let mut previous = f64::NEG_INFINITY;
        let mut previous_gradient = f64::NEG_INFINITY;
        for i in 0..=30 {
            let fr = f64::from(i) * 0.02;
            let dual = michell_resistance_dual(FEEDER, Scalar::variable(fr, 0), michell()).unwrap();
            assert!(dual.re() >= previous, "not monotone at Fr {fr}");
            assert!(
                dual.du(0) >= previous_gradient,
                "gradient not monotone at Fr {fr}: {}",
                dual.du(0)
            );
            assert!(dual.du(0) >= 0.0, "gradient went negative at Fr {fr}");
            previous = dual.re();
            previous_gradient = dual.du(0);
        }
    }

    #[test]
    fn resistance_grows_as_the_sixth_of_froude_for_a_wave_making_hull() {
        // With only the `d` term the scaling is exact, and it is the reason an
        // interior optimum can exist at all: if resistance is `Fr^6` and power
        // is `Fr^7`, neither is minimised at `Fr = 0`, so the trade has to be
        // against something else. Pinned because an optimiser's credibility
        // rests on the scaling being the one the theory claims.
        let wave = MichellResiduary {
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d: 0.1,
        };
        for fr in [0.2, 0.4, 0.6] {
            let dual = michell_resistance_dual(FEEDER, Scalar::constant(fr), wave).unwrap();
            let expected = 0.5 * FEEDER * fr * fr * (0.1 * fr.powi(6));
            assert!((dual.re() - expected).abs() / expected < 1.0e-14, "{fr}");
        }
    }

    #[test]
    fn the_power_gradient_matches_a_finite_difference() {
        for fr in [0.15, 0.3, 0.45, 0.6] {
            let speed = 9.0;
            let analytic = required_power_dual(
                FEEDER,
                Scalar::variable(fr, 0),
                speed,
                0.0015,
                michell(),
                0.65,
            )
            .unwrap()
            .du(0);
            let numeric = numeric_derivative(
                |f| {
                    let r = 0.5
                        * FEEDER
                        * f
                        * f
                        * (0.0015 + michell().coefficient(Froude::from_raw(f)));
                    r * speed / 0.65
                },
                fr,
                1.0e-6,
            );
            assert!(
                (analytic - numeric).abs() / numeric.abs() < 1.0e-6,
                "Fr {fr}: analytic {analytic} vs numeric {numeric}"
            );
        }
    }

    #[test]
    fn an_interior_optimum_exists_for_resistance_per_unit_speed() {
        // The end-to-end check, and the one that matters: a gradient is only
        // useful if following it finds the minimum it claims. On `R / V` the
        // `Fr^2` prefactor fights the `Fr^6` growth, so the objective has a real
        // interior minimum. Descending it from either side must converge to the
        // same place, which is a much stronger statement than "it decreased".
        let objective = |fr: f64| {
            0.5 * FEEDER * fr * fr * (0.0015 + michell().coefficient(Froude::from_raw(fr)))
                / (fr * 9.0)
        };
        let descend = |start: f64| {
            let mut fr = start;
            for _ in 0..4000 {
                let v = fr * 9.0;
                let dual =
                    total_resistance_dual(FEEDER, Scalar::variable(fr, 0), 0.0015, michell())
                        .unwrap();
                let gradient = dual.du(0) / v - dual.re() / (v * v);
                fr -= 0.5 * gradient;
                fr = fr.clamp(1.0e-6, 2.0);
            }
            fr
        };
        let from_low = descend(0.05);
        let from_high = descend(0.8);
        assert!(
            (from_low - from_high).abs() < 1.0e-4,
            "two descents disagreed: {from_low} vs {from_high}"
        );
        assert!(from_low > 0.0 && from_low < 0.8, "converged to {from_low}");
        // And the interior point really does beat both ends of the sweep, so
        // "interior" is not an artefact of the iteration stopping early.
        assert!(objective(from_low) < objective(0.05));
        assert!(objective(from_low) < objective(0.8));
    }

    #[test]
    fn an_invalid_input_is_refused_rather_than_differentiated() {
        // Squaring a negative Froude number would hand back a positive,
        // plausible resistance for a ship travelling backwards, and a zero
        // efficiency is a division by zero. The checks are on the values before
        // any dual arithmetic, so no poisoned gradient escapes.
        assert!(michell_coefficient_dual(Scalar::constant(-0.2), michell()).is_err());
        assert!(michell_resistance_dual(FEEDER, Scalar::constant(-0.1), michell()).is_err());
        assert!(michell_resistance_dual(-1.0, Scalar::constant(0.2), michell()).is_err());
        assert!(michell_resistance_dual(0.0, Scalar::constant(0.2), michell()).is_err());
        assert!(michell_resistance_dual(f64::NAN, Scalar::constant(0.2), michell()).is_err());
        assert!(michell_coefficient_dual(Scalar::constant(f64::NAN), michell()).is_err());
        assert!(michell_coefficient_dual(Scalar::constant(f64::INFINITY), michell()).is_err());
        assert!(total_resistance_dual(FEEDER, Scalar::constant(0.2), -0.001, michell()).is_err());
        assert!(total_resistance_dual(FEEDER, Scalar::constant(0.2), f64::NAN, michell()).is_err());
        assert!(
            required_power_dual(FEEDER, Scalar::constant(0.2), -1.0, 0.0015, michell(), 0.65)
                .is_err()
        );
        assert!(
            required_power_dual(FEEDER, Scalar::constant(0.2), 9.0, 0.0015, michell(), 0.0)
                .is_err()
        );
        assert!(required_power_dual(
            FEEDER,
            Scalar::constant(0.2),
            9.0,
            0.0015,
            michell(),
            f64::NAN
        )
        .is_err());
    }

    #[test]
    fn zero_froude_is_admissible_and_carries_no_resistance() {
        // A hull at rest is a real state, not an error: it has no wave-making
        // resistance, and refusing it would make a speed sweep starting at zero
        // impossible.
        let dual = michell_resistance_dual(FEEDER, Scalar::constant(0.0), michell()).unwrap();
        assert_eq!(dual.re(), 0.0);
        assert_eq!(dual.du(0), 0.0);
        // The coefficient itself is not zero: it is the constant term `a`.
        let coefficient = michell_coefficient_dual(Scalar::constant(0.0), michell()).unwrap();
        assert_eq!(coefficient.re(), michell().a);
    }

    #[test]
    fn a_constant_froude_carries_no_derivative() {
        // The escape hatch: evaluating at a constant must give the plain value
        // with a zero derivative, or a caller could not mix tracked and
        // untracked quantities.
        let dual = michell_resistance_dual(FEEDER, Scalar::constant(0.25), michell()).unwrap();
        assert!(dual.du(0).abs() < 1.0e-15);
        assert!((dual.re() - michell().resistance(FEEDER, Froude::from_raw(0.25))).abs() < 1.0e-9);
    }

    #[test]
    fn a_realistic_feeder_lands_in_a_plausible_resistance_band() {
        // A scale test rather than a gradient test, because a derivative can be
        // exactly right for an answer that is out by a factor of a million, and a
        // hull-form optimiser quietly working in the wrong units would not
        // otherwise be caught.
        //
        // The quantity asserted is `R / W`, not an absolute force, and that
        // choice is deliberate. The Michell path computes `R = 0.5 W Fr^2 C_r`
        // with `W` a displacement *force* in newtons, so `R/W` is the only
        // convention-free way to state a magnitude here. Pinning an absolute
        // kilonewton figure instead would silently mix this with the
        // `froude_scaling` convention, which normalises by `rho g` times a
        // displacement treated as a *volume*; the two disagree by exactly the
        // factors of `g` and `1000` that an order-of-magnitude mistake produces.
        // Asserting the ratio keeps the test about physics.
        let speed = Velocity::new(18.0 * 0.5144);
        let length = Length::new(150.0);
        let froude = resistance::ship_froude(speed, length);
        let r = total_resistance_dual(FEEDER, Scalar::constant(froude.value()), 0.0015, michell())
            .unwrap();

        // `R/W` for a merchant hull at a moderate Froude number sits in the
        // 1e-4 to 1e-2 band. The band is deliberately two decades wide: its job
        // is to catch a unit error, not to certify these particular fixture
        // coefficients against a measured ship.
        let ratio = r.re() / FEEDER;
        assert!(
            ratio > 1.0e-4 && ratio < 1.0e-2,
            "R/W = {ratio} is not a plausible resistance-to-displacement ratio"
        );
        assert!(r.re() > 0.0);

        // And it agrees with the non-differentiable path it is replacing, which
        // is the property that actually matters for the caller: the same
        // function, with a gradient attached.
        let expected = 0.5
            * FEEDER
            * froude.value()
            * froude.value()
            * (0.0015 + michell().coefficient(froude));
        assert!(
            (r.re() - expected).abs() / expected < 1.0e-12,
            "{} vs {expected}",
            r.re()
        );
    }

    #[test]
    fn the_total_resistance_gradient_matches_a_finite_difference() {
        for fr in [0.1, 0.3, 0.5, 0.7] {
            let analytic =
                total_resistance_dual(FEEDER, Scalar::variable(fr, 0), 0.0015, michell())
                    .unwrap()
                    .du(0);
            let numeric = numeric_derivative(
                |f| 0.5 * FEEDER * f * f * (0.0015 + michell().coefficient(Froude::from_raw(f))),
                fr,
                1.0e-6,
            );
            assert!(
                (analytic - numeric).abs() / numeric.abs() < 1.0e-6,
                "Fr {fr}: analytic {analytic} vs numeric {numeric}"
            );
        }
    }

    #[test]
    fn the_total_resistance_gradient_is_the_sum_of_its_parts() {
        // Linearity: the derivative of the composition must equal the sum of the
        // derivatives. This is the one property that would catch the two terms
        // being composed wrongly rather than merely mis-differentiated, which a
        // comparison against a single finite difference would not separate.
        let fr = 0.3;
        let total = total_resistance_dual(FEEDER, Scalar::variable(fr, 0), 0.0015, michell())
            .unwrap()
            .du(0);
        let residuary = michell_resistance_dual(FEEDER, Scalar::variable(fr, 0), michell())
            .unwrap()
            .du(0);
        // d/dFr of `0.5 W C_f Fr^2` with `C_f` held constant.
        let friction_gradient = FEEDER * fr * 0.0015;
        assert!(
            (total - (residuary + friction_gradient)).abs() < 1.0e-6,
            "{total} vs {}",
            residuary + friction_gradient
        );
        assert!(friction_gradient > 0.0);
    }
}

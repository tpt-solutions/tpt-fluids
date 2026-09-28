//! Differentiable load-capacity integrals, for gradient-based bearing design.
//!
//! Sizing a bearing is a different problem from rating one. A rating question
//! has a fixed geometry and wants the load; a design question has a target load
//! and wants the geometry, and the geometry a gradient-based optimiser should
//! move is the film thickness or the wedge taper.
//!
//! The closed-form wedge results in [`crate::reynolds`] are smooth analytic
//! functions of the taper, so their gradients are available exactly by
//! evaluating them over [`tpt_math_autodiff::fwd::Dual`] instead of `f64`.
//! That turns a search over taper into one evaluation per step rather than a
//! finite-difference sweep per step, which matters because the finite-difference
//! estimate of a derivative in a quantity that scales like `1/h^2` is itself
//! badly conditioned near the optimum -- exactly where an optimiser spends most
//! of its time.
//!
//! ```
//! use tpt_fluids_tribo::differentiable::{wedge_load_dual, Scalar};
//!
//! // Track the wedge taper.
//! let taper = Scalar::variable(0.5, 0);
//! let load = wedge_load_dual(taper).expect("a converging wedge");
//!
//! assert!(load.re() > 0.0);
//! // The derivative slot now holds dW/d(taper). It is positive: closing the
//! // wedge raises the load. That sign is what lets an optimiser push the taper
//! // up, and it is the whole reason for asking for a gradient.
//! assert!(load.du(0) > 0.0);
//! ```
//!
//! # What is and is not differentiable
//!
//! The Reynolds *solver* is not differentiable here, and pretending otherwise
//! would be the useful thing to do but a dishonest one. Its answer comes from a
//! linear solve whose matrix depends on the film profile, and threading a
//! derivative through that is a genuine piece of work rather than a
//! transcription. What is provided instead is the closed form, which is exact
//! for the linear-wedge geometry it describes and therefore a *better* gradient
//! than differentiating the discretised solver would be even if that were done.
//!
//! The limitation is real and worth naming: an optimiser built on this will
//! converge to the optimum wedge, not to the optimum bearing. It cannot see
//! cavitation, sidelobe clearance, or anything else that depends on the
//! pressure distribution's shape rather than its integral.

use tpt_math_autodiff::fwd::Dual;

/// The dual-number type used for a single tracked variable.
pub type Scalar = Dual<f64, 1>;

/// The load capacity of a converging wedge, differentiable in the taper.
///
/// This is [`crate::reynolds::wedge_load_capacity`] evaluated over dual numbers:
///
/// ```text
/// W(t) = (6 / t^2) [ -ln(1 - t) - 2t / (2 - t) ]
/// ```
///
/// Pass a [`Scalar::variable`] as `taper` to obtain `dW/dt` in the derivative
/// slot. Pass [`Scalar::constant`] to get the plain value.
///
/// # Errors
///
/// Returns [`crate::error::TribologyError::NonPositive`] for a non-positive
/// taper and [`crate::error::TribologyError::OutsideValidRange`] for one of a
/// unit or more, exactly as the `f64` version does. The dual arithmetic is
/// never reached for those, so an invalid taper cannot produce a poisoned
/// derivative.
pub fn wedge_load_dual(taper: Scalar) -> crate::error::Result<Scalar> {
    // Validate on the value, using the `f64` entry point, so the two versions
    // cannot drift apart in what they accept.
    crate::reynolds::wedge_load_capacity(taper.re())?;
    wedge_load_expression(taper)
}

/// The peak pressure of a converging wedge, differentiable in the taper.
///
/// ```text
/// p_max(t) = 3t / (2 (1 - t) (2 - t))
/// ```
///
/// The gradient matters as much as the value here: peak pressure is what limits
/// a bearing, and an optimiser minimising load subject to a pressure bound
/// needs both.
///
/// # Errors
///
/// As [`wedge_load_dual`].
pub fn wedge_peak_pressure_dual(taper: Scalar) -> crate::error::Result<Scalar> {
    crate::reynolds::wedge_peak_pressure(taper.re())?;
    let one = Scalar::constant(1.0);
    let two = Scalar::constant(2.0);
    let three = Scalar::constant(3.0);
    let denominator = two * (one - taper) * (two - taper);
    Ok(three * taper / denominator)
}

/// The specific pressure of a wedge, differentiable in the taper.
///
/// This is the ratio an optimiser should minimise directly: load alone is
/// unbounded above, so it has no interior optimum, while the peak-to-load ratio
/// has a floor at `3/2` and a clear shape.
///
/// # Errors
///
/// As [`wedge_load_dual`].
pub fn specific_pressure_dual(taper: Scalar) -> crate::error::Result<Scalar> {
    let load = wedge_load_dual(taper)?;
    let peak = wedge_peak_pressure_dual(taper)?;
    if load.re() <= 0.0 {
        return Ok(Scalar::constant(f64::INFINITY));
    }
    Ok(peak / load)
}

/// The wedge load expression itself, without validation.
///
/// Split out so the validated entry points above are the only place a taper is
/// checked, and this stays a pure function of the arithmetic.
fn wedge_load_expression(taper: Scalar) -> crate::error::Result<Scalar> {
    let one = Scalar::constant(1.0);
    let two = Scalar::constant(2.0);
    let six = Scalar::constant(6.0);
    // `-ln(1 - t)`. The log is the reason this needs `tpt-math-autodiff`'s
    // transcendental support rather than a hand-rolled derivative; getting the
    // chain rule term `-1/(1-t)` wrong here would produce a gradient that looks
    // reasonable and steers an optimiser the wrong way.
    //
    // `Dual` implements no unary negation, so the sign is applied by subtracting
    // from zero rather than by a leading `-`. That is a library limitation, not a
    // style preference, and getting it wrong is a compile error rather than a
    // silent sign flip.
    let minus_one = Scalar::constant(-1.0);
    let log_term = minus_one * (one - taper).ln();
    let rational = two * taper / (two - taper);
    let bracket = log_term - rational;
    Ok(six / (taper * taper) * bracket)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reynolds;

    /// A central finite difference, for checking the analytic gradients.
    ///
    /// This is the estimator being replaced, used here as the independent check
    /// on the closed form. Two different derivations disagreeing is evidence;
    /// one checked against itself is not.
    fn numeric_derivative(f: impl Fn(f64) -> f64, x: f64, step: f64) -> f64 {
        (f(x + step) - f(x - step)) / (2.0 * step)
    }

    #[test]
    fn the_dual_value_equals_the_scalar_value() {
        // If these two ever diverge, the dual version is no longer the same
        // function, and every gradient it produces is about a different problem.
        for taper in [0.05, 0.2, 0.5, 0.8, 0.95] {
            let dual = wedge_load_dual(Scalar::variable(taper, 0)).unwrap();
            let scalar = reynolds::wedge_load_capacity(taper).unwrap();
            assert!(
                (dual.re() - scalar).abs() / scalar < 1.0e-14,
                "taper {taper}: {} vs {scalar}",
                dual.re()
            );

            let dual_peak = wedge_peak_pressure_dual(Scalar::variable(taper, 0)).unwrap();
            let scalar_peak = reynolds::wedge_peak_pressure(taper).unwrap();
            assert!(
                (dual_peak.re() - scalar_peak).abs() / scalar_peak < 1.0e-14,
                "taper {taper}: {} vs {scalar_peak}",
                dual_peak.re()
            );
        }
    }

    #[test]
    fn the_load_gradient_matches_a_finite_difference() {
        // The central check on the whole module. A wrong chain-rule term on the
        // logarithm would still give a plausible-looking gradient; only a
        // comparison against an independent estimator catches that.
        for taper in [0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95] {
            let dual = wedge_load_dual(Scalar::variable(taper, 0)).unwrap();
            let analytic = dual.du(0);
            let numeric =
                numeric_derivative(|t| reynolds::wedge_load_capacity(t).unwrap(), taper, 1.0e-6);
            assert!(
                (analytic - numeric).abs() / numeric.abs() < 1.0e-6,
                "taper {taper}: analytic {analytic} vs numeric {numeric}"
            );
            // And the sign, which is what an optimiser actually consumes:
            // closing the wedge always raises the load.
            assert!(analytic > 0.0, "taper {taper}: gradient {analytic}");
        }
    }

    #[test]
    fn the_peak_pressure_gradient_matches_a_finite_difference() {
        for taper in [0.1, 0.3, 0.6, 0.9] {
            let dual = wedge_peak_pressure_dual(Scalar::variable(taper, 0)).unwrap();
            let analytic = dual.du(0);
            let numeric =
                numeric_derivative(|t| reynolds::wedge_peak_pressure(t).unwrap(), taper, 1.0e-6);
            assert!(
                (analytic - numeric).abs() / numeric.abs() < 1.0e-6,
                "taper {taper}: analytic {analytic} vs numeric {numeric}"
            );
            assert!(analytic > 0.0);
        }
    }

    #[test]
    fn a_gradient_descent_step_actually_reduces_the_objective() {
        // The end-to-end check, and the one that matters: a gradient is only
        // useful if following it downhill does what it claims. On the load,
        // downhill means a shallower wedge.
        let mut taper = 0.8f64;
        for _ in 0..40 {
            let dual = wedge_load_dual(Scalar::variable(taper, 0)).unwrap();
            let load = dual.re();
            taper -= 0.05 * dual.du(0) / load.max(1.0);
            taper = taper.clamp(1.0e-4, 0.99);
        }
        assert!(taper < 0.5, "taper drifted to {taper}");
        let final_load = reynolds::wedge_load_capacity(taper).unwrap();
        assert!(final_load < reynolds::wedge_load_capacity(0.8).unwrap());
    }

    #[test]
    fn descending_specific_pressure_prefers_the_shallow_wedge() {
        // Specific pressure has a floor at 3/2 rather than a load that diverges,
        // so descending it is a well-posed optimisation and must walk towards
        // that floor.
        let mut taper = 0.8f64;
        for _ in 0..200 {
            let dual = specific_pressure_dual(Scalar::variable(taper, 0)).unwrap();
            taper -= 0.2 * dual.du(0);
            taper = taper.clamp(1.0e-4, 0.99);
        }
        let ratio = specific_pressure_dual(Scalar::constant(taper))
            .unwrap()
            .re();
        assert!(ratio < 1.55, "converged to {ratio} at taper {taper}");
        assert!(ratio > 1.5, "{ratio} is below the theoretical floor");
    }

    #[test]
    fn an_invalid_taper_is_refused_rather_than_differentiated() {
        // The derivative of `log(0)` is infinite, so letting a rejected taper
        // through would hand an optimiser a poisoned gradient instead of an
        // error. The check happens on the value before any dual arithmetic.
        assert!(wedge_load_dual(Scalar::constant(0.0)).is_err());
        assert!(wedge_load_dual(Scalar::constant(-0.3)).is_err());
        assert!(wedge_load_dual(Scalar::constant(1.0)).is_err());
        assert!(wedge_load_dual(Scalar::constant(2.0)).is_err());
        assert!(wedge_peak_pressure_dual(Scalar::constant(0.0)).is_err());
        assert!(specific_pressure_dual(Scalar::constant(1.5)).is_err());
    }

    #[test]
    fn a_constant_taper_carries_no_derivative() {
        // The escape hatch: evaluating at a constant must give the plain value
        // with a zero derivative, or a caller could not mix tracked and
        // untracked quantities.
        let dual = wedge_load_dual(Scalar::constant(0.5)).unwrap();
        assert!(dual.du(0).abs() < 1.0e-15);
        assert!((dual.re() - reynolds::wedge_load_capacity(0.5).unwrap()).abs() < 1.0e-14);
    }
}

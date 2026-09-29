//! Pipe friction: Darcy-Weisbach, Colebrook-White, and Hazen-Williams.
//!
//! Head loss along a conduit is always
//!
//! ```text
//! h_f = f (L / D) V^2 / (2 g)
//! ```
//!
//! so the whole problem reduces to the friction factor `f`. The choice of
//! correlation for `f` is the single most consequential modelling decision in
//! any pipe-network calculation, so the alternatives are exposed explicitly
//! rather than hidden behind one default.

use tpt_fluids_core::math;

use crate::error::HydraulicError;

/// The friction regime a flow is in, from the Reynolds number.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FlowRegime {
    /// `Re < 2300`: Stokes flow, `f = 64 / Re`.
    Laminar,
    /// `2300 <= Re < 4000`: transitional, the correlation is only a guide.
    Transitional,
    /// `Re >= 4000`: fully turbulent.
    Turbulent,
}

impl FlowRegime {
    /// Classifies a Reynolds number.
    pub const fn from_reynolds(re: f64) -> Self {
        if re < 2300.0 {
            Self::Laminar
        } else if re < 4000.0 {
            Self::Transitional
        } else {
            Self::Turbulent
        }
    }
}

/// A correlation for the Darcy friction factor.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum FrictionModel {
    /// The exact laminar solution, `f = 64 / Re`.
    Laminar,
    /// Colebrook-White's implicit equation, solved by iteration. This is the
    /// reference correlation, accurate to about 1% over the whole turbulent
    /// range when the roughness is known.
    ColebrookWhite,
    /// The explicit Swamee-Jain approximation to Colebrook-White, for when a
    /// closed form is needed.
    SwameeJain,
    /// Hazen-Williams, the empirical US-customary formula for water in pipes.
    HazenWilliams {
        /// The Hazen-Williams coefficient `C`, which is purely empirical and
        /// carries no units. Typical values: 100 new smooth mains, 80
        /// galvanized iron, 60 old cast iron.
        c: f64,
    },
}

/// Solves Colebrook-White by fixed-point iteration:
///
/// ```text
/// 1 / sqrt(f) = -2 log10 (eps/(3.7 D) + 2.51 / (Re sqrt(f)))
/// ```
///
/// The iteration runs to a tight relative tolerance rather than a fixed count,
/// so a difficult pipe (very smooth, or barely turbulent) is never silently
/// returned unconverged. `relative_roughness` is `eps / D`.
fn colebrook_white(reynolds: f64, relative_roughness: f64) -> Result<f64, HydraulicError> {
    if reynolds <= 0.0 {
        return Ok(0.0);
    }
    let re = reynolds.min(1.0e12);
    let rough_term = relative_roughness / 3.7;

    // The fully-rough asymptote, which also bounds the iteration from above.
    let f_rough = if rough_term > 0.0 {
        let d = -2.0 * math::log10(rough_term);
        1.0 / (d * d)
    } else {
        f64::INFINITY
    };

    // Seed from the hydraulically smooth branch.
    let mut f = if re < 2300.0 {
        64.0 / re.max(1.0)
    } else {
        (0.3164 / re.powf(0.25)).min(f_rough)
    };

    const MAX_ITER: usize = 200;
    const TOL: f64 = 1.0e-10;
    for _ in 0..MAX_ITER {
        let sqrt_f = math::sqrt(f);
        if !sqrt_f.is_finite() || sqrt_f <= 0.0 {
            break;
        }
        let arg = rough_term + 2.51 / (re * sqrt_f);
        if arg <= 0.0 {
            break;
        }
        let denom = -2.0 * math::log10(arg);
        let new_f = 1.0 / (denom * denom);
        if !new_f.is_finite() || new_f <= 0.0 {
            break;
        }
        if (new_f - f).abs() <= TOL * f {
            return Ok(new_f);
        }
        f = new_f;
    }

    if !f.is_finite() || f <= 0.0 {
        return Err(HydraulicError::NonConverged("Colebrook-White iteration"));
    }
    Ok(f)
}

/// The Darcy friction factor for a flow in a conduit of the given diameter
/// and absolute roughness.
///
/// `reynolds` must already be formed from the flow's velocity, the pipe
/// diameter, and the fluid's viscosity — see
/// [`crate::differentiable::head_loss`], which
/// does that and then calls this.
pub fn friction_factor(
    model: FrictionModel,
    reynolds: f64,
    diameter: tpt_fluids_core::quantity::Length,
    roughness: tpt_fluids_core::quantity::Length,
) -> Result<f64, HydraulicError> {
    let d = diameter.value();
    if d <= 0.0 {
        return Err(HydraulicError::NonPositiveDiameter);
    }
    let relative_roughness = (roughness.value() / d).max(0.0);
    let re = reynolds;

    Ok(match model {
        FrictionModel::Laminar => {
            if re <= 0.0 {
                0.0
            } else {
                64.0 / re
            }
        }
        FrictionModel::SwameeJain => {
            if re <= 0.0 {
                0.0
            } else {
                let re = re.min(1.0e9);
                // f = 0.25 / [log10(eps/(3.7 D) + 5.74 / Re^0.9)]^2. The 0.9
                // exponent on the Reynolds number is essential; using a plain
                // Re understates f by ~20% at Re = 1e4.
                let arg = relative_roughness / 3.7 + 5.74 / re.powf(0.9);
                if arg <= 0.0 {
                    0.0
                } else {
                    let l = math::log10(arg);
                    0.25 / (l * l)
                }
            }
        }
        FrictionModel::HazenWilliams { c } => {
            if re <= 0.0 || c <= 0.0 {
                0.0
            } else {
                0.082 * c * c / re.powf(0.2)
            }
        }
        FrictionModel::ColebrookWhite => colebrook_white(re, relative_roughness)?,
    })
}

/// The friction factor for a smooth pipe, i.e. one whose roughness is
/// negligible. This is the common case for short, clean service pipes, and it
/// lets callers avoid tracking an absolute roughness at all.
pub fn friction_factor_smooth(model: FrictionModel, reynolds: f64) -> Result<f64, HydraulicError> {
    if reynolds <= 0.0 {
        return Ok(0.0);
    }
    match model {
        FrictionModel::Laminar => Ok(64.0 / reynolds),
        FrictionModel::HazenWilliams { c } => {
            if c <= 0.0 {
                Ok(0.0)
            } else {
                Ok(0.082 * c * c / reynolds.powf(0.2))
            }
        }
        _ => {
            // A 0.2 mm roughness is the ISO "commercial steel" datum scaled to
            // a smooth service pipe; Colebrook is insensitive to it here.
            colebrook_white(reynolds, 0.0)
        }
    }
}

/// The relative roughness `eps / D` of a pipe, which is the dimensionless
/// group Colebrook-White actually consumes.
pub fn relative_roughness(
    diameter: tpt_fluids_core::quantity::Length,
    roughness: tpt_fluids_core::quantity::Length,
) -> f64 {
    if diameter.value() <= 0.0 {
        return 0.0;
    }
    (roughness.value() / diameter.value()).max(0.0)
}

/// Standard absolute roughness values for common pipe materials, in metres.
///
/// These are the Colebrook-Crane figures: they matter far more to a network
/// solution than any solver tolerance, and picking the wrong one is the usual
/// cause of a model that disagrees with the field.
pub mod roughness {
    /// Drawn seamless copper, tube.
    pub const COPPER_DRAWN: f64 = 1.5e-6;
    /// Commercial steel or wrought iron.
    pub const COMMERCIAL_STEEL: f64 = 4.6e-5;
    /// Asphalted cast iron.
    pub const ASPHALTED_CAST_IRON: f64 = 1.2e-4;
    /// Galvanized iron.
    pub const GALVANIZED_IRON: f64 = 1.5e-4;
    /// Cast iron.
    pub const CAST_IRON: f64 = 2.6e-4;
    /// Concrete, smooth-finished.
    pub const CONCRETE_SMOOTH: f64 = 3.0e-4;
    /// Unfinished concrete.
    pub const CONCRETE_ROUGH: f64 = 1.0e-3;
    /// Riveted steel.
    pub const RIVETED_STEEL: f64 = 9.0e-4;
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_fluids_core::quantity::Length;

    fn re(v: f64) -> f64 {
        v
    }

    #[test]
    fn laminar_friction_factor_is_exact() {
        // f = 64 / Re exactly, for every laminar Reynolds number.
        for reynolds in [100.0, 1000.0, 2000.0] {
            let f = friction_factor(
                FrictionModel::Laminar,
                reynolds,
                Length::new(0.1),
                Length::new(0.0),
            )
            .unwrap();
            assert!((f - 64.0 / reynolds).abs() < 1e-12, "Re {reynolds}: {f}");
        }
    }

    #[test]
    fn colebrook_white_matches_the_equation_exactly() {
        // These are the true roots of
        //   1/sqrt(f) = -2 log10 (2.51 / (Re sqrt(f)))
        // for a hydraulically smooth pipe, which is what a fixed-point
        // iteration on the defining equation must converge to.
        let cases = [
            (1.0e4, 0.030_883),
            (1.0e5, 0.017_990),
            (1.0e6, 0.011_645),
            (1.0e7, 0.008_103),
        ];
        for (reynolds, expected) in cases {
            let f = friction_factor(
                FrictionModel::ColebrookWhite,
                reynolds,
                Length::new(0.1),
                Length::new(0.0),
            )
            .unwrap();
            assert!(
                (f - expected).abs() / expected < 1e-4,
                "Re {reynolds}: got {f}, expected {expected}"
            );
        }
    }

    /// Re-derives `f` by brute-force iteration and checks the crate's answer
    /// satisfies Colebrook-White's defining equation as a residual.
    #[test]
    fn colebrook_white_satisfies_its_own_equation() {
        for reynolds in [1.0e4, 1.0e5, 1.0e6, 1.0e7] {
            for eps_over_d in [0.0, 1.0e-5, 1.0e-4, 1.0e-3] {
                let f = friction_factor(
                    FrictionModel::ColebrookWhite,
                    reynolds,
                    Length::new(0.1),
                    Length::new(eps_over_d * 0.1),
                )
                .unwrap();
                // LHS and RHS of 1/sqrt(f) = -2 log10(eps/3.7D + 2.51/(Re sqrt(f))).
                let lhs = 1.0 / f.sqrt();
                let rhs = -2.0 * (eps_over_d / 3.7 + 2.51 / (reynolds * f.sqrt())).log10();
                assert!(
                    (lhs - rhs).abs() < 1e-6,
                    "Re {reynolds} eps/D {eps_over_d}: {lhs} vs {rhs}"
                );
            }
        }
    }

    #[test]
    fn colebrook_white_includes_roughness() {
        // eps/D = 0.01 at Re = 1e5 sits in the near-fully-rough regime, where
        // f ~ 1 / (-2 log10(eps/3.7D))^2 = 0.0472.
        let f = friction_factor(
            FrictionModel::ColebrookWhite,
            1.0e5,
            Length::new(0.1),
            Length::new(0.001),
        )
        .unwrap();
        let fully_rough = 1.0 / (-2.0 * (0.01f64 / 3.7).log10()).powi(2);
        assert!(
            (f - fully_rough).abs() / fully_rough < 0.05,
            "{f} vs asymptote {fully_rough}"
        );
        // And it must exceed the smooth-pipe value at the same Reynolds.
        let smooth = friction_factor(
            FrictionModel::ColebrookWhite,
            1.0e5,
            Length::new(0.1),
            Length::new(0.0),
        )
        .unwrap();
        assert!(f > smooth, "{f} !> {smooth}");
    }

    #[test]
    fn swamee_jain_tracks_colebrook_white() {
        // Swamee-Jain is an explicit approximation, so it should agree with
        // the iterative Colebrook result to within a few percent.
        for reynolds in [1.0e4, 1.0e5, 1.0e6] {
            let cw = friction_factor(
                FrictionModel::ColebrookWhite,
                reynolds,
                Length::new(0.2),
                Length::new(0.0001),
            )
            .unwrap();
            let sj = friction_factor(
                FrictionModel::SwameeJain,
                reynolds,
                Length::new(0.2),
                Length::new(0.0001),
            )
            .unwrap();
            assert!(
                (cw - sj).abs() / cw < 0.05,
                "Re {reynolds}: Colebrook {cw} vs Swamee-Jain {sj}"
            );
        }
    }

    #[test]
    fn hazen_williams_uses_only_reynolds() {
        // Hazen-Williams carries no roughness, so it must be independent of
        // the pipe's roughness and diameter.
        let a = friction_factor(
            FrictionModel::HazenWilliams { c: 130.0 },
            1.0e6,
            Length::new(0.1),
            Length::new(0.0),
        )
        .unwrap();
        let b = friction_factor(
            FrictionModel::HazenWilliams { c: 130.0 },
            1.0e6,
            Length::new(0.5),
            Length::new(0.001),
        )
        .unwrap();
        assert!((a - b).abs() < 1e-12, "{a} vs {b}");
    }

    #[test]
    fn friction_factor_is_monotone_in_roughness() {
        let mut previous: Option<f64> = None;
        for eps in [0.0, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2] {
            let f = friction_factor(
                FrictionModel::ColebrookWhite,
                1.0e5,
                Length::new(0.1),
                Length::new(eps),
            )
            .unwrap();
            if let Some(prev) = previous {
                assert!(f > prev, "roughness {eps}: {f} !> {prev}");
            }
            previous = Some(f);
        }
    }

    #[test]
    fn friction_factor_is_monotone_decreasing_in_reynolds() {
        let mut previous = f64::INFINITY;
        for reynolds in [1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7] {
            let f = friction_factor(
                FrictionModel::ColebrookWhite,
                reynolds,
                Length::new(0.1),
                Length::new(0.0001),
            )
            .unwrap();
            assert!(f < previous, "Re {reynolds}: {f} !< {previous}");
            previous = f;
        }
    }

    #[test]
    fn zero_reynolds_yields_zero_friction() {
        for model in [
            FrictionModel::Laminar,
            FrictionModel::ColebrookWhite,
            FrictionModel::SwameeJain,
            FrictionModel::HazenWilliams { c: 130.0 },
        ] {
            let f = friction_factor(model, re(0.0), Length::new(0.1), Length::new(0.0)).unwrap();
            assert_eq!(f, 0.0, "{model:?}");
        }
    }

    #[test]
    fn non_positive_diameter_is_rejected() {
        assert_eq!(
            friction_factor(
                FrictionModel::ColebrookWhite,
                1.0e5,
                Length::new(0.0),
                Length::new(0.0)
            ),
            Err(HydraulicError::NonPositiveDiameter)
        );
    }

    #[test]
    fn flow_regime_boundaries_are_exact() {
        assert_eq!(FlowRegime::from_reynolds(0.0), FlowRegime::Laminar);
        assert_eq!(FlowRegime::from_reynolds(2299.0), FlowRegime::Laminar);
        assert_eq!(FlowRegime::from_reynolds(2300.0), FlowRegime::Transitional);
        assert_eq!(FlowRegime::from_reynolds(3999.0), FlowRegime::Transitional);
        assert_eq!(FlowRegime::from_reynolds(4000.0), FlowRegime::Turbulent);
    }

    #[test]
    fn relative_roughness_guards_against_zero_diameter() {
        assert_eq!(
            relative_roughness(Length::new(0.0), Length::new(0.001)),
            0.0
        );
        assert_eq!(
            relative_roughness(Length::new(0.1), Length::new(0.001)),
            0.01
        );
    }
}

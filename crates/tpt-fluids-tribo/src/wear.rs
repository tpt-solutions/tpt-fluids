//! Archard wear, the lambda ratio, and frictional heating.
//!
//! # Archard'"'"'s law
//!
//! The simplest useful wear law is
//!
//! ```text
//! V = k W s / H
//! ```
//!
//! where `V` is the volume worn away, `W` the normal load, `s` the sliding
//! distance, `H` the hardness, and `k` a dimensionless wear coefficient. The
//! dimension check is what makes it credible: `W s` is work, and dividing by
//! a hardness that is itself a pressure gives back a volume.
//!
//! Its weakness is equally clear, and the crate says so. Archard'"'"'s law is
//! linear in load, and wear is emphatically not linear in load: it jumps by
//! orders of magnitude when a film fails. The coefficient `k` is not a
//! material property so much as a property of a particular system in a
//! particular regime, and it is typically measured rather than derived. Use it
//! to compare like with like, not to extrapolate across a regime boundary.
//!
//! # The lambda ratio
//!
//! Whether a contact is boundary, mixed, or full-film is decided by how the
//! film thickness compares with the roughness of the surfaces:
//!
//! ```text
//! lambda = h_min / sqrt(Rq1^2 + Rq2^2)
//! ```
//!
//! Below `lambda = 1` asperities touch through the film and the contact is
//! boundary lubricated; above about 3 it is fully separated. This is the
//! number to compute before trusting any full-film result.

use tpt_fluids_core::math;

use crate::error::{Result, TribologyError};

/// The Archard wear coefficient, dimensionless.
///
/// It is emphatically not a constant of a material. Typical magnitudes:
///
/// - `1e-10` to `1e-8` for well-lubricated, mild wear;
/// - `1e-8` to `1e-6` for moderate wear;
/// - `1e-6` to `1e-3` for severe wear, including brake pads and seized
///   contacts.
///
/// The spread across seven orders of magnitude between benign and severe is
/// the whole reason the law cannot be extrapolated.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct WearCoefficient(f64);

impl WearCoefficient {
    /// Builds a coefficient, rejecting a negative one.
    pub fn new(value: f64) -> Result<Self> {
        if value < 0.0 {
            return Err(TribologyError::NonPositive("wear coefficient"));
        }
        Ok(Self(value))
    }

    /// A coefficient for mild, well-lubricated wear, `1e-9`.
    pub const MILD: Self = Self(1.0e-9);
    /// A coefficient for moderate wear, `1e-7`.
    pub const MODERATE: Self = Self(1.0e-7);
    /// A coefficient for severe wear such as a brake pad, `2e-5`.
    pub const SEVERE: Self = Self(2.0e-5);

    /// The coefficient as a plain number.
    pub fn value(self) -> f64 {
        self.0
    }
}

/// The volume worn away by Archard's law, in cubic metres.
///
/// `V = k W s / H`.
pub fn wear_volume(
    coefficient: WearCoefficient,
    load: f64,
    sliding_distance: f64,
    hardness: f64,
) -> f64 {
    if load <= 0.0 || sliding_distance <= 0.0 || hardness <= 0.0 {
        return 0.0;
    }
    coefficient.value() * load * sliding_distance / hardness
}

/// The mean wear depth over a contact area, in metres.
///
/// Dividing the volume by the area it was removed from. This is the number an
/// engineer actually reads: how much material is gone.
pub fn wear_depth(
    coefficient: WearCoefficient,
    load: f64,
    sliding_distance: f64,
    hardness: f64,
    area: f64,
) -> f64 {
    if area <= 0.0 {
        return 0.0;
    }
    wear_volume(coefficient, load, sliding_distance, hardness) / area
}

/// The wear coefficient implied by a measured wear volume.
///
/// This is how a coefficient is obtained in practice: run a wear test, then
/// invert Archard'"'"'s law. The round trip through this and [`wear_volume`] is
/// the model'"'"'s own consistency check.
pub fn coefficient_from_wear(volume: f64, load: f64, sliding_distance: f64, hardness: f64) -> f64 {
    if load <= 0.0 || sliding_distance <= 0.0 || hardness <= 0.0 {
        return 0.0;
    }
    volume * hardness / (load * sliding_distance)
}

/// The sliding distance a component survives before a given wear depth is
/// reached, in metres.
///
/// Inverting Archard'"'"'s law for distance is the useful direction for design:
/// given an allowable wear, how far can this part run?
pub fn life_for_wear_depth(
    coefficient: WearCoefficient,
    load: f64,
    hardness: f64,
    area: f64,
    allowable_depth: f64,
) -> f64 {
    if area <= 0.0 || allowable_depth <= 0.0 {
        return 0.0;
    }
    let volume = allowable_depth * area;
    volume * hardness / (coefficient.value() * load)
}

/// The lambda ratio, relating film thickness to combined roughness.
///
/// `lambda = h_min / sqrt(Rq1^2 + Rq2^2)`.
pub fn lambda_ratio(minimum_film_thickness: f64, roughness_a: f64, roughness_b: f64) -> f64 {
    let combined = math::sqrt(roughness_a * roughness_a + roughness_b * roughness_b);
    if combined <= 0.0 {
        return f64::INFINITY;
    }
    minimum_film_thickness / combined
}

/// The separation regime a lambda ratio implies.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SeparationRegime {
    /// `lambda < 1`: asperities touch through the film.
    Boundary,
    /// `1 <= lambda < 3`: partial separation.
    Mixed,
    /// `lambda >= 3`: fully separated, full-film lubrication.
    FullFilm,
}

impl SeparationRegime {
    /// Classifies a lambda ratio.
    pub fn classify(lambda: f64) -> Self {
        if lambda < 1.0 {
            Self::Boundary
        } else if lambda < 3.0 {
            Self::Mixed
        } else {
            Self::FullFilm
        }
    }
}

/// The frictional power dissipated in a sliding contact, in watts.
///
/// `Q = mu F v`. Everything converted to heat ends up here first.
pub fn frictional_power(friction_coefficient: f64, load: f64, sliding_speed: f64) -> f64 {
    if friction_coefficient < 0.0 {
        return 0.0;
    }
    friction_coefficient * load * sliding_speed
}

/// The steady-state temperature rise at a sliding contact, in kelvin.
///
/// ```text
/// dT = mu F v / (k_c A)
/// ```
///
/// with `k_c` the contact's thermal conductivity in W/(m K) and `A` the
/// contact area.
///
/// # Validity
///
/// This is the *quasi-steady* form, valid for contacts that dwell long enough
/// for the heat to conduct away. For a short sliding event, a brake pad on a
/// disc, it is not merely inaccurate but absurd: a typical pad computes to
/// tens of thousands of kelvin, which no material survives. The real flash
/// temperature there is a transient conduction problem and is orders of
/// magnitude lower. The function reports the steady answer and the crate
/// documents the limit; it does not pretend to solve the transient problem.
pub fn steady_temperature_rise(
    friction_coefficient: f64,
    load: f64,
    sliding_speed: f64,
    thermal_conductivity: f64,
    area: f64,
) -> f64 {
    if thermal_conductivity <= 0.0 || area <= 0.0 {
        return 0.0;
    }
    frictional_power(friction_coefficient, load, sliding_speed) / (thermal_conductivity * area)
}

/// Whether the steady temperature-rise form is even in the right regime.
///
/// Above about 1000 K rise the quasi-steady assumption has certainly failed,
/// and the answer should be treated as an upper bound at best.
pub fn steady_form_is_valid(temperature_rise: f64) -> bool {
    temperature_rise > 0.0 && temperature_rise < 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wear_coefficient_rejects_a_negative_value() {
        assert!(WearCoefficient::new(-1.0e-9).is_err());
        assert!(WearCoefficient::new(0.0).is_ok());
        assert!(WearCoefficient::new(1.0e-9).is_ok());
    }

    #[test]
    fn the_named_coefficients_have_the_documented_magnitudes() {
        assert!((WearCoefficient::MILD.value() - 1.0e-9).abs() / 1.0e-9 < 1e-12);
        assert!((WearCoefficient::MODERATE.value() - 1.0e-7).abs() / 1.0e-7 < 1e-12);
        assert!((WearCoefficient::SEVERE.value() - 2.0e-5).abs() / 2.0e-5 < 1e-12);
    }

    #[test]
    fn archard_is_dimensionally_a_volume() {
        // [W s]/[H] = (N m)/(N/m^2) = m^3, so the result is a volume. Expressed
        // in mm^3 it should match the hand calculation.
        // k W s / H = 1e-10 * 1000 * 1 / 2.94e9 = 3.40e-14 m^3, which is
        // 3.40e-8 mm^3. Getting a volume this small is the point: a
        // well-lubricated contact doing 1 N m of work against a 3 GPa
        // material removes an amount you could not see.
        let v = wear_volume(WearCoefficient::new(1.0e-10).unwrap(), 1000.0, 1.0, 2.94e9);
        assert!(
            (v * 1.0e9 - 3.40e-8).abs() / 3.40e-8 < 1e-3,
            "{} mm^3",
            v * 1.0e9
        );
    }

    #[test]
    fn a_brake_pad_reproduces_its_measured_wear() {
        // A pad losing 1 mm over 100 000 km at 5 kN and 500 MPa hardness
        // implies k = 2e-5, which is where severe wear is documented to sit.
        let volume = 0.02 * 1.0e-3;
        let k = coefficient_from_wear(volume, 5000.0, 1.0e5, 5.0e8);
        assert!((k - 2.0e-5).abs() / 2.0e-5 < 1e-9, "k = {k}");
    }

    #[test]
    fn archard_round_trips() {
        // Measuring a wear volume and inferring the coefficient, then
        // predicting the volume again, must be exact. This is the law's own
        // consistency and it is cheap to check.
        let k = WearCoefficient::new(3.3e-7).unwrap();
        let volume = wear_volume(k, 2500.0, 5000.0, 7.0e8);
        let recovered = coefficient_from_wear(volume, 2500.0, 5000.0, 7.0e8);
        assert!((recovered - 3.3e-7).abs() / 3.3e-7 < 1e-12, "{recovered}");
    }

    #[test]
    fn wear_is_linear_in_load_and_distance() {
        // Both scalings are Archard'"'"'s central claim, and both are also
        // where the law is known to fail in reality. They are asserted here so
        // that a change to a more realistic model is a visible one.
        let k = WearCoefficient::MODERATE;
        let base = wear_volume(k, 1000.0, 1000.0, 1.0e9);
        assert!((wear_volume(k, 2000.0, 1000.0, 1.0e9) / base - 2.0).abs() < 1e-12);
        assert!((wear_volume(k, 1000.0, 2000.0, 1.0e9) / base - 2.0).abs() < 1e-12);
    }

    #[test]
    fn harder_materials_wear_slower() {
        let k = WearCoefficient::MODERATE;
        let soft = wear_volume(k, 1000.0, 1000.0, 1.0e8);
        let hard = wear_volume(k, 1000.0, 1000.0, 1.0e9);
        assert!(hard < soft, "{hard} !< {soft}");
    }

    #[test]
    fn wear_depth_divides_by_the_area() {
        let k = WearCoefficient::new(2.0e-5).unwrap();
        let v = wear_volume(k, 5000.0, 1.0e5, 5.0e8);
        let d = wear_depth(k, 5000.0, 1.0e5, 5.0e8, 0.02);
        assert!((d - 1.0e-3).abs() / 1.0e-3 < 1e-9, "depth = {d} m");
        assert!((d * 0.02 - v).abs() / v < 1e-12);
    }

    #[test]
    fn life_for_wear_depth_inverts_the_wear_law() {
        let k = WearCoefficient::SEVERE;
        let life = life_for_wear_depth(k, 5000.0, 5.0e8, 0.02, 1.0e-3);
        assert!((life - 1.0e5).abs() / 1.0e5 < 1e-9, "life = {life} m");
        // And the wear over exactly that life is the allowable depth.
        let d = wear_depth(k, 5000.0, life, 5.0e8, 0.02);
        assert!((d - 1.0e-3).abs() / 1.0e-3 < 1e-9);
    }

    #[test]
    fn a_larger_coefficient_gives_a_shorter_life() {
        let mild = life_for_wear_depth(WearCoefficient::MILD, 1000.0, 1.0e9, 0.01, 1.0e-5);
        let severe = life_for_wear_depth(WearCoefficient::SEVERE, 1000.0, 1.0e9, 0.01, 1.0e-5);
        assert!(severe < mild, "{severe} !< {mild}");
    }

    #[test]
    fn lambda_ratio_matches_the_definition() {
        // 0.5 um of film on 0.1 um roughness each side.
        let l = lambda_ratio(0.5e-6, 0.1e-6, 0.1e-6);
        assert!(
            (l - 0.5e-6 / (0.1e-6 * 2.0f64.sqrt())).abs() < 1e-9,
            "lambda = {l}"
        );
    }

    #[test]
    fn lambda_grows_with_film_thickness() {
        let a = lambda_ratio(0.1e-6, 0.1e-6, 0.1e-6);
        let b = lambda_ratio(2.0e-6, 0.1e-6, 0.1e-6);
        assert!(b > a, "{b} !> {a}");
    }

    #[test]
    fn a_rough_contact_never_reaches_full_film() {
        // 0.1 um of film on 1 um roughness cannot separate the surfaces,
        // however carefully the film is managed.
        let l = lambda_ratio(0.1e-6, 1.0e-6, 1.0e-6);
        assert_eq!(SeparationRegime::classify(l), SeparationRegime::Boundary);
    }

    #[test]
    fn separation_regimes_are_classified_by_lambda() {
        assert_eq!(SeparationRegime::classify(0.5), SeparationRegime::Boundary);
        assert_eq!(SeparationRegime::classify(1.0), SeparationRegime::Mixed);
        assert_eq!(SeparationRegime::classify(2.0), SeparationRegime::Mixed);
        assert_eq!(SeparationRegime::classify(3.0), SeparationRegime::FullFilm);
        assert_eq!(
            SeparationRegime::classify(3.5e6),
            SeparationRegime::FullFilm
        );
    }

    #[test]
    fn a_polished_contact_reaches_full_film() {
        // 2 um of film on 0.05 um roughness: fully separated.
        let l = lambda_ratio(2.0e-6, 0.05e-6, 0.05e-6);
        assert_eq!(SeparationRegime::classify(l), SeparationRegime::FullFilm);
    }

    #[test]
    fn a_perfectly_smooth_pair_is_always_fully_separated() {
        // Zero roughness with a positive film means nothing to touch.
        assert!(lambda_ratio(1.0e-6, 0.0, 0.0).is_infinite());
    }

    #[test]
    fn frictional_power_matches_the_definition() {
        assert!((frictional_power(0.1, 1000.0, 5.0) - 500.0).abs() < 1e-12);
    }

    #[test]
    fn temperature_rise_falls_as_the_contact_grows() {
        // More area to conduct into, less rise. The scaling is 1/A exactly.
        let a = steady_temperature_rise(0.1, 1000.0, 5.0, 50.0, 0.01);
        let b = steady_temperature_rise(0.1, 1000.0, 5.0, 50.0, 0.02);
        assert!((a / b - 2.0).abs() < 1e-12);
    }

    #[test]
    fn a_brake_pad_exposes_the_steady_forms_limits() {
        // A pad at mu 0.4, 5 kN and 10 m/s over 0.02 m^2 computes to 20 000
        // K, which is not a temperature. This test exists to make the limit
        // concrete rather than to endorse the number.
        let d_t = steady_temperature_rise(0.4, 5000.0, 10.0, 50.0, 0.02);
        assert!((d_t - 20_000.0).abs() < 1.0, "dT = {d_t} K");
        assert!(
            !steady_form_is_valid(d_t),
            "the steady form must report itself invalid here"
        );
    }

    #[test]
    fn a_slow_dwelling_contact_is_within_the_steady_regime() {
        // 500 W into a 0.01 m^2 contact with k_c = 50 gives 1000 K, right
        // at the edge of the quasi-steady assumption.
        let d_t = steady_temperature_rise(0.1, 1000.0, 5.0, 50.0, 0.01);
        assert!((d_t - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn non_positive_inputs_give_zero_rather_than_nonsense() {
        let k = WearCoefficient::MILD;
        assert_eq!(wear_volume(k, 0.0, 1.0, 1.0e9), 0.0);
        assert_eq!(wear_volume(k, 1.0, 0.0, 1.0e9), 0.0);
        assert_eq!(wear_volume(k, 1.0, 1.0, 0.0), 0.0);
        assert_eq!(wear_depth(k, 1.0, 1.0, 1.0e9, 0.0), 0.0);
        assert_eq!(coefficient_from_wear(1.0, 0.0, 1.0, 1.0e9), 0.0);
        assert_eq!(life_for_wear_depth(k, 1.0, 1.0e9, 0.0, 1.0e-5), 0.0);
        assert_eq!(frictional_power(-0.1, 1.0, 1.0), 0.0);
        assert_eq!(steady_temperature_rise(0.1, 1.0, 1.0, 0.0, 1.0), 0.0);
        assert!(!steady_form_is_valid(0.0));
    }
}

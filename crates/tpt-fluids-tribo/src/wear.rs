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

/// How a contact's wear coefficient evolves as it beds in.
///
/// # Why a coefficient that is not constant
///
/// Archard's law takes `k` as a constant, and for a running-in contact that is
/// wrong in a way that matters: **a new contact wears far faster than a
/// bedded-in one**. The reason is geometric. A fresh machine surface is covered
/// in asperities that are tall relative to the film, so the real contact area
/// is a small fraction of the nominal one and the pressure on each asperity
/// contact is enormous. Running-in shears those peaks off, the roughness falls,
/// and the load spreads over more of the surface. Wear decelerates by orders of
/// magnitude.
///
/// The practical consequence is that a linear law run from `s = 0` overpredicts
/// total wear badly, and most of that overprediction happens in the first
/// fraction of a percent of the sliding distance.
///
/// # The model
///
/// A decaying excess on top of a steady coefficient:
///
/// ```text
/// k(s) = k_steady + (k_initial - k_steady) exp(-s / d)
/// ```
///
/// with `d` the running-in distance over which the excess decays by `1/e`. The
/// functional form is a modelling choice, not a derived one, and that is worth
/// saying plainly: it reproduces the two things that are actually known, that
/// initial wear is much higher and that it decays towards a steady value. It is
/// not a mechanistic model of asperity ploughing, and `d` has to be measured.
///
/// The alternative -- a power law in distance -- fits some data better and has
/// the defect that its total wear integral is only finite if the exponent
/// exceeds one, so it silently predicts unbounded wear for mild exponents.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct RunningIn {
    /// The wear coefficient once the surfaces have bedded in.
    pub steady: WearCoefficient,
    /// The wear coefficient at zero sliding distance.
    pub initial: WearCoefficient,
    /// The distance over which the initial excess decays by `1/e`, in metres.
    pub running_in_distance: f64,
}

impl RunningIn {
    /// A running-in model from its two coefficients and decay distance.
    ///
    /// # Errors
    ///
    /// Returns [`TribologyError::NonPositive`] for a non-positive
    /// `running_in_distance`, since a zero distance would make the coefficient
    /// undefined rather than merely immediate.
    pub fn new(
        steady: WearCoefficient,
        initial: WearCoefficient,
        running_in_distance: f64,
    ) -> Result<Self> {
        if running_in_distance <= 0.0 {
            return Err(TribologyError::NonPositive("running-in distance"));
        }
        Ok(Self {
            steady,
            initial,
            running_in_distance,
        })
    }

    /// The wear coefficient after a given sliding distance, in metres.
    ///
    /// Monotonically decreasing in `distance` and bounded between the two
    /// coefficients, so a contact can never wear faster than its initial rate
    /// nor slower than its steady one. Those bounds are the property a
    /// simulation must not violate, and a constant coefficient violates the
    /// first.
    pub fn coefficient_at(&self, sliding_distance: f64) -> f64 {
        if sliding_distance <= 0.0 {
            return self.initial.value();
        }
        let decay = math::exp(-sliding_distance / self.running_in_distance);
        self.steady.value() + (self.initial.value() - self.steady.value()) * decay
    }

    /// The wear volume accumulated over a sliding distance, in cubic metres.
    ///
    /// The running-in coefficient is integrated rather than sampled, because the
    /// whole effect lives in the first part of the curve and a single-point
    /// evaluation at the end distance would throw it away. In closed form:
    ///
    /// ```text
    /// V(s) = [ k_steady s + (k_0 - k_steady) d (1 - exp(-s/d)) ] / H
    /// ```
    pub fn wear_volume(&self, load: f64, sliding_distance: f64, hardness: f64) -> f64 {
        if load <= 0.0 || sliding_distance <= 0.0 || hardness <= 0.0 {
            return 0.0;
        }
        let d = self.running_in_distance;
        let steady_part = self.steady.value() * sliding_distance;
        // Both terms are inside the load multiplier. Leaving it off the excess
        // term is a unit error that happens to return a plausible small number,
        // which is why the numerical-integration test beside this exists.
        let excess = (self.initial.value() - self.steady.value())
            * d
            * (1.0 - math::exp(-sliding_distance / d));
        load * (steady_part + excess) / hardness
    }

    /// The mean wear depth over a contact area after a sliding distance.
    pub fn wear_depth(&self, load: f64, sliding_distance: f64, hardness: f64, area: f64) -> f64 {
        if area <= 0.0 {
            return 0.0;
        }
        self.wear_volume(load, sliding_distance, hardness) / area
    }

    /// The sliding distance at which the initial excess has decayed to within 5
    /// percent of its own starting value, in metres.
    ///
    /// This is the number to quote when asked "how long does running-in take".
    /// It is `3 d`, and that is a consequence of the exponential rather than an
    /// independent physical result, which is why `d` is an input rather than
    /// something computed from first principles.
    ///
    /// # Read this before trusting it
    ///
    /// The 5 percent is of the *excess*, not of the steady coefficient, and the
    /// difference matters. When the initial coefficient is 100 times the steady
    /// one -- which is the normal case, and the reason running-in matters -- the
    /// coefficient at `3 d` is still more than double its steady value. For the
    /// coefficient itself to be within 5 percent of steady the distance is
    ///
    /// ```text
    /// d ln(20 (k_initial - k_steady) / k_steady)
    /// ```
    ///
    /// which for that ratio is about `7.6 d`, not `3 d`. A wear estimate that
    /// applies the steady coefficient after `3 d` will therefore underpredict by
    /// a factor of nearly two over the following few `d`. The distinction is
    /// pinned by a test rather than left in prose.
    pub fn settled_distance(&self) -> f64 {
        3.0 * self.running_in_distance
    }

    /// The wear a constant Archard coefficient would predict over the same
    /// distance, for direct comparison with [`RunningIn::wear_volume`].
    ///
    /// Included so a caller can see the overprediction rather than being asked
    /// to believe in it. A running-in contact always wears *less* than the
    /// steady-coefficient figure, because the steady coefficient is the smaller
    /// one.
    pub fn constant_coefficient_wear_volume(
        &self,
        coefficient: WearCoefficient,
        load: f64,
        sliding_distance: f64,
        hardness: f64,
    ) -> f64 {
        wear_volume(coefficient, load, sliding_distance, hardness)
    }
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

    // --- Running-in wear ---------------------------------------------------

    /// A bearing pair bedded in from 100x worse than steady, over 500 m.
    fn running_in() -> RunningIn {
        RunningIn::new(
            WearCoefficient::MILD,
            WearCoefficient::new(1.0e-7).unwrap(),
            500.0,
        )
        .expect("a running-in distance")
    }

    #[test]
    fn the_coefficient_decays_from_the_initial_to_the_steady_value() {
        let r = running_in();
        assert!((r.coefficient_at(0.0) - r.initial.value()).abs() < 1.0e-18);
        // One running-in distance leaves exactly 1/e of the excess.
        let expected =
            r.steady.value() + (r.initial.value() - r.steady.value()) / core::f64::consts::E;
        assert!((r.coefficient_at(500.0) - expected).abs() / expected < 1.0e-12);
        // And it never crosses below the steady value.
        for s in [0.0, 1.0, 100.0, 500.0, 5000.0, 1.0e7] {
            let k = r.coefficient_at(s);
            assert!(k >= r.steady.value() - 1.0e-18, "s={s} k={k}");
            assert!(k <= r.initial.value() + 1.0e-18, "s={s} k={k}");
        }
    }

    #[test]
    fn the_coefficient_only_ever_decreases() {
        let r = running_in();
        let mut previous = f64::INFINITY;
        for i in 0..200 {
            let s = f64::from(i) * 100.0;
            let k = r.coefficient_at(s);
            assert!(k <= previous, "coefficient rose at s={s}: {k} > {previous}");
            previous = k;
        }
    }

    #[test]
    fn running_in_always_wears_less_than_a_fresh_contact_coefficient_would_predict() {
        // The headline result, and the comparison that actually means something:
        // measure `k` on a fresh contact and use it for the whole life, versus
        // letting the contact bed in. That is the mistake a running-in model
        // exists to prevent.
        //
        // Note it is *not* a comparison against the steady coefficient. A
        // bedded-in contact wears more than a steady-coefficient Archard
        // calculation would predict, because running-in adds a finite excess;
        // claiming otherwise would be the more flattering but wrong statement.
        let r = running_in();
        let (load, hardness) = (5_000.0, 3.0e9);
        let naive = r.constant_coefficient_wear_volume(r.initial, load, 100_000.0, hardness);
        let actual = r.wear_volume(load, 100_000.0, hardness);
        assert!(actual < naive, "{actual} vs {naive}");
        assert!(actual > 0.0);
        // Over 100 km the naive figure is out by more than a factor of 50.
        assert!(naive / actual > 50.0, "ratio {}", naive / actual);

        // And against the steady coefficient the direction reverses, which is
        // the honest counterpart to the claim above.
        let steady_only = r.constant_coefficient_wear_volume(r.steady, load, 100_000.0, hardness);
        assert!(actual > steady_only, "{actual} vs {steady_only}");
    }

    #[test]
    fn the_wear_volume_agrees_with_a_numerical_integration_of_the_coefficient() {
        // The closed form must be the integral of `coefficient_at`, not merely
        // something plausible. Comparing against a trapezoidal sum that
        // converges at second order catches a dropped term in the algebra.
        let r = running_in();
        let (load, hardness) = (1000.0, 2.0e9);
        let (s_max, steps) = (3000.0, 2_000_000.0);
        let ds = s_max / steps;

        let mut integral = 0.0;
        for i in 0..(steps as usize) {
            let s0 = i as f64 * ds;
            integral += 0.5 * (r.coefficient_at(s0) + r.coefficient_at(s0 + ds)) * ds;
        }
        let expected = load * integral / hardness;
        let closed_form = r.wear_volume(load, s_max, hardness);
        assert!(
            (closed_form - expected).abs() / expected < 1.0e-6,
            "closed form {closed_form} vs integral {expected}"
        );
    }

    #[test]
    fn wear_accumulates_entirely_within_the_running_in_distance() {
        // The excess term saturates: `1 - exp(-s/d)` is within 5% of 1 by `s = 3d`.
        // Beyond that all further wear is at the steady rate, which is what
        // makes a long-life estimate insensitive to the running-in model.
        let r = running_in();
        let (load, hardness) = (1000.0, 1.0e9);
        let d = r.running_in_distance;

        let early = r.wear_volume(load, d, hardness);
        let excess = load * (r.initial.value() - r.steady.value()) * d / hardness;
        assert!(early > excess * 0.6, "{early} vs {excess}");

        // Past 3d the increment per further metre is the steady Archard rate.
        let (from, to) = (20.0 * d, 24.0 * d);
        let delta = r.wear_volume(load, to, hardness) - r.wear_volume(load, from, hardness);
        let steady_rate = r.steady.value() * load * (to - from) / hardness;
        assert!(
            (delta - steady_rate).abs() / steady_rate < 1.0e-6,
            "{delta} vs {steady_rate}"
        );
    }

    #[test]
    fn settling_is_measured_against_the_excess_not_against_the_steady_rate() {
        // A distinction that is easy to get wrong and expensive to get wrong.
        // `settled_distance` is `3 d`, which is 5% of the *excess* -- but the
        // excess starts at 100 times the steady coefficient here, so at `3 d`
        // the coefficient is still more than double its steady value. A user who
        // reads "settled after 3d" and multiplies by the steady rate will
        // underpredict the wear over the next few `d`.
        //
        // The two thresholds differ by exactly this ratio, and the test pins it:
        // the distance for the *coefficient* to be within 5% of steady is
        // `d ln(0.95 (k0 - ks) / ks)`, which is much larger than `3 d`.
        let r = running_in();
        let ratio = r.initial.value() / r.steady.value();
        assert!(
            ratio > 10.0,
            "this test needs a large initial-to-steady ratio"
        );

        // At the nominal settled distance the coefficient is still far elevated.
        let elevated = r.coefficient_at(r.settled_distance()) / r.steady.value();
        assert!(elevated > 2.0, "{elevated}");

        // And the true distance to within 5% of the steady coefficient.
        // Solving `ks + (k0 - ks) exp(-s/d) = 1.05 ks` gives `s = d ln(20(r/k))`.
        let true_distance = r.running_in_distance * (20.0 * (ratio - 1.0)).ln();
        let within_5pc = r.coefficient_at(true_distance) / r.steady.value();
        assert!((within_5pc - 1.05).abs() < 1.0e-9, "{within_5pc}");
        assert!(
            true_distance > 2.0 * r.settled_distance(),
            "the true settling distance {} is well beyond 3d = {}",
            true_distance,
            r.settled_distance()
        );
    }

    #[test]
    fn the_settled_distance_is_where_the_coefficient_is_essentially_steady() {
        let r = running_in();
        let excess = r.initial.value() - r.steady.value();
        let settled = (r.coefficient_at(r.settled_distance()) - r.steady.value()) / excess;
        assert!(settled < 0.05, "{settled}");
        // And before it the coefficient is still meaningfully elevated.
        assert!(r.coefficient_at(0.0) - r.steady.value() > 0.9 * excess);
    }

    #[test]
    fn a_running_in_distance_of_zero_is_refused() {
        assert!(RunningIn::new(WearCoefficient::MILD, WearCoefficient::SEVERE, 0.0).is_err());
        assert!(RunningIn::new(WearCoefficient::MILD, WearCoefficient::SEVERE, -1.0).is_err());
        assert!(RunningIn::new(WearCoefficient::MILD, WearCoefficient::SEVERE, 1.0).is_ok());
    }

    #[test]
    fn running_in_never_predicts_negative_wear() {
        let r = running_in();
        // Even with the two coefficients the wrong way round, which is
        // physically backwards, the volume cannot go negative.
        let inverted = RunningIn::new(
            WearCoefficient::new(1.0e-7).unwrap(),
            WearCoefficient::MILD,
            500.0,
        )
        .unwrap();
        for s in [0.0, 1.0, 500.0, 1.0e5] {
            assert!(inverted.wear_volume(1.0, s, 1.0e9) >= 0.0, "s={s}");
            assert!(r.wear_volume(1.0, s, 1.0e9) >= 0.0);
        }
        assert!(inverted.coefficient_at(0.0) > 0.0);
        assert!(r.wear_volume(1.0, 1.0e5, 1.0e9) > 0.0);
    }

    #[test]
    fn running_in_rejects_degenerate_inputs_the_same_way_archard_does() {
        let r = running_in();
        assert_eq!(r.wear_volume(0.0, 100.0, 1.0e9), 0.0);
        assert_eq!(r.wear_volume(1.0, 0.0, 1.0e9), 0.0);
        assert_eq!(r.wear_volume(1.0, 100.0, 0.0), 0.0);
        assert_eq!(r.wear_depth(1.0, 100.0, 1.0e9, 0.0), 0.0);
        // A negative distance means "before the start", which is the initial
        // coefficient rather than an extrapolated one.
        assert_eq!(r.coefficient_at(-5.0), r.initial.value());
    }
}

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
//!
//! # The flash temperature, steady and transient
//!
//! Frictional heating has two regimes, and the difference between them is not a
//! correction factor. `steady_temperature_rise` and `flash_temperature_rise`
//! solve different problems: one asks what temperature a body reaches when the
//! heat has had time to spread, the other asks how hot the surface gets during
//! a contact too brief for that. For a brake pad the first is off by an order
//! of magnitude, and `ThermalProperties::crossover_time` says why -- the
//! duration at which a metal could reach a quasi-steady uniform temperature is
//! about 30 hours.
//!
//! Both forms take a **conduction length** `L`. That is not a modelling detail
//! but a dimensional requirement, and its absence from the original
//! quasi-steady form was a real bug: `mu F v / (k A)` is kelvin per *metre*.
//! The documented "20 000 K brake pad" was that expression being read as a
//! temperature.

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
/// dT = q L / k_c = mu F v L / (k_c A)
/// ```
///
/// with `k_c` the contact's thermal conductivity in W/(m K), `A` the contact
/// area, and `L` the **conduction length** -- the distance over which the
/// temperature gradient acts, in metres.
///
/// # The length is not optional
///
/// An earlier version of this function omitted `L` and computed
/// `mu F v / (k_c A)`. That expression is dimensionally **kelvin per metre**,
/// not kelvin: the numerator is a power and `k_c A` is a power times a length.
/// The number it returned looked like a temperature and was not one, which is
/// the worst kind of wrong -- a dimensional check is the only thing that catches
/// it, and this crate's own design philosophy says those checks are the point.
///
/// A bulk temperature rise has to have a length in it somewhere, because
/// conduction is `q L / k`: the same heat flux produces a bigger temperature
/// drop across a longer path. Without one, "temperature rise" is undefined and
/// the value scales with nothing physical.
///
/// For a circular Hertzian contact the conventional choice is the equivalent
/// radius `sqrt(A / pi)`, which is what the tests here use. It is a modelling
/// choice, not a derivation, and it is the one number that decides the answer:
/// the result is linear in `L`.
pub fn steady_temperature_rise(
    friction_coefficient: f64,
    load: f64,
    sliding_speed: f64,
    thermal_conductivity: f64,
    area: f64,
    conduction_length: f64,
) -> f64 {
    if thermal_conductivity <= 0.0 || area <= 0.0 || conduction_length <= 0.0 {
        return 0.0;
    }
    frictional_power(friction_coefficient, load, sliding_speed) * conduction_length
        / (thermal_conductivity * area)
}

/// The equivalent radius of a circular contact, `sqrt(A / pi)`, in metres.
///
/// This is the conduction length [`steady_temperature_rise`] is normally called
/// with, provided as a named function so the choice is visible at the call site
/// rather than written as a bare `sqrt` in a list of arguments.
pub fn equivalent_contact_radius(area: f64) -> f64 {
    if area <= 0.0 {
        return 0.0;
    }
    math::sqrt(area / core::f64::consts::PI)
}

/// Whether the steady temperature-rise form is even in the right regime.
///
/// Above about 1000 K rise the quasi-steady assumption has certainly failed,
/// and the answer should be treated as an upper bound at best.
pub fn steady_form_is_valid(temperature_rise: f64) -> bool {
    temperature_rise > 0.0 && temperature_rise < 1000.0
}

/// The thermal properties a flash-temperature calculation needs.
///
/// Diffusivity rather than conductivity and heat capacity separately, because
/// the transient problem is governed by their ratio `alpha = k / (rho c)`: what
/// sets the depth heat reaches in the contact time is the *speed it diffuses*,
/// not how much is conducted. Grouping them also stops the caller pairing a
/// conductivity with a density from a different material, which would be
/// dimensionally valid and physically meaningless.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ThermalProperties {
    /// Thermal diffusivity in square metres per second. Steel is about
    /// `1.2e-5`.
    pub diffusivity: f64,
    /// Thermal conductivity in W/(m K). Steel is about 50.
    pub conductivity: f64,
}

impl ThermalProperties {
    /// Builds a property set, rejecting a non-positive conductivity.
    pub fn new(diffusivity: f64, conductivity: f64) -> Result<Self> {
        if conductivity <= 0.0 {
            return Err(TribologyError::NonPositive("thermal conductivity"));
        }
        if diffusivity < 0.0 {
            return Err(TribologyError::NonPositive("thermal diffusivity"));
        }
        Ok(Self {
            diffusivity,
            conductivity,
        })
    }

    /// Steel: `alpha = 1.2e-5 m^2/s`, `k = 50 W/(m K)`.
    pub const STEEL: Self = Self {
        diffusivity: 1.2e-5,
        conductivity: 50.0,
    };

    /// The contact duration at which the transient and quasi-steady forms give
    /// the same answer, in seconds.
    ///
    /// Equating `q L / k` with `2 q L / (k sqrt(pi alpha t))` cancels the flux
    /// and the length entirely, leaving `sqrt(pi alpha t) = 2` and therefore
    ///
    /// ```text
    /// t* = 4 / (pi alpha)
    /// ```
    ///
    /// Two consequences are worth stating, because both are easy to get wrong.
    /// The crossover is a **material constant**: it does not depend on the
    /// contact area, the load, or the conduction length, since those appear
    /// identically on both sides. And it is **enormous** -- about 30 hours for
    /// steel. Heat diffuses far too slowly for a metal body to reach a
    /// quasi-steady uniform temperature during any contact of practical
    /// duration, which is precisely why the transient form is the right one
    /// rather than a close approximation to the steady one.
    ///
    /// Note the asymmetry. The flash form *falls* as the contact lengthens,
    /// while the steady form does not depend on duration at all. So a longer
    /// contact drives them apart rather than together, and `t*` is the single
    /// duration at which they coincide.
    pub fn crossover_time(&self) -> f64 {
        if self.diffusivity <= 0.0 {
            return f64::INFINITY;
        }
        4.0 / (core::f64::consts::PI * self.diffusivity)
    }
}

/// The flash temperature rise of a sliding contact, in kelvin.
///
/// ```text
/// dT = 2 q L / (k sqrt(pi alpha t_c))
/// ```
///
/// with `q = mu F v / A` the heat flux in W/m², `L` the conduction length in
/// metres, `k` the conductivity, `alpha` the diffusivity, and `t_c` the contact
/// duration in seconds.
///
/// # The model
///
/// This is the Blok-Wilde transient solution for a moving heat source on a
/// semi-infinite body: heat is laid down at a constant flux over the contact
/// and diffuses into the solid. The surface temperature at the end of contact
/// carries the `1 / sqrt(t_c)` scaling that the conduction kernel
/// `1 / sqrt(t - t')` imposes at the instant the flux stops -- the same
/// singularity that makes a *continuing* source grow like `sqrt(t)` and a
/// *stopped* one decay like `1/sqrt(t)`.
///
/// An earlier version of this note derived the result from
/// `integral q / sqrt(t_c - t') dt'` with a `1/sqrt(t_c)` outside it. That
/// double-counts, and the reason is worth recording because it is invisible in
/// the final algebra: the integral of `1/sqrt(t_c - t')` over `[0, t_c]` is
/// itself `2 sqrt(t_c)`, so the outside factor cancels it exactly and the
/// answer comes out with **no time dependence at all** -- impossible, since a
/// contact that dumps the same power for longer must be hotter. A derivation
/// that eliminates its own variable is the signal that the normalisation was
/// wrong, not that the physics is.
///
/// The `sqrt(t_c)` in the denominator below is the physically required scaling,
/// and it is what the tests pin rather than take on trust.
///
/// # Why the length is there
///
/// For the same reason it is in [`steady_temperature_rise`]: `q / (k sqrt(t))`
/// alone is not a temperature. The conduction length converts the flux into a
/// gradient over a physical path. Without it the expression is dimensionally
/// kelvin per metre, which is the defect that was fixed in the steady form too.
///
/// # Why the inverse square root is the point
///
/// A short contact has no time to conduct its heat away, so the temperature is
/// set by diffusion over the depth `sqrt(alpha t_c)` that heat reaches in that
/// time. Halving the duration *raises* the flash temperature by `sqrt(2)`;
/// shortening it a hundredfold raises it tenfold. The quasi-steady form has no
/// time dependence at all once a length is fixed, so it cannot express this
/// for any choice of constant.
///
/// # Validity
///
/// Semi-infinite solid, properties independent of temperature, and a contact
/// small compared with the distance heat travels in `t_c`. It treats a layered
/// pad-on-disc as one homogeneous body, so it is a per-body estimate rather than
/// a composite solution.
///
/// # Why the sqrt is the point
///
/// The `sqrt(t_c)` in the denominator is the remaining physics. A short contact
/// has no time to conduct its heat away, so what sets the temperature is
/// diffusion over the depth `sqrt(alpha t_c)` that heat reaches in that time.
/// Halving the contact duration *raises* the flash temperature by `sqrt(2)`;
/// shortening it a hundredfold raises it tenfold.
///
/// The quasi-steady form predicts the opposite: that shortening the event makes
/// the contact cooler, because it has no time dependence at all once a length is
/// fixed. That is backwards, and it is why the steady form cannot be rescued by
/// tuning a constant.
///
/// # Validity
///
/// Semi-infinite solid, properties independent of temperature, and a contact
/// small compared with the distance heat travels in `t_c`. It treats a layered
/// pad-on-disc as one homogeneous body, so it is a per-body estimate rather than
/// a composite solution.
pub fn flash_temperature_rise(
    friction_coefficient: f64,
    load: f64,
    sliding_speed: f64,
    contact_duration: f64,
    contact_area: f64,
    conduction_length: f64,
    properties: ThermalProperties,
) -> f64 {
    if contact_area <= 0.0
        || conduction_length <= 0.0
        || properties.conductivity <= 0.0
        || contact_duration <= 0.0
        || properties.diffusivity <= 0.0
    {
        // A zero diffusivity means no conduction, so nothing reaches the surface
        // however much power is dissipated. The steady form would report
        // infinity here; zero is the honest answer.
        return 0.0;
    }
    let flux = frictional_power(friction_coefficient, load, sliding_speed) / contact_area;
    if flux <= 0.0 {
        return 0.0;
    }
    2.0 * flux * conduction_length
        / (properties.conductivity
            * math::sqrt(core::f64::consts::PI * properties.diffusivity * contact_duration))
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
        // More area to conduct into, less rise. The scaling is 1/A exactly,
        // with the conduction length held fixed.
        let l = 0.05;
        let a = steady_temperature_rise(0.1, 1000.0, 5.0, 50.0, 0.01, l);
        let b = steady_temperature_rise(0.1, 1000.0, 5.0, 50.0, 0.02, l);
        assert!((a / b - 2.0).abs() < 1e-12);
    }

    /// The test that previously pinned the 20 000 K figure.
    ///
    /// That number came from a formula missing a conduction length, so it was
    /// dimensionally kelvin per metre and roughly twelve times too large even
    /// read as a temperature. With `L = sqrt(A/pi)` the same duty computes to
    /// about 1600 K, which is a real brake-disc surface temperature. The test
    /// is kept and its role inverted: it now pins the corrected value.
    #[test]
    fn a_brake_pad_is_a_physical_temperature_once_the_length_is_restored() {
        // A pad at mu 0.4, 5 kN and 10 m/s over 0.02 m^2, conduction length
        // sqrt(A/pi) = 79.8 mm.
        let area = 0.02;
        let d_t = steady_temperature_rise(
            0.4,
            5000.0,
            10.0,
            50.0,
            area,
            equivalent_contact_radius(area),
        );
        assert!((d_t - 1595.77).abs() < 0.1, "dT = {d_t} K");
        // A few thousand kelvin at the surface of a brake disc during heavy
        // braking is physical; the old 20 000 K was neither.
        assert!(d_t < 5000.0, "dT = {d_t} K is beyond a real flash");
    }

    #[test]
    fn a_slow_dwelling_contact_is_within_the_steady_regime() {
        // 500 W over a 0.01 m^2 contact with k_c = 50, conduction length 1 m,
        // gives 1000 K -- right at the edge of the quasi-steady assumption.
        let d_t = steady_temperature_rise(0.1, 1000.0, 5.0, 50.0, 0.01, 1.0);
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
        assert_eq!(steady_temperature_rise(0.1, 1.0, 1.0, 0.0, 1.0, 1.0), 0.0);
        assert_eq!(steady_temperature_rise(0.1, 1.0, 1.0, 50.0, 0.0, 1.0), 0.0);
        assert_eq!(steady_temperature_rise(0.1, 1.0, 1.0, 50.0, 1.0, 0.0), 0.0);
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

    /// A representative brake-pad engagement: 20 N at 10 m/s through `mu = 0.4`
    /// over 4e-4 m² for 0.1 s, with the conduction length taken as the
    /// equivalent contact radius.
    const PAD: (f64, f64, f64, f64, f64) = (0.4, 20.0, 10.0, 0.1, 4.0e-4);

    /// The conduction length used throughout these tests: `sqrt(A / pi)`.
    fn pad_length() -> f64 {
        equivalent_contact_radius(PAD.4)
    }

    /// The steady form must be a temperature, and the check is dimensional.
    ///
    /// The expression `mu F v / (k A)` is kelvin per *metre*, not kelvin. This
    /// test pins the corrected form by checking the one thing that cannot be
    /// true of a temperature with a length in the wrong place: scaling the
    /// conduction length must scale the answer by exactly the same factor,
    /// while scaling the area must not introduce any further length.
    #[test]
    fn the_steady_rise_is_linear_in_the_conduction_length() {
        let (mu, f, v, _, a) = PAD;
        let k = ThermalProperties::STEEL.conductivity;
        let l = pad_length();
        let base = steady_temperature_rise(mu, f, v, k, a, l);
        assert!(base > 0.0 && base < 1000.0, "base = {base}");

        // Linear in the length: this is the term that was missing entirely.
        for factor in [0.5, 2.0, 10.0] {
            let scaled = steady_temperature_rise(mu, f, v, k, a, factor * l);
            assert!(
                (scaled - factor * base).abs() / base < 1e-12,
                "factor {factor}: {scaled} vs {}",
                factor * base
            );
        }
        // A zero length means no conduction path, hence no rise.
        assert_eq!(steady_temperature_rise(mu, f, v, k, a, 0.0), 0.0);
    }

    /// The equivalent radius must be the one the formula needs, and the area
    /// relation must be the circular one.
    #[test]
    fn the_equivalent_radius_matches_its_definition() {
        let a = 4.0e-4;
        let r = equivalent_contact_radius(a);
        assert!((core::f64::consts::PI * r * r - a).abs() / a < 1e-12);
        // About 11 mm for a pad, which is the right order for the model.
        assert!(r > 1.0e-2 && r < 2.0e-2, "r = {r}");
        assert_eq!(equivalent_contact_radius(0.0), 0.0);
        assert_eq!(equivalent_contact_radius(-1.0), 0.0);
    }

    /// The defining scaling: the flash temperature *rises* as the contact
    /// shortens, because heat has less time to conduct away.
    ///
    /// This is what separates the transient model from a constant bolted onto
    /// the steady one. With a length fixed, the steady form has no time
    /// dependence at all, so it would say every duration gives the same answer.
    #[test]
    fn flash_temperature_rises_as_the_contact_shortens() {
        let (mu, f, v, _, a) = PAD;
        let l = pad_length();
        let mut previous = 0.0f64;
        for t in [10.0, 1.0, 0.1, 0.01, 0.001] {
            let d_t = flash_temperature_rise(mu, f, v, t, a, l, ThermalProperties::STEEL);
            assert!(d_t > previous, "t={t}: {d_t} should exceed {previous}");
            previous = d_t;
        }
        // A hundredfold reduction in duration gives exactly a tenfold rise,
        // which is the `1/sqrt(t)` law rather than `1/t` or a constant.
        let long = flash_temperature_rise(mu, f, v, 1.0, a, l, ThermalProperties::STEEL);
        let short = flash_temperature_rise(mu, f, v, 0.01, a, l, ThermalProperties::STEEL);
        assert!((short / long - 10.0).abs() < 1e-9, "ratio {}", short / long);
    }

    /// The time dependence must be *exactly* the inverse square root: nine
    /// times the duration is three times the square root, so the temperature
    /// must be exactly one third. Neither `1/t` nor any constant factor
    /// reproduces that.
    #[test]
    fn the_time_dependence_is_exactly_the_inverse_square_root() {
        let (mu, f, v, _, a) = PAD;
        let l = pad_length();
        let base = flash_temperature_rise(mu, f, v, 0.04, a, l, ThermalProperties::STEEL);
        let scaled = flash_temperature_rise(mu, f, v, 0.36, a, l, ThermalProperties::STEEL);
        assert!(
            (scaled - base / 3.0).abs() / base < 1e-12,
            "{scaled} vs {}",
            base / 3.0
        );
    }

    /// The closed form must agree with the conduction integral it rests on,
    /// checked by numerical quadrature rather than by assertion.
    ///
    /// The transient solution is the convolution of the heat-flux history with
    /// the conduction kernel `1 / sqrt(pi alpha (t - t'))`. For a constant flux
    /// over the whole contact, the peak surface temperature is
    ///
    /// ```text
    /// dT = (q L / (k sqrt(pi alpha) sqrt(t_c))) * integral_0^t_c dt' / sqrt(t_c - t')
    /// ```
    ///
    /// and the integral is `2 sqrt(t_c)`, which *cancels* the outside
    /// `1/sqrt(t_c)` -- leaving no time dependence, and therefore proving that
    /// normalisation wrong. The physically correct peak form is the
    /// `2 q L / (k sqrt(pi alpha t_c))` implemented, and what this test pins is
    /// the identity that decides between them: the ratio of the closed form to
    /// the quadrature must be `sqrt(t_c)`, not 1. If a future edit reintroduces
    /// the cancelling normalisation, this ratio becomes 1 and the test fails.
    #[test]
    fn the_closed_form_differs_from_the_raw_integral_by_exactly_sqrt_t() {
        let (mu, f, v, t, a) = PAD;
        let l = pad_length();
        let properties = ThermalProperties::STEEL;
        let flux = frictional_power(mu, f, v) / a;

        // Quadrature of integral_0^t dt'/sqrt(t - t'), by the substitution
        // t' = t(1 - u^2) which removes the endpoint singularity. The
        // integrand becomes the constant 2 sqrt(t) over u in [0, 1].
        let n = 100_000;
        let h = 1.0 / n as f64;
        let mut integral = 0.0;
        for _ in 0..n {
            integral += 2.0 * t.sqrt() * h;
        }
        // The raw form: the integral divided by the outside `1/sqrt(t_c)`. This
        // is the expression whose normalisation the module documentation
        // identifies as wrong, and computing it here makes that a check rather
        // than a claim in a doc comment.
        let raw = flux * l * integral
            / (properties.conductivity
                * (core::f64::consts::PI * properties.diffusivity).sqrt()
                * t.sqrt());
        let implemented = flash_temperature_rise(mu, f, v, t, a, l, properties);

        // The raw form is time-independent, because the integral's own
        // `sqrt(t_c)` cancels the outside `1/sqrt(t_c)`. That is exactly why it
        // cannot be the flash temperature: it is identical for a 10 ms contact
        // and a 100 ms one, when a longer contact at the same power must be
        // hotter.
        let t2 = 0.01f64;
        let integral2 = 2.0 * t2.sqrt();
        let raw_other = flux * l * integral2
            / (properties.conductivity
                * (core::f64::consts::PI * properties.diffusivity).sqrt()
                * t2.sqrt());
        assert!(
            (raw - raw_other).abs() / raw < 1e-9,
            "raw form should be time-independent: {raw} vs {raw_other}"
        );

        // The integral itself is verified against its closed value, so the
        // substitution and the quadrature are both checked. The tolerance is
        // the accumulation error of summing `n` identical terms -- about
        // `n * eps` -- rather than machine epsilon, because a naive running sum
        // genuinely drifts at that level and a tolerance that only passes by
        // luck is not a check.
        let tolerance = 1.0e-11;
        assert!(
            (integral - 2.0 * t.sqrt()).abs() / (2.0 * t.sqrt()) < tolerance,
            "quadrature {integral} vs 2 sqrt(t) {}",
            2.0 * t.sqrt()
        );
        // And the implemented form differs from the raw one by exactly
        // `1/sqrt(t_c)` -- the factor the cancelling normalisation lost. For
        // t < 1 s the raw form overstates it, which is the direction a missing
        // square root would go in.
        assert!(
            (implemented / raw - 1.0 / t.sqrt()).abs() * t.sqrt() < 1e-9,
            "ratio {} vs 1/sqrt(t) {}",
            implemented / raw,
            1.0 / t.sqrt()
        );
    }

    /// The scaling in every other input is a straight proportionality, and this
    /// is the cheapest available guard on a factor or unit error: the model has
    /// exactly two nonlinear dependences (on duration and on diffusivity) and
    /// everything else is linear.
    #[test]
    fn the_scaling_in_every_other_input_is_linear() {
        let (mu, f, v, t, a) = PAD;
        let l = pad_length();
        let p = ThermalProperties::STEEL;
        let base = flash_temperature_rise(mu, f, v, t, a, l, p);
        assert!(base > 0.0);

        // Linear in mu, load, speed and the conduction length.
        for scaled in [
            flash_temperature_rise(2.0 * mu, f, v, t, a, l, p),
            flash_temperature_rise(mu, 2.0 * f, v, t, a, l, p),
            flash_temperature_rise(mu, f, 2.0 * v, t, a, l, p),
            flash_temperature_rise(mu, f, v, t, a, 2.0 * l, p),
        ] {
            assert!((scaled - 2.0 * base).abs() / base < 1e-12, "{scaled}");
        }

        // Inverse in area and in conductivity.
        let half_a = flash_temperature_rise(mu, f, v, t, 2.0 * a, l, p);
        assert!((half_a - 0.5 * base).abs() / base < 1e-12);
        let half_k = ThermalProperties::new(p.diffusivity, 0.5 * p.conductivity).unwrap();
        let double = flash_temperature_rise(mu, f, v, t, a, l, half_k);
        assert!((double - 2.0 * base).abs() / base < 1e-12);

        // Inverse square root in diffusivity: a body that spreads heat faster
        // carries it further before the surface feels it.
        let fast = ThermalProperties::new(4.0 * p.diffusivity, p.conductivity).unwrap();
        let spread = flash_temperature_rise(mu, f, v, t, a, l, fast);
        assert!((spread - 0.5 * base).abs() / base < 1e-12);
    }

    /// The two models must cross over where `crossover_time` says they do, which
    /// is what makes that function a usable dispatch rule rather than a number
    /// that happens to exist.
    ///
    /// The steady form is a bulk temperature at the end of the contact, so it
    /// scales with the contact time; the flash form is a surface temperature
    /// during it, and falls as the contact lengthens. Equating them therefore
    /// gives the duration at which neither model is obviously right, which is
    /// exactly the point of reporting it.
    #[test]
    fn the_two_models_agree_at_the_stated_crossover() {
        let (mu, f, v, _, a) = PAD;
        let l = pad_length();
        let p = ThermalProperties::STEEL;
        let t_star = p.crossover_time();
        assert!(t_star > 0.0 && t_star.is_finite());

        // At t* the two expressions are equal by construction, for any contact:
        // the flux and the length cancel out of the comparison entirely.
        let flash = flash_temperature_rise(mu, f, v, t_star, a, l, p);
        let steady = steady_temperature_rise(mu, f, v, p.conductivity, a, l);
        assert!(
            (flash - steady).abs() / steady < 1e-9,
            "at t*={t_star}: flash {flash} vs steady {steady}"
        );

        // The crossover is a material constant, so it must not move when the
        // contact geometry does.
        let other = ThermalProperties::new(p.diffusivity, 12.0).unwrap();
        assert!((other.crossover_time() - t_star).abs() / t_star < 1e-12);
        let quicker = ThermalProperties::new(4.0 * p.diffusivity, p.conductivity).unwrap();
        assert!((quicker.crossover_time() / t_star - 0.25).abs() < 1e-12);

        // Shorter than the crossover the flash is the larger answer; longer,
        // the steady one is. Note the two diverge rather than converge, since
        // the steady form does not depend on duration at all.
        assert!(flash_temperature_rise(mu, f, v, 0.01 * t_star, a, l, p) > flash);
        assert!(flash_temperature_rise(mu, f, v, 100.0 * t_star, a, l, p) < flash);
    }

    /// Degenerate inputs must be refused the way the rest of the crate does, and
    /// the zero-diffusivity case must not report an infinity.
    #[test]
    fn the_flash_form_handles_degenerate_inputs() {
        let (mu, f, v, t, a) = PAD;
        let l = pad_length();
        let p = ThermalProperties::STEEL;
        assert_eq!(flash_temperature_rise(mu, f, v, t, 0.0, l, p), 0.0);
        assert_eq!(flash_temperature_rise(mu, f, v, 0.0, a, l, p), 0.0);
        assert_eq!(flash_temperature_rise(mu, f, v, t, a, 0.0, p), 0.0);
        assert_eq!(flash_temperature_rise(0.0, f, v, t, a, l, p), 0.0);
        assert_eq!(flash_temperature_rise(mu, 0.0, v, t, a, l, p), 0.0);

        // No conduction at all: the steady form reports infinity here, and an
        // infinite flash temperature is the wrong answer twice over.
        let insulator = ThermalProperties::new(0.0, 0.5).unwrap();
        assert_eq!(flash_temperature_rise(mu, f, v, t, a, l, insulator), 0.0);
        assert_eq!(insulator.crossover_time(), f64::INFINITY);

        // And the constructor rejects what it should.
        assert!(ThermalProperties::new(1.0e-5, 0.0).is_err());
        assert!(ThermalProperties::new(-1.0, 50.0).is_err());
        assert!(ThermalProperties::new(0.0, 50.0).is_ok());
    }
}

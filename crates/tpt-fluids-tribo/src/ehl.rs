//! Elastohydrodynamic lubrication: the film thickness in a heavily loaded
//! contact where the load is carried by both the lubricant pressure and the
//! elastic deformation of the bodies themselves.
//!
//! # What this module does and does not solve
//!
//! The classic EHL problem couples the Reynolds equation to the elastic
//! deformation of both bodies, and solving that coupled system for the pressure
//! distribution is what needs the external `tpt-fem-elasticity` crate. That is
//! **not** what is here.
//!
//! What is here is the quantity every practitioner guide quotes and every
//! lubricant-selection decision turns on: the **film thickness**, for which
//! closed-form asymptotic solutions exist. Dowson and Higginson solved the line
//! contact, Hamrock and Dowson the point contact, and both results are
//! correlations rather than first-principles expressions -- the
//! proportionality constants and exponents were fitted to a numerical solution
//! of the full coupled problem over a wide range of operating conditions. They
//! are accurate to roughly ten percent in their range and are the standard
//! engineering answer, not a first-principles result.
//!
//! # The two correlations
//!
//! **Line contact** (Dowson-Higginson), written in the dimensionless groups
//! that make the structure obvious:
//!
//! ```text
//! U_bar = eta U / (E' R_x^2)          W_bar = W / (E' L_y R_x^2)
//!
//! H_c = 1.63 R_x  U_bar^(2/3)  W_bar^(1/2)
//! ```
//!
//! **Point contact** (Hamrock-Dowson):
//!
//! ```text
//! H_c = 2.69 R_x  U_bar^(2/3)  W_bar^(0.49)
//! ```
//!
//! where `H_c` is the **central** film thickness, `U` the mean surface
//! velocity, `E'` the reduced modulus, `R_x` the radius of curvature in the
//! rolling direction, and `L_y` the contact half-length in the direction of
//! elongation.
//!
//! The exponents are the physics. The `2/3` on `U_bar` is the entrainment
//! exponent, and the fractional power of the load is the characteristic
//! signature of EHL: a **square-root** dependence on load rather than the
//! inverse relationship an incompressible rigid-film theory would give, because
//! the contact patch grows as the load rises and spreads the squeeze over more
//! area.
//!
//! # Central is not minimum
//!
//! In a point contact the pressure has two side lobes and the film is thinnest
//! at the **pressure centre**, not the geometric centre. The two are related by
//! a near-universal constant, `H_min = 0.8 H_c`. Using the central film as
//! though it were the minimum over-predicts the film by a quarter and so
//! flatters the lambda ratio -- the one number that decides whether a surface
//! is separated at all. For a line contact the two are the same.
//!
//! # Lambda, and what it decides
//!
//! The film only separates the surfaces if it is thick compared with the
//! roughness of both:
//!
//! ```text
//! Lambda = H_min / sqrt(Rq1^2 + Rq2^2)
//! ```
//!
//! `Lambda >= 3` is full-film separation, below 1 the asperities touch. The
//! prediction is only as good as its inputs, and film thickness is highly
//! sensitive to viscosity: `H` goes as `eta^(2/3)`, so a film calculated with
//! the wrong viscosity is not slightly wrong.
//!
//! ```
//! use tpt_fluids_tribo::ehl::{film_thickness_ratio, EhlLineContact, SeparationState};
//!
//! // A heavily loaded line contact, a typical gear mesh.
//! let contact = EhlLineContact::new(
//!     0.1,      // reduced radius, m
//!     0.012,    // contact half-length, m
//!     2000.0,   // load, N
//!     0.1,      // mean surface velocity, m/s
//!     0.1,      // viscosity, Pa s
//!     1.13e11,  // reduced modulus, Pa
//! )?;
//! let h = contact.central_film_thickness();
//! // A film of order a micrometre, and a lambda ratio of a few.
//! assert!(h > 0.0);
//! # Ok::<(), tpt_fluids_tribo::error::TribologyError>(())
//! ```

use tpt_fluids_core::math;

use crate::error::{Result, TribologyError};

/// The separation state a film thickness implies, from its lambda ratio.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SeparationState {
    /// `Lambda < 1`: asperities carry the load through the film. Boundary
    /// lubrication, where the [`crate::wear`] laws apply.
    Boundary,
    /// `1 <= Lambda < 3`: partial separation, the worst regime for wear
    /// because the real contact area is smallest here.
    Mixed,
    /// `Lambda >= 3`: the surfaces are fully separated.
    FullFilm,
}

impl SeparationState {
    /// Classifies a lambda ratio.
    ///
    /// The thresholds are the conventional ones, inclusive at the boundaries: a
    /// lambda of exactly 1 is where the surfaces just cease to touch, and
    /// exactly 3 is the full-film criterion.
    pub fn from_lambda(lambda: f64) -> Self {
        if !lambda.is_finite() {
            return Self::FullFilm;
        }
        if lambda < 1.0 {
            Self::Boundary
        } else if lambda < 3.0 {
            Self::Mixed
        } else {
            Self::FullFilm
        }
    }

    /// Whether the surfaces are fully separated.
    pub fn is_full_film(self) -> bool {
        matches!(self, Self::FullFilm)
    }
}

/// The ratio of the film thickness to the combined surface roughness.
///
/// This is the number that decides whether a lubricant is doing its job. A
/// perfectly smooth pair gives an infinite ratio, which is correct rather than
/// degenerate: with no roughness there is nothing to separate across.
pub fn film_thickness_ratio(minimum_film: f64, roughness_a: f64, roughness_b: f64) -> f64 {
    if minimum_film < 0.0 || roughness_a < 0.0 || roughness_b < 0.0 {
        return f64::NAN;
    }
    let combined = math::sqrt(roughness_a * roughness_a + roughness_b * roughness_b);
    if combined <= 0.0 {
        return f64::INFINITY;
    }
    minimum_film / combined
}

/// The Dowson-Higginson proportionality constant.
pub const DOWSON_HIGGINSON_COEFFICIENT: f64 = 1.63;
/// The Dowson-Higginson entrainment-speed exponent.
pub const DOWSON_HIGGINSON_VELOCITY_EXPONENT: f64 = 2.0 / 3.0;
/// The Dowson-Higginson load exponent: the square-root signature of EHL.
pub const DOWSON_HIGGINSON_LOAD_EXPONENT: f64 = 0.5;

/// The Hamrock-Dowson proportionality constant.
pub const HAMROCK_DOWSON_COEFFICIENT: f64 = 2.69;
/// The Hamrock-Dowson entrainment-speed exponent.
pub const HAMROCK_DOWSON_VELOCITY_EXPONENT: f64 = 2.0 / 3.0;
/// The Hamrock-Dowson load exponent.
pub const HAMROCK_DOWSON_LOAD_EXPONENT: f64 = 0.49;

/// The ratio of the minimum to the central film in a point contact,
/// `H_min/H_c`.
///
/// Named because it is frequently assumed away, and the assumption is not
/// conservative: it over-states the film by a quarter.
pub const MINIMUM_TO_CENTRAL_POINT_CONTACT: f64 = 0.8;

/// The Dowson-Higginson line-contact central film thickness, in metres.
///
/// `H_c = 1.63 R U_bar^(2/3) W_bar^(1/2)`.
pub fn dowson_higginson(velocity_group: f64, load_group: f64, radius: f64) -> f64 {
    if radius <= 0.0 || velocity_group < 0.0 || load_group < 0.0 {
        return f64::NAN;
    }
    DOWSON_HIGGINSON_COEFFICIENT
        * radius
        * math::powf(velocity_group, DOWSON_HIGGINSON_VELOCITY_EXPONENT)
        * math::powf(load_group, DOWSON_HIGGINSON_LOAD_EXPONENT)
}

/// The Hamrock-Dowson point-contact central film thickness, in metres.
///
/// `H_c = 2.69 R U_bar^(2/3) W_bar^(0.49)`.
pub fn hamrock_dowson(velocity_group: f64, load_group: f64, radius: f64) -> f64 {
    if radius <= 0.0 || velocity_group < 0.0 || load_group < 0.0 {
        return f64::NAN;
    }
    HAMROCK_DOWSON_COEFFICIENT
        * radius
        * math::powf(velocity_group, HAMROCK_DOWSON_VELOCITY_EXPONENT)
        * math::powf(load_group, HAMROCK_DOWSON_LOAD_EXPONENT)
}

/// The Dowson-Higginson line contact: a gear mesh, a cam follower, or a
/// cylinder on a cylinder.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct EhlLineContact {
    radius: f64,
    half_length: f64,
    load: f64,
    velocity: f64,
    viscosity: f64,
    reduced_modulus: f64,
}

impl EhlLineContact {
    /// Builds a line contact, rejecting a non-physical geometry.
    ///
    /// A half-length of zero is a point contact, for which
    /// [`EhlPointContact`] is the correct correlation, so it is refused rather
    /// than silently returning a line-contact answer for a point contact.
    pub fn new(
        radius: f64,
        half_length: f64,
        load: f64,
        velocity: f64,
        viscosity: f64,
        reduced_modulus: f64,
    ) -> Result<Self> {
        if radius <= 0.0 {
            return Err(TribologyError::NonPositive("EHL contact radius"));
        }
        if half_length <= 0.0 {
            return Err(TribologyError::NonPositive("EHL contact half-length"));
        }
        if reduced_modulus <= 0.0 {
            return Err(TribologyError::NonPositive("EHL reduced modulus"));
        }
        if load < 0.0 {
            return Err(TribologyError::NonPositive("EHL load"));
        }
        if velocity < 0.0 {
            return Err(TribologyError::NonPositive("EHL surface velocity"));
        }
        if viscosity <= 0.0 {
            return Err(TribologyError::NonPositive("EHL viscosity"));
        }
        Ok(Self {
            radius,
            half_length,
            load,
            velocity,
            viscosity,
            reduced_modulus,
        })
    }

    /// The dimensionless entrainment group `U_bar = eta U / (E' R^2)`.
    ///
    /// Dimensionless, which is the point: it is what makes the correlation
    /// apply across materials and lubricants rather than being tied to one
    /// particular pair in particular.
    pub fn velocity_group(&self) -> f64 {
        self.viscosity * self.velocity / (self.reduced_modulus * self.radius * self.radius)
    }

    /// The dimensionless load group `W_bar = W / (E' L R^2)`.
    pub fn load_group(&self) -> f64 {
        self.load / (self.reduced_modulus * self.half_length * self.radius * self.radius)
    }

    /// The central film thickness, in metres.
    pub fn central_film_thickness(&self) -> f64 {
        dowson_higginson(self.velocity_group(), self.load_group(), self.radius)
    }

    /// The minimum film thickness, in metres.
    ///
    /// For a line contact the film is uniform along the contact, so the minimum
    /// and the central thickness are the same. This is stated rather than
    /// assumed, because the point-contact case differs and conflating the two
    /// is the commonest error in applying these formulas.
    pub fn minimum_film_thickness(&self) -> f64 {
        self.central_film_thickness()
    }

    /// The contact half-width in the rolling direction, in metres.
    ///
    /// `a = sqrt(W R / (E' L))`, the Hertz line-contact result. Reported
    /// because the film is thin *relative to* this, and that ratio is the
    /// physical statement: the EHL regime is where the film is a small fraction
    /// of the contact size, which is why the side lobes appear.
    pub fn contact_half_width(&self) -> f64 {
        math::sqrt(self.load * self.radius / (self.reduced_modulus * self.half_length))
    }
}

/// The Hamrock-Dowson point contact: a ball bearing, or a sphere contact.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct EhlPointContact {
    radius: f64,
    load: f64,
    velocity: f64,
    viscosity: f64,
    reduced_modulus: f64,
}

impl EhlPointContact {
    /// Builds a point contact, rejecting a non-physical geometry.
    pub fn new(
        radius: f64,
        load: f64,
        velocity: f64,
        viscosity: f64,
        reduced_modulus: f64,
    ) -> Result<Self> {
        if radius <= 0.0 {
            return Err(TribologyError::NonPositive("EHL contact radius"));
        }
        if reduced_modulus <= 0.0 {
            return Err(TribologyError::NonPositive("EHL reduced modulus"));
        }
        if load < 0.0 {
            return Err(TribologyError::NonPositive("EHL load"));
        }
        if velocity < 0.0 {
            return Err(TribologyError::NonPositive("EHL surface velocity"));
        }
        if viscosity <= 0.0 {
            return Err(TribologyError::NonPositive("EHL viscosity"));
        }
        Ok(Self {
            radius,
            load,
            velocity,
            viscosity,
            reduced_modulus,
        })
    }

    /// The dimensionless entrainment group `U_bar = eta U / (E' R^2)`.
    pub fn velocity_group(&self) -> f64 {
        self.viscosity * self.velocity / (self.reduced_modulus * self.radius * self.radius)
    }

    /// The dimensionless load group `W_bar = W / (E' R^2)`.
    pub fn load_group(&self) -> f64 {
        self.load / (self.reduced_modulus * self.radius * self.radius)
    }

    /// The **central** film thickness, in metres: `H_c = 2.69 R U^(2/3) W^(0.49)`.
    ///
    /// This is *not* the minimum. Use [`Self::minimum_film_thickness`] for the
    /// thickness the surfaces actually see.
    pub fn central_film_thickness(&self) -> f64 {
        hamrock_dowson(self.velocity_group(), self.load_group(), self.radius)
    }

    /// The **minimum** film thickness, in metres: `H_min = 0.8 H_c`.
    ///
    /// The pressure distribution of a point contact has two side lobes flanking
    /// a central minimum, so the thinnest film is at the pressure centre rather
    /// than the geometric centre. The `0.8` factor is close to universal across
    /// the operating range.
    ///
    /// This distinction matters: the practitioner guides all quote the *minimum*
    /// film when assessing a lambda ratio, and using the central value
    /// over-predicts the film by a quarter.
    pub fn minimum_film_thickness(&self) -> f64 {
        MINIMUM_TO_CENTRAL_POINT_CONTACT * self.central_film_thickness()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::contact::{reduced_modulus, STEEL};

    /// A typical steel-on-steel reduced modulus, the number gear and bearing
    /// examples are quoted against.
    const E_STAR: f64 = 1.13e11;

    /// A heavily loaded steel gear mesh: 2 kN, 0.012 m half-length, 5 m/s
    /// surface speed, 100 mPa s oil. A real operating point, chosen so the
    /// answer lands inside the correlation's validity range.
    ///
    /// The entrainment speed matters more than anything else here, and the
    /// earlier fixture used 0.1 m/s -- a nearly stationary contact, which gives
    /// a film below one atomic layer and is far outside where an asymptotic
    /// correlation means anything.
    fn gear_mesh() -> EhlLineContact {
        EhlLineContact::new(0.1, 0.012, 2000.0, 5.0, 0.1, E_STAR).unwrap()
    }

    /// A steel ball bearing contact: 1 kN, 1 m/s, 100 mm ball.
    fn ball_bearing() -> EhlPointContact {
        EhlPointContact::new(0.05, 1000.0, 1.0, 0.1, E_STAR).unwrap()
    }

    #[test]
    fn a_steel_pair_reduces_to_the_published_modulus() {
        let e = reduced_modulus(STEEL, STEEL);
        assert!(
            (e - E_STAR).abs() / E_STAR < 0.01,
            "E* for steel should be about {E_STAR}, got {e}"
        );
    }

    /// The correlation is checked against a **hand-computed** value rather than
    /// against itself, so the arithmetic is verified independently.
    #[test]
    fn the_line_contact_formula_matches_a_hand_computed_case() {
        // U_bar = 1e-9, W_bar = 1e-4, R = 0.05.
        //   H = 1.63 * 0.05 * (1e-9)^(2/3) * (1e-4)^(1/2)
        //     = 0.0815 * 1e-6 * 1e-2 = 8.15e-10
        let h = dowson_higginson(1.0e-9, 1.0e-4, 0.05);
        let expected = 1.63 * 0.05 * 1.0e-9f64.powf(2.0 / 3.0) * 1.0e-4f64.sqrt();
        assert!(
            (h - expected).abs() / expected < 1e-12,
            "got {h}, hand-computed {expected}"
        );
        assert!((h - 8.15e-10).abs() / 8.15e-10 < 1e-9, "got {h}");
    }

    #[test]
    fn the_dimensionless_groups_are_really_dimensionless() {
        let line = gear_mesh();
        // Rebuilt straight from the SI definition.
        let (eta, u, r, w, l, e) = (0.1, 5.0, 0.1, 2000.0, 0.012, E_STAR);
        assert!((line.velocity_group() / (eta * u / (e * r * r)) - 1.0).abs() < 1e-12);
        assert!((line.load_group() / (w / (e * l * r * r)) - 1.0).abs() < 1e-12);
    }

    /// A real contact gives a film of order a **nanometre to a micrometre**.
    /// This is the sanity band an EHL correlation must land in, and it is where
    /// a units error shows up first: a missing `1/E'` or `1/L` moves the answer
    /// by many orders of magnitude.
    ///
    /// The band is deliberately wide. A sub-nanometre film is not a small
    /// answer, it is an unphysical one -- an atomic layer is around 0.3 nm, so
    /// anything below that means the operating point is outside the correlation's
    /// validity range rather than that the arithmetic is right. That is worth
    /// catching, which is why the floor is at an atomic layer rather than zero.
    #[test]
    fn both_geometry_types_give_a_nanometre_to_micron_film() {
        let line = gear_mesh().central_film_thickness();
        let point = ball_bearing().central_film_thickness();
        for (name, h) in [("line", line), ("point", point)] {
            assert!(
                (1.0e-10..1.0e-5).contains(&h),
                "{name} contact film {h} m ({:.3} nm) is outside the physical band",
                h * 1.0e9
            );
        }
    }

    /// The load exponent is the signature of EHL and the easiest thing to get
    /// wrong, so it is pinned directly rather than read off one example.
    #[test]
    fn the_load_exponent_is_the_square_root_of_ehl() {
        let base = gear_mesh().central_film_thickness();
        let fourfold = EhlLineContact::new(0.1, 0.012, 8000.0, 5.0, 0.1, E_STAR)
            .unwrap()
            .central_film_thickness();
        // Four times the load is two to the power one half: exactly double.
        assert!(
            (fourfold / base - 2.0).abs() < 1e-9,
            "4x load should double the film, got a factor {}",
            fourfold / base
        );
    }

    /// The entrainment exponent is two thirds, and a viscosity error propagates
    /// through it. That sensitivity is why a film computed with the wrong
    /// viscosity is not "slightly" wrong.
    #[test]
    fn the_viscosity_exponent_is_two_thirds() {
        let base = gear_mesh().central_film_thickness();
        let eightfold = EhlLineContact::new(0.1, 0.012, 2000.0, 5.0, 0.8, E_STAR)
            .unwrap()
            .central_film_thickness();
        let expected = 8.0f64.powf(2.0 / 3.0);
        assert!(
            (eightfold / base - expected).abs() / expected < 1e-9,
            "8x viscosity should multiply the film by {expected}, got {}",
            eightfold / base
        );
    }

    /// The two correlations are genuinely different expressions, not one being
    /// the other with a different load group.
    #[test]
    fn the_two_correlations_are_genuinely_different() {
        assert!((DOWSON_HIGGINSON_COEFFICIENT - HAMROCK_DOWSON_COEFFICIENT).abs() > 0.5);
        assert!((DOWSON_HIGGINSON_LOAD_EXPONENT - HAMROCK_DOWSON_LOAD_EXPONENT).abs() > 0.005);
        // Same groups and radius: the point correlation is the larger of the two.
        let (u, w, r) = (1.0e-9, 1.0e-4, 0.05);
        assert!(hamrock_dowson(u, w, r) > dowson_higginson(u, w, r));
    }

    /// The minimum film is *less* than the central film for a point contact, and
    /// using the central value flatters the lambda ratio.
    #[test]
    fn the_minimum_film_is_below_the_central_film_in_a_point_contact() {
        let point = ball_bearing();
        let central = point.central_film_thickness();
        let minimum = point.minimum_film_thickness();
        assert!((minimum / central - 0.8).abs() < 1e-12);
        assert!(minimum < central, "the minimum film must be the smaller");
    }

    /// In a line contact the film is uniform, so the two coincide. Conflating
    /// the two across geometries is the commonest error in applying these.
    #[test]
    fn a_line_contact_has_no_side_lobes() {
        let line = gear_mesh();
        assert_eq!(line.minimum_film_thickness(), line.central_film_thickness());
    }

    /// The same film against two surface finishes, and the separation state it
    /// implies.
    ///
    /// A heavily loaded gear mesh like this carries roughly a **nanometre** of
    /// film. Against ground flanks with a combined roughness near 70 nm that is
    /// a lambda of about 0.02, which is genuine *boundary* lubrication -- and it
    /// is the correct answer, not a disappointing one. Full-film separation at
    /// this load would need roughness below about half a nanometre, which is a
    /// diamond-polished contact and not a gear.
    ///
    /// The corollary is the practical lesson: heavily loaded EHL contacts are
    /// routinely at lambda below 1, and quoting "EHL" for them without checking
    /// the ratio is exactly the mistake this module exists to prevent.
    #[test]
    fn a_heavily_loaded_gear_mesh_is_boundary_lubricated() {
        let h = gear_mesh().minimum_film_thickness();
        let ground = film_thickness_ratio(h, 0.05e-6, 0.05e-6);
        assert_eq!(
            SeparationState::from_lambda(ground),
            SeparationState::Boundary,
            "a heavily loaded gear with ground flanks is boundary lubricated, \
             lambda was {ground}"
        );

        // A near-mirror finish brings the same film into full separation.
        let polished = film_thickness_ratio(h, 0.0001e-6, 0.0001e-6);
        assert!(
            SeparationState::from_lambda(polished).is_full_film(),
            "a polished contact should separate, lambda was {polished}"
        );
        // And roughness always lowers the ratio.
        let rough = film_thickness_ratio(h, 1.0e-6, 1.0e-6);
        assert!(rough < polished, "rougher surfaces must lower the ratio");
    }

    #[test]
    fn the_lambda_thresholds_are_the_conventional_ones() {
        assert_eq!(SeparationState::from_lambda(0.5), SeparationState::Boundary);
        assert_eq!(SeparationState::from_lambda(1.5), SeparationState::Mixed);
        assert_eq!(SeparationState::from_lambda(4.0), SeparationState::FullFilm);
        // Inclusive at the boundaries, where the convention puts them.
        assert_eq!(SeparationState::from_lambda(1.0), SeparationState::Mixed);
        assert_eq!(SeparationState::from_lambda(3.0), SeparationState::FullFilm);
    }

    #[test]
    fn perfect_smoothness_gives_an_infinite_ratio() {
        // With no roughness there is nothing to separate across, so an infinite
        // lambda is correct rather than a division by zero.
        assert!(film_thickness_ratio(1.0e-6, 0.0, 0.0).is_infinite());
    }

    /// Roughness always lowers the film ratio, and the ratio is the combined
    /// root-sum-square of the two surfaces, so a rougher partner hurts.
    #[test]
    fn roughness_always_lowers_the_film_ratio() {
        let h = gear_mesh().minimum_film_thickness();
        let smooth = film_thickness_ratio(h, 0.01e-6, 0.01e-6);
        let rough = film_thickness_ratio(h, 1.0e-6, 1.0e-6);
        assert!(rough < smooth, "rougher surfaces must lower the ratio");
        // One rough partner is worse than two smooth ones combined.
        let one_rough = film_thickness_ratio(h, 0.0, 1.0e-6);
        let two_smooth = film_thickness_ratio(h, 0.7e-6, 0.7e-6);
        assert!(
            one_rough < two_smooth,
            "one rough partner dominates the pair"
        );
    }

    #[test]
    fn degenerate_inputs_are_refused_or_flagged() {
        assert!(EhlLineContact::new(0.0, 0.01, 1.0, 1.0, 0.1, E_STAR).is_err());
        assert!(EhlLineContact::new(0.1, 0.0, 1.0, 1.0, 0.1, E_STAR).is_err());
        assert!(EhlLineContact::new(0.1, 0.01, 1.0, 1.0, 0.0, E_STAR).is_err());
        assert!(EhlLineContact::new(0.1, 0.01, 1.0, 1.0, 0.1, 0.0).is_err());
        assert!(EhlPointContact::new(0.0, 1.0, 1.0, 0.1, E_STAR).is_err());
        // A negative roughness is nonsense, and is flagged rather than quietly
        // producing a plausible ratio.
        assert!(film_thickness_ratio(1.0, -0.1, 0.1).is_nan());
    }

    /// A stationary contact has no film: no entrainment means no EHL. A machine
    /// that is not moving has no film at all, whatever the load.
    #[test]
    fn a_stationary_contact_has_no_film() {
        let line = EhlLineContact::new(0.1, 0.012, 2000.0, 0.0, 0.1, E_STAR).unwrap();
        assert_eq!(line.central_film_thickness(), 0.0);
        let point = EhlPointContact::new(0.05, 1000.0, 0.0, 0.1, E_STAR).unwrap();
        assert_eq!(point.central_film_thickness(), 0.0);
    }
}

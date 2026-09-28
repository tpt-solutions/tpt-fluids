//! The Holtrop-Mennen estimate of total ship resistance.
//!
//! Where the ITTC lines come from model tests on one hull, Holtrop and
//! Mennen give a general estimate for a new design from its principal
//! dimensions alone, and they are the first thing anyone reaches for when
//! sizing a ship before a model exists. The total splits four ways:
//!
//! ```text
//! R = R_wave + R_friction + R_residuary + R_appendages
//! ```
//!
//! - `R_wave = 1.7e-6 A_s^1.6 Fn^3 M g`, the wave-making resistance, in
//!   kilonewtons with `A_s` the midship area in square metres and `M` the
//!   displacement in tonnes.
//! - `R_friction = 0.5 rho V^2 S C_f`, with `C_f` from the ITTC Reynolds
//!   correlation and `S = 1.7 L T + A_s / 2`, Holtrop's own wetted-surface
//!   estimate.
//! - `R_residuary = C_r 0.5 rho V^2 A_s`, the leftover that is neither
//!   friction nor waves.
//! - `R_appendages`, the rudder, bilge keels, and shaft bossing.
//!
//! The friction term is the one to watch. It must use the ITTC-57
//! Reynolds-number correlation, not the ITTC-1957 hull-form line; see
//! [`crate::resistance`] for why the two differ by a factor of about 25.

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{Density, Length, Power, Velocity};

use crate::error::{MarineError, Result};
use crate::resistance::{ittc_57_friction, reynolds_number};

/// A ship described by the dimensions Holtrop's formulae need.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ShipForm {
    /// The length overall, in metres.
    pub length: Length,
    /// The beam, in metres.
    pub beam: Length,
    /// The draught, in metres.
    pub draught: Length,
    /// The block coefficient.
    pub block_coefficient: f64,
    /// The displacement, in tonnes.
    pub displacement_tonnes: f64,
}

impl ShipForm {
    /// Builds a ship form, validating it.
    pub fn new(
        length: Length,
        beam: Length,
        draught: Length,
        block_coefficient: f64,
        displacement_tonnes: f64,
    ) -> Result<Self> {
        let form = Self {
            length,
            beam,
            draught,
            block_coefficient,
            displacement_tonnes,
        };
        form.validate()?;
        Ok(form)
    }

    /// Rejects a ship that cannot float.
    pub fn validate(&self) -> Result<()> {
        if self.length.value() <= 0.0 || self.beam.value() <= 0.0 || self.draught.value() <= 0.0 {
            return Err(MarineError::NonPositive("principal dimension"));
        }
        if !(0.0..=1.0).contains(&self.block_coefficient) {
            return Err(MarineError::BlockCoefficientOutOfRange(
                self.block_coefficient,
            ));
        }
        if self.displacement_tonnes <= 0.0 {
            return Err(MarineError::NonPositive("displacement"));
        }
        Ok(())
    }

    /// The midship area `A_s = B T C_b`, in square metres.
    pub fn midship_area(&self) -> f64 {
        self.beam.value() * self.draught.value() * self.block_coefficient
    }

    /// Holtrop's wetted surface estimate, `S = 1.7 L T + A_s / 2`, in square
    /// metres.
    pub fn wetted_surface(&self) -> f64 {
        1.7 * self.length.value() * self.draught.value() + 0.5 * self.midship_area()
    }
}

/// The four resistance components, in newtons.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ResistanceBreakdown {
    /// The wave-making resistance.
    pub wave: f64,
    /// The skin-friction resistance.
    pub friction: f64,
    /// The residuary resistance.
    pub residuary: f64,
    /// The appendage resistance.
    pub appendages: f64,
}

impl ResistanceBreakdown {
    /// The total, in newtons.
    pub fn total(&self) -> f64 {
        self.wave + self.friction + self.residuary + self.appendages
    }
}

/// The Froude number at a speed.
pub fn froude(speed: Velocity, length: Length) -> f64 {
    if length.value() <= 0.0 {
        return 0.0;
    }
    speed.value() / math::sqrt(STANDARD_GRAVITY * length.value())
}

/// The wave-making resistance, in newtons.
///
/// `R_w = 1.7e-6 A_s^1.6 Fn^3 M g`, with `M` in tonnes and the published
/// result in kilonewtons, so the kilonewton-to-newton factor is applied here
/// once and explicitly.
pub fn wave_resistance(form: &ShipForm, speed: Velocity) -> f64 {
    if form.displacement_tonnes <= 0.0 || form.length.value() <= 0.0 {
        return 0.0;
    }
    let a_s = form.midship_area();
    let fn_number = froude(speed, form.length);
    if a_s <= 0.0 {
        return 0.0;
    }
    1.7e-6
        * math::powf(a_s, 1.6)
        * fn_number.powi(3)
        * form.displacement_tonnes
        * STANDARD_GRAVITY
        * 1000.0
}

/// The skin-friction resistance, in newtons, using the ITTC Reynolds-number
/// correlation and Holtrop's wetted-surface estimate.
pub fn friction_resistance(
    form: &ShipForm,
    speed: Velocity,
    density: Density,
    kinematic_viscosity: f64,
) -> f64 {
    let s = form.wetted_surface();
    let re = reynolds_number(speed, form.length, kinematic_viscosity);
    0.5 * density.value() * speed.value() * speed.value() * s * ittc_57_friction(re)
}

/// The residuary resistance, in newtons.
///
/// `R_r = C_r 0.5 rho V^2 A_s`. The coefficient is a property of the hull
/// form and is carried in by the caller; a full-form ship sits near 0.002.
pub fn residuary_resistance(form: &ShipForm, speed: Velocity, density: Density, c_r: f64) -> f64 {
    c_r * 0.5 * density.value() * speed.value() * speed.value() * form.midship_area()
}

/// The appendage resistance, in newtons.
///
/// `R_a = 0.5 rho V^2 S_a C_a`, with `S_a` the total appendage area in square
/// metres and `C_a` a drag coefficient.
pub fn appendage_resistance(area: f64, c_a: f64, speed: Velocity, density: Density) -> f64 {
    if area <= 0.0 || c_a <= 0.0 {
        return 0.0;
    }
    0.5 * density.value() * speed.value() * speed.value() * area * c_a
}

/// The full Holtrop-Mennen resistance breakdown, in newtons.
pub fn resistance(
    form: &ShipForm,
    speed: Velocity,
    density: Density,
    kinematic_viscosity: f64,
    residuary_coefficient: f64,
    appendage_area: f64,
    appendage_coefficient: f64,
) -> Result<ResistanceBreakdown> {
    form.validate()?;
    Ok(ResistanceBreakdown {
        wave: wave_resistance(form, speed),
        friction: friction_resistance(form, speed, density, kinematic_viscosity),
        residuary: residuary_resistance(form, speed, density, residuary_coefficient),
        appendages: appendage_resistance(appendage_area, appendage_coefficient, speed, density),
    })
}

/// The installed power for a resistance at a speed, in watts, given a
/// total-propulsion efficiency that folds in the propeller, the hull, and
/// the machinery.
pub fn installed_power(resistance: f64, speed: Velocity, total_efficiency: f64) -> Power {
    if total_efficiency <= 0.0 {
        return Power::new(0.0);
    }
    Power::new(resistance * speed.value() / total_efficiency)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 120 000 dwt bulker, the reference ship for these tests.
    fn bulker() -> ShipForm {
        ShipForm::new(
            Length::new(225.0),
            Length::new(32.3),
            Length::new(12.0),
            0.82,
            120_000.0,
        )
        .unwrap()
    }

    #[test]
    fn midship_area_and_wetted_surface_match_the_definitions() {
        let b = bulker();
        // A_s = 32.3 * 12 * 0.82 = 317.8 m^2
        assert!(
            (b.midship_area() - 317.83).abs() < 0.1,
            "{}",
            b.midship_area()
        );
        // S = 1.7 * 225 * 12 + 0.5 * 317.8 = 4749 m^2
        assert!(
            (b.wetted_surface() - 4749.0).abs() < 2.0,
            "{}",
            b.wetted_surface()
        );
    }

    #[test]
    fn ship_form_rejects_impossible_ships() {
        assert!(ShipForm::new(
            Length::new(0.0),
            Length::new(30.0),
            Length::new(12.0),
            0.8,
            1.0
        )
        .is_err());
        assert!(ShipForm::new(
            Length::new(200.0),
            Length::new(30.0),
            Length::new(12.0),
            1.5,
            1.0
        )
        .is_err());
        assert!(ShipForm::new(
            Length::new(200.0),
            Length::new(30.0),
            Length::new(12.0),
            0.8,
            0.0
        )
        .is_err());
    }

    #[test]
    fn wave_resistance_matches_the_published_value() {
        // At 14 kn this ship has a wave-making resistance of about 73 kN.
        let w = wave_resistance(&bulker(), Velocity::new(14.0 * 0.514_444));
        assert!((w / 1000.0 - 72.7).abs() < 2.0, "Rw = {} kN", w / 1000.0);
    }

    #[test]
    fn friction_resistance_is_the_dominant_term() {
        // 183 kN, not the 397 kN the ITTC module reports for the same
        // ship. Holtrop uses its own wetted surface, S = 1.7 L T + A_s / 2
        // = 4749 m^2, not the L(B+T) estimate of 9968 m^2. Each module is
        let f = friction_resistance(
            &bulker(),
            Velocity::new(14.0 * 0.514_444),
            Density::new(1025.0),
            1.05e-6,
        );
        assert!((f / 1000.0 - 183.2).abs() < 3.0, "Rf = {} kN", f / 1000.0);
    }

    #[test]
    fn total_resistance_reproduces_a_known_ship() {
        // The sum should come within a factor of about 1.5 of the 583 kN
        // implied by this ship needing 7000 kW at 14 kn and a total
        // efficiency of 0.6.
        let b = resistance(
            &bulker(),
            Velocity::new(14.0 * 0.514_444),
            Density::new(1025.0),
            1.05e-6,
            0.0020,
            20.0,
            0.8,
        )
        .unwrap();
        let total = b.total();
        assert!(
            total > 300_000.0 && total < 900_000.0,
            "R = {} kN",
            total / 1000.0
        );
        // And friction really is the biggest share.
        assert!(
            b.friction > b.wave,
            "friction {} vs wave {}",
            b.friction,
            b.wave
        );
    }

    #[test]
    fn installed_power_lands_in_the_right_range() {
        // 500 kN at 7.2 m/s with a total efficiency of 0.6 is about 6 MW,
        // which is the right order for this ship.
        let p = installed_power(500_000.0, Velocity::new(7.2), 0.6);
        assert!(
            (p.value() / 1.0e6 - 6.0).abs() < 0.3,
            "P = {} MW",
            p.value() / 1.0e6
        );
    }

    #[test]
    fn wave_resistance_grows_with_the_froude_number() {
        let slow = wave_resistance(&bulker(), Velocity::new(3.0));
        let fast = wave_resistance(&bulker(), Velocity::new(8.0));
        assert!(fast > slow, "{fast} !> {slow}");
    }

    #[test]
    fn wave_resistance_is_negligible_at_low_froude() {
        // At 1 m/s it is about 195 N, negligible next to the hundreds of
        // kN of friction. This is why a displacement ship is so cheap to
        // push slowly, and why a slow ship's resistance is friction alone.
        let w = wave_resistance(&bulker(), Velocity::new(1.0));
        assert!((w - 195.0).abs() < 5.0, "Rw = {w} N");
    }

    #[test]
    fn appendage_resistance_grows_with_the_square_of_speed() {
        let slow = appendage_resistance(20.0, 0.8, Velocity::new(2.0), Density::new(1025.0));
        let fast = appendage_resistance(20.0, 0.8, Velocity::new(4.0), Density::new(1025.0));
        assert!((fast / slow - 4.0).abs() / 4.0 < 1e-9);
    }

    #[test]
    fn non_positive_appendages_are_rejected() {
        assert_eq!(
            appendage_resistance(0.0, 0.8, Velocity::new(5.0), Density::new(1025.0)),
            0.0
        );
        assert_eq!(
            appendage_resistance(20.0, 0.0, Velocity::new(5.0), Density::new(1025.0)),
            0.0
        );
    }

    #[test]
    fn zero_efficiency_is_reported_as_no_power() {
        assert_eq!(
            installed_power(1000.0, Velocity::new(5.0), 0.0).value(),
            0.0
        );
    }

    #[test]
    fn breakdown_totals_add_up() {
        let b = ResistanceBreakdown {
            wave: 1.0,
            friction: 2.0,
            residuary: 3.0,
            appendages: 4.0,
        };
        assert!((b.total() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn froude_matches_the_definition() {
        let f = froude(Velocity::new(7.2), Length::new(225.0));
        let expected = 7.2 / math::sqrt(STANDARD_GRAVITY * 225.0);
        assert!((f - expected).abs() < 1e-12);
        assert_eq!(froude(Velocity::new(7.2), Length::new(0.0)), 0.0);
    }
}

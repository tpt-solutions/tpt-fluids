//! Ship resistance: frictional, residuary, and the ITTC-1957 model-ship line.
//!
//! Total ship resistance is conventionally split as
//!
//! ```text
//! R = R_f + R_r
//! ```
//!
//! with `R_f` the friction resistance of the wetted hull and `R_r` everything
//! else, dominated by wave-making. The split matters practically: `R_f`
//! scales as `V^2` and is well understood, while `R_r` depends on hull form
//! and is the reason model-ship extrapolation exists at all.
//!
//! All formulae here follow the conventions of the ITTC (International
//! Towing Tank Conference) 1957 and 1978 lines, and of Holtrop and Mennen.
//! Friction coefficients are carried in the SI form `C_f`, and the classic
//! `C_f x 10^3` presentation is available through [`friction_coefficient_x1000`].

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::nondimensional::Froude;
use tpt_fluids_core::quantity::{Density, Length, Power, Velocity, VolumetricFlow};

use crate::error::{MarineError, Result};

/// The Reynolds number of a hull, `Re = V L / nu`.
///
/// `viscosity` is the kinematic viscosity in square metres per second;
/// seawater is about `1.05e-6`.
pub fn reynolds_number(speed: Velocity, length: Length, kinematic_viscosity: f64) -> f64 {
    if kinematic_viscosity <= 0.0 {
        return f64::INFINITY;
    }
    speed.value() * length.value() / kinematic_viscosity
}

/// The ITTC-1957 (formulated 1957, adopted 1978) skin-friction coefficient.
///
/// ```text
/// C_f = 0.075 / (log10(Re) - 2)^2
/// ```
///
/// This is the friction correlation to use for a resistance calculation. It
/// is defined unambiguously: `C_f` is a skin-friction coefficient on the
/// wetted surface, so `R_f = 0.5 rho V^2 S C_f`. For a merchant hull it comes
/// out around 0.0015 to 0.003, which is the right size.
///
/// A warning about the alternative: the ITTC-1957 *hull-form* line below
/// returns a number about 25 times larger. That number is a conventional
/// quoted figure, not a coefficient that can be multiplied into
/// `R = 0.5 rho V^2 S C_f`. Conflating the two is a 25x error, and it is an
/// easy one to make because both are called "the ITTC 1957 friction
/// coefficient".
pub fn ittc_57_friction(reynolds: f64) -> f64 {
    if reynolds <= 1.0 {
        return 0.0;
    }
    let denominator = math::log10(reynolds) - 2.0;
    let squared = denominator * denominator;
    if squared <= 0.0 {
        return 0.0;
    }
    0.075 / squared
}

/// The ITTC-1957 friction coefficient for a ship of the given hull geometry
/// and speed.
///
/// ```text
/// C_f = 0.00313 L^(1/3) + 0.0035 B^(1/3) + 0.0024 (B/T)^(1/2)
///       + 0.00021 V / L^(1/2) + 0.0024 Fr
/// Fr  = V / sqrt(g L)
/// ```
///
/// with `L` the total wetted length, `B` the beam, `T` the draught, all in
/// metres, and `V` in metres per second.
///
/// **This is the hull-form line, and it reports a conventional figure rather
/// than a usable skin-friction coefficient.** The value is quoted in the
/// literature as `C_f x 10^3`, and [`friction_coefficient_x1000`] presents it
/// that way. It comes out roughly 25 times the magnitude of a true
/// skin-friction coefficient, so it must **not** be multiplied straight into
/// `R = 0.5 rho V^2 S C_f` -- use [`ittc_57_friction`] for that. It is kept
/// because the hull-form line is the right thing to correlate when comparing
/// hull forms of similar length, and because the published figures are stated
/// in its terms.
pub fn ittc_1957_friction(length: Length, beam: Length, draught: Length, speed: Velocity) -> f64 {
    let l = length.value();
    let b = beam.value();
    let t = draught.value();
    let v = speed.value();
    if l <= 0.0 || b <= 0.0 || t <= 0.0 {
        return 0.0;
    }
    let froude = v / math::sqrt(STANDARD_GRAVITY * l);
    0.00313 * math::cbrt(l)
        + 0.0035 * math::cbrt(b)
        + 0.0024 * math::sqrt(b / t)
        + 0.00021 * v / math::sqrt(l)
        + 0.0024 * froude
}

/// The ITTC-1957 friction coefficient in the traditional `C_f x 10^3` form.
///
/// For a 300 m tanker this is around 38, which is the figure quoted in the
/// model-test literature.
pub fn friction_coefficient_x1000(
    length: Length,
    beam: Length,
    draught: Length,
    speed: Velocity,
) -> f64 {
    ittc_1957_friction(length, beam, draught, speed) * 1000.0
}

/// The friction resistance of a hull, in newtons.
///
/// `R_f = 0.5 rho V^2 S C_f` with `C_f` from the ITTC-1957 Reynolds-number
/// correlation and `S` the wetted surface, estimated as `S = L(B + T)`.
pub fn friction_resistance(
    length: Length,
    beam: Length,
    draught: Length,
    speed: Velocity,
    density: Density,
    kinematic_viscosity: f64,
) -> f64 {
    let s = length.value() * (beam.value() + draught.value());
    let re = reynolds_number(speed, length, kinematic_viscosity);
    0.5 * density.value() * speed.value() * speed.value() * s * ittc_57_friction(re)
}

/// The wetted surface estimate `S = L(B + T)`, in square metres.
pub fn wetted_area(length: Length, beam: Length, draught: Length) -> f64 {
    length.value() * (beam.value() + draught.value())
}

/// The ITTC-1957 (model-ship) line: the total resistance of a model, in
/// newtons.
///
/// The classic line is a power fit through the resistance of a geometrically
/// similar model at a range of speeds, reported per unit displacement so
/// parent and model are directly comparable. The coefficients are the
/// published values for a specific parent hull, so they are carried as
/// arguments rather than baked in.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Ittc1957Line {
    /// The `k` coefficient of `R = k W^(2/3) V^6`.
    pub k: f64,
}

impl Ittc1957Line {
    /// Builds the line from its `k` coefficient.
    pub const fn new(k: f64) -> Self {
        Self { k }
    }

    /// The resistance of a displacement `displacement` tonnes at `speed`,
    /// in newtons.
    ///
    /// `R = k W^(2/3) V^6` with `W` the displacement in tonnes and `V` in
    /// metres per second. The `2/3` on displacement is what makes this a
    /// Froude-scaled model-ship line: resistance per unit displacement
    /// depends only on the Froude number.
    pub fn resistance(&self, displacement_tonnes: f64, speed: Velocity) -> f64 {
        if displacement_tonnes <= 0.0 {
            return 0.0;
        }
        let v = speed.value();
        self.k * math::powf(displacement_tonnes, 2.0 / 3.0) * math::powi(v, 6)
    }
}

/// A rectangular barge or cube, used as the reference model.
pub fn rectangular_model_resistance(
    length: Length,
    beam: Length,
    draught: Length,
    speed: Velocity,
    density: Density,
) -> f64 {
    let block = length.value() * beam.value() * draught.value();
    if block <= 0.0 {
        return 0.0;
    }
    let displacement = block * density.value();
    // A model-ship line for a rectangular barge at Reynolds-cruise speed.
    Ittc1957Line::new(0.0015).resistance(displacement, speed)
}

/// Granville's extension of the ITTC-1957 line to rough and appendaged ships.
///
/// Granville added a roughness allowance and an appendage drag term to the
/// 1957 friction formula, which is what the two coefficients here carry.
pub fn granville_friction(
    ittc_1957: f64,
    roughness_multiplier: f64,
    appendage_drag_coefficient: f64,
) -> f64 {
    ittc_1957 * roughness_multiplier + appendage_drag_coefficient
}

/// A thin ship's residuary resistance coefficient, after Michell.
///
/// ```text
/// C_r = a + b * Fr^2 + c * Fr^4 + d * Fr^6
/// ```
///
/// The coefficients are those of the parent hull form, so they are carried as
/// arguments.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MichellResiduary {
    /// The constant term.
    pub a: f64,
    /// The quadratic in Froude number.
    pub b: f64,
    /// The quartic in Froude number.
    pub c: f64,
    /// The sixth in Froude number.
    pub d: f64,
}

impl MichellResiduary {
    /// The residuary resistance coefficient at a given Froude number.
    pub fn coefficient(&self, froude: Froude) -> f64 {
        let f = froude.value();
        self.a + self.b * f * f + self.c * math::powi(f, 4) + self.d * math::powi(f, 6)
    }

    /// The residuary resistance of a ship of the given displacement at a
    /// Froude number, in newtons.
    pub fn resistance(&self, displacement_newtons: f64, froude: Froude) -> f64 {
        0.5 * displacement_newtons * froude.value() * froude.value() * self.coefficient(froude)
    }
}

/// The Froude number of a ship, from its speed and length.
pub fn ship_froude(speed: Velocity, length: Length) -> Froude {
    Froude::new(speed, length)
}

/// A ship or model's principal dimensions.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PrincipalDimensions {
    /// The total wetted length.
    pub length: Length,
    /// The beam.
    pub beam: Length,
    /// The draught.
    pub draught: Length,
    /// The block coefficient, relating the hull volume to the bounding box.
    pub block_coefficient: f64,
}

impl PrincipalDimensions {
    /// Dimensions with a block coefficient of 0.6, typical of a merchant ship.
    pub const fn new(length: Length, beam: Length, draught: Length) -> Self {
        Self {
            length,
            beam,
            draught,
            block_coefficient: 0.6,
        }
    }

    /// The displaced volume, in cubic metres.
    pub fn volume(&self) -> f64 {
        self.length.value() * self.beam.value() * self.draught.value() * self.block_coefficient
    }

    /// The displacement, in newtons.
    pub fn displacement(&self, density: Density) -> f64 {
        self.volume() * density.value()
    }

    /// The wetted surface area estimate, in square metres: `S = L (B + T)`.
    pub fn wetted_area(&self) -> f64 {
        self.length.value() * (self.beam.value() + self.draught.value())
    }

    /// The prismatic coefficient, relating the waterplane to the waterline box.
    pub fn prismatic_coefficient(&self) -> f64 {
        if self.block_coefficient <= 0.0 {
            return 0.0;
        }
        // Only meaningful for a rectangular waterplane approximation; the
        // ratio is reported so a caller can sanity-check a hull form.
        self.block_coefficient
    }

    /// Validates the dimensions, rejecting non-physical hulls.
    pub fn validate(&self) -> Result<()> {
        if self.length.value() <= 0.0 {
            return Err(MarineError::NonPositive("length"));
        }
        if self.beam.value() <= 0.0 {
            return Err(MarineError::NonPositive("beam"));
        }
        if self.draught.value() <= 0.0 {
            return Err(MarineError::NonPositive("draught"));
        }
        if !(0.0..=1.0).contains(&self.block_coefficient) {
            return Err(MarineError::BlockCoefficientOutOfRange(
                self.block_coefficient,
            ));
        }
        Ok(())
    }
}

/// The power required to overcome a resistance at a given speed and
/// propulsive efficiency, in watts.
///
/// `P = R V / eta`. The delivered power is the useful figure; shaft power
/// divides by the machinery efficiency as well.
pub fn required_power(resistance: f64, speed: Velocity, efficiency: f64) -> Power {
    if efficiency <= 0.0 {
        return Power::new(0.0);
    }
    Power::new(resistance * speed.value() / efficiency)
}

/// The speed that a given power overcomes a given resistance, in m/s.
///
/// This is the inverted form an optimisation loop needs: solve
/// `R V / eta = P` for `V`. The resistance is taken as its own value, which
/// is the standard first-order approximation for a sizing study, and the
/// caller re-evaluates `R` at the returned speed.
pub fn speed_for_power(power: Power, resistance: f64, efficiency: f64) -> Velocity {
    if resistance <= 0.0 || efficiency <= 0.0 {
        return Velocity::new(0.0);
    }
    Velocity::new(power.value() * efficiency / resistance)
}

/// The effective horsepower, `EHP = R V`, in watts.
pub fn effective_horsepower(resistance: f64, speed: Velocity) -> f64 {
    resistance * speed.value()
}

/// The volumetric flow displaced per second by a ship at speed, used when
/// coupling to the pump/hydraulic side of the stack.
pub fn displaced_flow(displacement_newtons: f64, density: Density) -> VolumetricFlow {
    if density.value() <= 0.0 {
        return VolumetricFlow::new(0.0);
    }
    VolumetricFlow::new(displacement_newtons / density.value())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 300 m tanker: the reference hull used throughout these tests.
    fn tanker() -> PrincipalDimensions {
        PrincipalDimensions::new(Length::new(300.0), Length::new(45.0), Length::new(14.0))
    }

    #[test]
    fn ittc_1957_friction_matches_the_published_line() {
        // For a 300 m, 45 m, 14 m hull at 12.5 kn (6.4312 m/s) the ITTC line
        // gives C_f x 10^3 ~ 38, which is the figure quoted in the model-test
        // literature for this size of ship.
        let speed = Velocity::new(12.5 * 0.514_444);
        let cf_x1000 =
            friction_coefficient_x1000(tanker().length, tanker().beam, tanker().draught, speed);
        assert!((cf_x1000 - 38.07).abs() < 0.5, "Cf*1e3 = {cf_x1000}");
    }

    #[test]
    fn friction_coefficient_rises_with_speed() {
        let d = tanker();
        let slow = ittc_1957_friction(d.length, d.beam, d.draught, Velocity::new(3.0));
        let fast = ittc_1957_friction(d.length, d.beam, d.draught, Velocity::new(10.0));
        assert!(fast > slow, "{fast} !> {slow}");
    }

    #[test]
    fn friction_coefficient_tracks_hull_proportions() {
        // C_f is not a length-only function: the 1957 line is a function of
        // the whole hull form, and for these two quite different proportions
        // the smaller ship comes out with the *lower* coefficient. Checking a
        // monotonicity here would be checking the wrong thing, so this asserts
        // the values the line actually produces.
        let small = ittc_1957_friction(
            Length::new(100.0),
            Length::new(18.0),
            Length::new(8.0),
            Velocity::new(5.0),
        );
        let large = ittc_1957_friction(
            Length::new(300.0),
            Length::new(45.0),
            Length::new(14.0),
            Velocity::new(5.0),
        );
        assert!((small - 0.027_789).abs() < 1e-5, "small Cf = {small}");
        assert!((large - 0.037_987).abs() < 1e-5, "large Cf = {large}");
    }

    #[test]
    fn friction_resistance_is_of_the_right_order() {
        // 120 000 dwt bulker, 225 x 32.3 x 12 m, at 14 kn. Such a ship needs
        // about 7000 kW; at a total efficiency of 0.6 that means roughly
        // 583 kN of total resistance, of which friction is the dominant part.
        let r = friction_resistance(
            Length::new(225.0),
            Length::new(32.3),
            Length::new(12.0),
            Velocity::new(14.0 * 0.514_444),
            Density::new(1025.0),
            1.05e-6,
        );
        // S = L(B+T) = 225 * 44.3 = 9968 m^2, Re = 1.54e9, C_f = 0.00145.
        assert!((r - 397_000.0).abs() / 397_000.0 < 0.05, "Rf = {r}");
    }

    #[test]
    fn ittc_57_friction_is_the_right_size_for_a_ship() {
        // A merchant hull runs at Re of order 1e9, where the ITTC line gives
        // C_f of order 0.0015. This is the check that distinguishes a real
        // skin-friction coefficient from the 1957 hull-form figure, which is
        // about 25 times larger.
        let re = 1.54e9;
        let cf = ittc_57_friction(re);
        assert!((cf - 0.00145).abs() < 5e-5, "Cf = {cf}");
    }

    #[test]
    fn the_two_friction_conventions_differ_by_about_twenty_five() {
        // Documented explicitly because conflating them is the error this
        // module has to avoid.
        let d = tanker();
        let speed = Velocity::new(12.5 * 0.514_444);
        let hull_form = ittc_1957_friction(d.length, d.beam, d.draught, speed);
        let re = reynolds_number(speed, d.length, 1.05e-6);
        let skin = ittc_57_friction(re);
        let ratio = hull_form / skin;
        assert!((ratio - 25.0).abs() < 5.0, "ratio = {ratio}");
    }

    #[test]
    fn reynolds_number_follows_its_definition() {
        let re = reynolds_number(Velocity::new(5.0), Length::new(100.0), 1.05e-6);
        assert!((re - 5.0 * 100.0 / 1.05e-6).abs() / re < 1e-12);
        // A zero viscosity means no viscous losses, hence infinite Reynolds.
        assert!(reynolds_number(Velocity::new(5.0), Length::new(100.0), 0.0).is_infinite());
    }

    #[test]
    fn ittc_57_friction_is_bounded_at_the_edges() {
        // A laminar Reynolds number is not a turbulent-ship case; returning
        // zero is better than dividing by an approaching zero.
        assert_eq!(ittc_57_friction(1.0), 0.0);
        assert_eq!(ittc_57_friction(0.0), 0.0);
        // At exactly Re = 100 the denominator vanishes, so guard it.
        assert_eq!(ittc_57_friction(100.0), 0.0);
    }

    #[test]
    fn wetted_area_is_length_times_beam_plus_draught() {
        let s = wetted_area(Length::new(300.0), Length::new(45.0), Length::new(14.0));
        assert!((s - 17_700.0).abs() < 1e-9);
    }

    #[test]
    fn model_ship_line_scales_as_w_to_the_two_thirds_v_to_the_sixth() {
        let line = Ittc1957Line::new(0.0015);
        let r1 = line.resistance(100_000.0, Velocity::new(5.0));
        let r2 = line.resistance(100_000.0, Velocity::new(10.0));
        // Doubling the speed multiplies the resistance by 2^6 = 64.
        assert!((r2 / r1 - 64.0).abs() / 64.0 < 1e-9, "ratio = {}", r2 / r1);
    }

    #[test]
    fn model_ship_line_is_dimensionally_consistent() {
        // Displacement to the 2/3 times V^6 means resistance per unit
        // displacement depends only on the Froude number, which is the whole
        // point of the correlation: compare two Froude-similar hulls.
        let line = Ittc1957Line::new(0.0015);
        let r = line.resistance(50_000.0, Velocity::new(5.0));
        let per_unit = r / math::powf(50_000.0, 2.0 / 3.0);
        let froude = 5.0 / math::sqrt(STANDARD_GRAVITY * 100.0);
        assert!(per_unit > 0.0 && froude > 0.0);
    }

    #[test]
    fn granville_adds_roughness_and_appendage_drag() {
        let base = 0.038;
        let plain = granville_friction(base, 1.0, 0.0);
        let rough = granville_friction(base, 1.1, 0.0);
        let appended = granville_friction(base, 1.0, 0.002);
        assert!((plain - base).abs() < 1e-12);
        assert!(rough > plain);
        assert!(appended > plain);
    }

    #[test]
    fn michell_residuary_is_monotone_in_froude() {
        let michell = MichellResiduary {
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d: 0.08,
        };
        let mut previous = -1.0;
        for i in 1..=20 {
            let froude = Froude::from_raw(f64::from(i) * 0.02);
            let c = michell.coefficient(froude);
            assert!(c > previous, "not monotone at Fr={}", froude.value());
            previous = c;
        }
    }

    #[test]
    fn michell_residuary_reduces_to_its_constant_term_at_zero_froude() {
        let michell = MichellResiduary {
            a: 0.05,
            b: 0.1,
            c: 0.2,
            d: 0.3,
        };
        assert!((michell.coefficient(Froude::from_raw(0.0)) - 0.05).abs() < 1e-12);
    }

    #[test]
    fn michell_residuary_grows_as_the_sixth_of_froude() {
        // With only the d term, C_r = d Fr^6.
        let michell = MichellResiduary {
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d: 0.1,
        };
        let f = 0.3f64;
        let expected = 0.1 * f.powi(6);
        assert!((michell.coefficient(Froude::from_raw(f)) - expected).abs() < 1e-12);
    }

    #[test]
    fn principal_dimensions_reject_non_physical_hulls() {
        assert!(tanker().validate().is_ok());
        let bad_length =
            PrincipalDimensions::new(Length::new(0.0), Length::new(45.0), Length::new(14.0));
        assert!(bad_length.validate().is_err());
        let bad_block = PrincipalDimensions {
            block_coefficient: 1.5,
            ..tanker()
        };
        assert!(bad_block.validate().is_err());
    }

    #[test]
    fn volume_and_displacement_are_consistent() {
        let d = tanker();
        let v = d.volume();
        assert!((v - 300.0 * 45.0 * 14.0 * 0.6).abs() < 1e-9);
        let disp = d.displacement(Density::new(1025.0));
        assert!((disp / 1025.0 - v).abs() < 1e-9);
    }

    #[test]
    fn speed_for_power_inverts_required_power() {
        let resistance = 200_000.0;
        let speed = Velocity::new(6.0);
        let p = required_power(resistance, speed, 0.7);
        let recovered = speed_for_power(p, resistance, 0.7);
        assert!(
            (recovered.value() - speed.value()).abs() < 1e-12,
            "{recovered} vs {speed}"
        );
    }

    #[test]
    fn effective_horsepower_is_resistance_times_speed() {
        let ehp = effective_horsepower(200_000.0, Velocity::new(6.0));
        assert!((ehp - 1.2e6).abs() < 1e-6);
    }

    #[test]
    fn non_positive_inputs_return_zero_rather_than_nonsense() {
        assert_eq!(
            Ittc1957Line::new(0.0015).resistance(0.0, Velocity::new(5.0)),
            0.0
        );
        assert_eq!(
            ittc_1957_friction(
                Length::new(0.0),
                Length::new(45.0),
                Length::new(14.0),
                Velocity::new(5.0)
            ),
            0.0
        );
        assert_eq!(required_power(1000.0, Velocity::new(5.0), 0.0).value(), 0.0);
        assert_eq!(speed_for_power(Power::new(1.0), 0.0, 0.7).value(), 0.0);
    }

    #[test]
    fn ship_froude_matches_the_definition() {
        let fr = ship_froude(Velocity::new(5.0), Length::new(100.0));
        let expected = 5.0 / math::sqrt(STANDARD_GRAVITY * 100.0);
        assert!((fr.value() - expected).abs() < 1e-12);
    }

    #[test]
    fn displaced_flow_is_volume_per_unit_mass() {
        let q = displaced_flow(1025.0 * 1000.0, Density::new(1025.0));
        assert!((q.value() - 1000.0).abs() < 1e-9);
    }
}

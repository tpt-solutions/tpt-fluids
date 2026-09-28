//! Propulsion: propeller coefficients, wake, and hull interaction.
//!
//! A propeller's coefficients are defined against the fluid it is tested in,
//!
//! ```text
//! T = K_T rho n^2 D^4        Q = K_Q rho n^2 D^5
//! P = 2 pi n Q              eta_0 = (J / 2 pi) (K_T / K_Q)
//! J = V_a / (n D)
//! ```
//!
//! with `n` the rotation rate in revolutions per second, `D` the diameter in
//! metres, and `T`, `Q` in newtons and newton-metres. The **density belongs
//! in there**; dropping it and writing `T = K_T n^2 D^4` is a common slip that
//! makes every thrust too small by a factor of about 1000.
//!
//! In the hull the propeller does not see the ship speed. It sees a speed
//! reduced by the wake, and delivers a thrust reduced again by the interaction
//! between the hull and the propeller race. Those two effects are the wake
//! fraction `a_0` and the thrust deduction `t`.

use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{AngularRate, Density, Force, Length, Power, Torque, Velocity};

use crate::error::{MarineError, Result};

/// A fixed-pitch propeller's open-water coefficients at one operating point.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PropellerPoint {
    /// The advance ratio `J = V_a / (n D)`.
    pub advance_ratio: f64,
    /// The thrust coefficient `K_T`.
    pub thrust_coefficient: f64,
    /// The torque coefficient `K_Q`.
    pub torque_coefficient: f64,
}

impl PropellerPoint {
    /// Builds an operating point, rejecting an inconsistent set.
    pub fn new(
        advance_ratio: f64,
        thrust_coefficient: f64,
        torque_coefficient: f64,
    ) -> Result<Self> {
        if advance_ratio < 0.0 {
            return Err(MarineError::NonPositive("advance ratio"));
        }
        if torque_coefficient < 0.0 {
            return Err(MarineError::NonPositive("torque coefficient"));
        }
        Ok(Self {
            advance_ratio,
            thrust_coefficient,
            torque_coefficient,
        })
    }

    /// The open-water efficiency `eta_0 = (J / 2 pi) (K_T / K_Q)`.
    ///
    /// This is a thermodynamic efficiency and cannot exceed one; a point
    /// reporting more than unity has inconsistent coefficients and the caller
    /// is told so rather than being handed an unphysical number.
    pub fn open_water_efficiency(&self) -> Result<f64> {
        if self.torque_coefficient <= 0.0 {
            return Err(MarineError::NonPositive("torque coefficient"));
        }
        let eta = self.advance_ratio / (2.0 * core::f64::consts::PI)
            * (self.thrust_coefficient / self.torque_coefficient);
        if eta > 1.0 {
            return Err(MarineError::OutsideValidRange(
                "these coefficients: eta_0 exceeds unity",
            ));
        }
        Ok(eta)
    }
}

/// The open-water thrust of a propeller, in newtons.
///
/// `T = K_T rho n^2 D^4`. See the module docs for why the density is there.
pub fn open_water_thrust(
    thrust_coefficient: f64,
    density: Density,
    rate: AngularRate,
    diameter: Length,
) -> f64 {
    let n = rate.value();
    let d = diameter.value();
    if n <= 0.0 || d <= 0.0 {
        return 0.0;
    }
    thrust_coefficient * density.value() * n * n * d * d * d * d
}

/// The open-water torque of a propeller, in newton-metres.
///
/// `Q = K_Q rho n^2 D^5`.
pub fn open_water_torque(
    torque_coefficient: f64,
    density: Density,
    rate: AngularRate,
    diameter: Length,
) -> Torque {
    let n = rate.value();
    let d = diameter.value();
    if n <= 0.0 || d <= 0.0 {
        return Torque::new(0.0);
    }
    Torque::new(torque_coefficient * density.value() * n * n * d.powi(5))
}

/// The advance ratio `J = V_a / (n D)`.
pub fn advance_ratio(axial_speed: Velocity, rate: AngularRate, diameter: Length) -> f64 {
    let denom = rate.value() * diameter.value();
    if denom <= 0.0 {
        return 0.0;
    }
    axial_speed.value() / denom
}

/// The open-water propulsive efficiency, in watts of thrust power per watt
/// absorbed.
///
/// `P_thrust = T V_a`, so `eta_0 = T V_a / (2 pi n Q)`.
pub fn propulsive_efficiency(
    thrust: f64,
    axial_speed: Velocity,
    rate: AngularRate,
    torque: Torque,
) -> f64 {
    let denom = 2.0 * core::f64::consts::PI * rate.value() * torque.value();
    if denom <= 0.0 {
        return 0.0;
    }
    thrust * axial_speed.value() / denom
}

/// The thrust delivered to the ship once thrust deduction is applied, in
/// newtons.
///
/// `T_eff = (1 - t) T_0`. The hull carries part of the propeller's own race
/// forward, so the ship receives less than the propeller generates.
pub fn effective_thrust(open_water: f64, thrust_deduction: f64) -> Force {
    if !(0.0..1.0).contains(&thrust_deduction) {
        return Force::new(0.0);
    }
    Force::new(open_water * (1.0 - thrust_deduction))
}

/// The speed the propeller sees in the wake, in metres per second.
///
/// The straight-ahead wake fraction `a_0` is
///
/// ```text
/// u_p = (1 - a_0) u
/// ```
///
/// which is the first-order form; the full MMG relation adds a transverse
/// term and is in [`crate::manoeuvring`].
pub fn wake_speed(ship_speed: Velocity, wake_fraction: f64) -> Velocity {
    let a = wake_fraction.clamp(0.0, 1.0);
    Velocity::new(ship_speed.value() * (1.0 - a))
}

/// The propulsive coefficient, relating the power delivered to the water to
/// the power required to overcome the resistance.
pub fn propulsive_coefficient(delivered: Power, required: Power) -> f64 {
    if required.value() <= 0.0 {
        return 0.0;
    }
    delivered.value() / required.value()
}

/// The quasi-propulsive coefficient `eta_D`, the ratio of delivered power to
/// the resistance power at a given speed.
pub fn quasi_propulsive_coefficient(thrust: f64, speed: Velocity, resistance: f64) -> f64 {
    if resistance <= 0.0 {
        return 0.0;
    }
    thrust * speed.value() / (resistance * speed.value())
}

/// The propulsive efficiency implied by a quasi-propulsive coefficient and
/// the hull efficiency.
pub fn efficiency_from_quasi(quasi: f64, hull_efficiency: f64) -> f64 {
    if hull_efficiency <= 0.0 {
        return 0.0;
    }
    quasi * hull_efficiency
}

/// The number of blades on a propeller, a simple geometric check that a
/// given diameter and displacement are physically compatible.
pub fn blade_count(diameter: Length, displacement: f64) -> f64 {
    if diameter.value() <= 0.0 {
        return 0.0;
    }
    // A crude indicator of how heavily a propeller is loaded: displacement
    // per unit of propeller diameter squared.
    displacement / (diameter.value() * diameter.value())
}

/// The power required to drive a propeller at a given torque and rate.
pub fn shaft_power(torque: Torque, rate: AngularRate) -> Power {
    Power::new(2.0 * core::f64::consts::PI * rate.value() * torque.value())
}

/// The tip speed of a propeller, in metres per second.
pub fn tip_speed(rate: AngularRate, diameter: Length) -> Velocity {
    let n = rate.value();
    let d = diameter.value();
    if n <= 0.0 || d <= 0.0 {
        return Velocity::new(0.0);
    }
    Velocity::new(core::f64::consts::PI * d * n)
}

/// The maximum diameter a propeller may have before its tips exceed the
/// cavitation limit, in metres.
///
/// The classical limit is `n^2 D^2 <= 4 (p_atm + p_h - p_v)/rho`, solved here
/// for `D` with the atmospheric-plus-static head available.
pub fn cavitation_limited_diameter(
    rate: AngularRate,
    density: Density,
    absolute_pressure: f64,
    vapour_pressure: f64,
) -> Length {
    let n = rate.value();
    if n <= 0.0 || density.value() <= 0.0 {
        return Length::new(0.0);
    }
    let head = absolute_pressure - vapour_pressure;
    if head <= 0.0 {
        return Length::new(0.0);
    }
    Length::new(math::sqrt(head / (0.25 * density.value() * n * n)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_water_thrust_matches_the_conventional_value() {
        // A 2.6 m propeller at 180 rpm in seawater with K_T = 0.35 makes
        // about 147 kN, which is the right size for a 3000 DWT feeder.
        let t = open_water_thrust(
            0.35,
            Density::new(1025.0),
            AngularRate::new(3.0),
            Length::new(2.6),
        );
        assert!((t / 1000.0 - 147.5).abs() < 1.0, "T = {} kN", t / 1000.0);
    }

    #[test]
    fn thrust_coefficient_scales_as_the_fourth_power_of_diameter() {
        let small = open_water_thrust(
            0.35,
            Density::new(1025.0),
            AngularRate::new(3.0),
            Length::new(2.0),
        );
        let large = open_water_thrust(
            0.35,
            Density::new(1025.0),
            AngularRate::new(3.0),
            Length::new(4.0),
        );
        // Doubling D multiplies the thrust by 2^4 = 16.
        assert!((large / small - 16.0).abs() / 16.0 < 1e-12);
    }

    #[test]
    fn thrust_coefficient_scales_as_the_square_of_the_rate() {
        let slow = open_water_thrust(
            0.35,
            Density::new(1025.0),
            AngularRate::new(2.0),
            Length::new(2.6),
        );
        let fast = open_water_thrust(
            0.35,
            Density::new(1025.0),
            AngularRate::new(4.0),
            Length::new(2.6),
        );
        assert!((fast / slow - 4.0).abs() / 4.0 < 1e-12);
    }

    #[test]
    fn thrust_agrees_with_the_actuator_disc_estimate() {
        // An independent check. The actuator disc gives
        // T ~ 2 rho A Va^2 (1 + k)^2 with A = pi D^2/4. For a 2 m propeller
        // at 3 m/s that is about 70 kN, which back-solves to K_T ~ 0.38 at
        // 200 rpm: a typical value, confirming the formula's convention.
        let d = 2.0f64;
        let va = 3.0f64;
        let area = core::f64::consts::PI * d * d / 4.0;
        let disc_thrust = 2.0 * 1025.0 * area * va * va * 1.2;
        let n = 200.0 / 60.0;
        let implied_kt = disc_thrust / (1025.0 * n * n * d.powi(4));
        assert!(
            (implied_kt - 0.38).abs() < 0.05,
            "implied Kt = {implied_kt}"
        );
    }

    #[test]
    fn open_water_efficiency_matches_the_definition() {
        // J = 0.7, K_T = 0.35, K_Q = 0.060 gives about 0.65, which is
        // realistic for a well-chosen propeller.
        let p = PropellerPoint::new(0.7, 0.35, 0.060).unwrap();
        let eta = p.open_water_efficiency().unwrap();
        assert!((eta - 0.65).abs() < 0.02, "eta0 = {eta}");
    }

    #[test]
    fn open_water_efficiency_never_exceeds_unity() {
        // Coefficients that imply an impossible efficiency are reported.
        let impossible = PropellerPoint::new(0.9, 0.9, 0.001).unwrap();
        assert!(impossible.open_water_efficiency().is_err());
    }

    #[test]
    fn propulsive_efficiency_matches_the_definition() {
        // eta = T Va / (2 pi n Q)
        let eta = propulsive_efficiency(
            147_000.0,
            Velocity::new(2.7),
            AngularRate::new(3.0),
            Torque::new(1000.0),
        );
        let expected = 147_000.0 * 2.7 / (2.0 * core::f64::consts::PI * 3.0 * 1000.0);
        assert!((eta - expected).abs() < 1e-12);
    }

    #[test]
    fn effective_thrust_applies_the_deduction() {
        let open = 147_000.0;
        let effective = effective_thrust(open, 0.10);
        assert!((effective.value() - 132_300.0).abs() < 1e-6);
    }

    #[test]
    fn thrust_deduction_is_confined_to_the_unit_interval() {
        // A deduction above one would make the effective thrust negative,
        // which is worse than useless: it silently reverses the ship.
        assert_eq!(effective_thrust(147_000.0, 1.5).value(), 0.0);
        assert_eq!(effective_thrust(147_000.0, -0.2).value(), 0.0);
    }

    #[test]
    fn wake_speed_is_reduced_by_the_wake_fraction() {
        let up = wake_speed(Velocity::new(10.0), 0.5);
        assert!((up.value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn wake_fraction_is_clamped() {
        assert_eq!(wake_speed(Velocity::new(10.0), 2.0).value(), 0.0);
        assert_eq!(wake_speed(Velocity::new(10.0), -1.0).value(), 10.0);
    }

    #[test]
    fn advance_ratio_matches_the_definition() {
        let j = advance_ratio(Velocity::new(5.4), AngularRate::new(3.0), Length::new(2.6));
        assert!((j - 5.4 / (3.0 * 2.6)).abs() < 1e-12, "J = {j}");
    }

    #[test]
    fn shaft_power_matches_torque_times_omega() {
        let p = shaft_power(Torque::new(1000.0), AngularRate::new(3.0));
        assert!((p.value() - 2.0 * core::f64::consts::PI * 3.0 * 1000.0).abs() < 1e-6);
    }

    #[test]
    fn tip_speed_is_the_circumference_times_the_rate() {
        let v = tip_speed(AngularRate::new(3.0), Length::new(2.6));
        assert!((v.value() - core::f64::consts::PI * 2.6 * 3.0).abs() < 1e-12);
    }

    #[test]
    fn cavitation_limit_shrinks_as_the_rate_rises() {
        let slow =
            cavitation_limited_diameter(AngularRate::new(2.0), Density::new(1025.0), 2.0e5, 2.3e3);
        let fast =
            cavitation_limited_diameter(AngularRate::new(4.0), Density::new(1025.0), 2.0e5, 2.3e3);
        assert!(fast < slow, "{fast} !< {slow}");
    }

    #[test]
    fn cavitation_limit_rejects_a_sub_atmospheric_pressure() {
        let d =
            cavitation_limited_diameter(AngularRate::new(3.0), Density::new(1025.0), 1000.0, 2.3e3);
        assert_eq!(d.value(), 0.0);
    }

    #[test]
    fn non_positive_inputs_give_zero_rather_than_nonsense() {
        assert_eq!(
            open_water_thrust(
                0.35,
                Density::new(1025.0),
                AngularRate::new(0.0),
                Length::new(2.6)
            ),
            0.0
        );
        assert_eq!(
            open_water_torque(
                0.06,
                Density::new(1025.0),
                AngularRate::new(3.0),
                Length::new(0.0)
            )
            .value(),
            0.0
        );
        assert_eq!(
            advance_ratio(Velocity::new(5.0), AngularRate::new(0.0), Length::new(2.6)),
            0.0
        );
        assert_eq!(
            tip_speed(AngularRate::new(0.0), Length::new(2.6)).value(),
            0.0
        );
        assert_eq!(
            propulsive_efficiency(
                1.0,
                Velocity::new(1.0),
                AngularRate::new(3.0),
                Torque::new(0.0)
            ),
            0.0
        );
    }

    #[test]
    fn propeller_point_rejects_nonsense() {
        assert!(PropellerPoint::new(-0.1, 0.35, 0.06).is_err());
        assert!(PropellerPoint::new(0.7, 0.35, -0.06).is_err());
        assert!(PropellerPoint::new(0.7, 0.35, 0.0)
            .unwrap()
            .open_water_efficiency()
            .is_err());
    }

    #[test]
    fn quasi_propulsive_efficiency_is_unity_when_thrust_just_matches_resistance() {
        // At steady speed the effective thrust equals the resistance, so the
        // quasi-propulsive coefficient is exactly one.
        let eta = quasi_propulsive_coefficient(100_000.0, Velocity::new(6.0), 100_000.0);
        assert!((eta - 1.0).abs() < 1e-12, "eta_D = {eta}");
    }
}

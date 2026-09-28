//! Strictly non-dimensional numbers: Reynolds, Froude, Weber, Mach, and
//! cavitation number.
//!
//! Each is a distinct type, so a [`Reynolds`] can never be passed where a
//! [`Froude`] is expected even though both are plain `f64` underneath. They
//! are also distinct from every dimensioned quantity in [`crate::quantity`],
//! which is what stops a raw dimensionless value from being used as, say, a
//! velocity.
//!
//! Each constructor takes the *physical* quantities the number is built from,
//! so the formula that produced it is always visible at the call site.

use core::fmt;

use crate::consts;
use crate::quantity::{
    AbsoluteTemperature, Density, DynamicViscosity, KinematicViscosity, Length, Power, Pressure,
    SurfaceTension, Time, Velocity, Volume, VolumetricFlow,
};

/// Declares a non-dimensional number type with its symbol for display.
macro_rules! nondimensional {
    (
        $(#[$meta:meta])*
        $name:ident, $symbol:literal
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, PartialOrd, Default)]
        pub struct $name(f64);

        impl $name {
            /// Wraps an already-computed non-dimensional value.
            #[inline]
            pub const fn from_raw(value: f64) -> Self {
                Self(value)
            }

            /// The value of this number.
            #[inline]
            pub const fn value(self) -> f64 {
                self.0
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, concat!(stringify!($name), "({})"), self.0)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{} ({})", self.0, $symbol)
            }
        }
    };
}

nondimensional!(
    /// The Reynolds number `Re = V L / nu`, comparing inertial to viscous
    /// forces. `Re < 2300` is laminar flow; above roughly 4000 the flow is
    /// turbulent. This is the number that selects the friction correlation in
    /// `tpt-fluids-hydraulic`.
    Reynolds, "Re"
);

nondimensional!(
    /// The Froude number `Fr = V / sqrt(g L)`, comparing inertial to gravity
    /// forces. `Fr < 1` describes a subcritical free-surface flow; the
    /// ship-resistance scaling in `tpt-fluids-marine` is formulated entirely
    /// in terms of this number.
    Froude, "Fr"
);

nondimensional!(
    /// The Weber number `We = rho V^2 L / sigma`, comparing inertial to
    /// surface-tension forces. Large values mean surface tension is
    /// negligible; small values mean it dominates.
    Weber, "We"
);

nondimensional!(
    /// The Mach number `Ma = V / c`, comparing flow speed to the local speed
    /// of sound.
    Mach, "Ma"
);

nondimensional!(
    /// The cavitation number `sigma = (p - p_v) / (0.5 rho V^2)`, measuring
    /// how far the local pressure sits above the fluid's vapour pressure.
    /// Values approaching zero mark incipient cavitation.
    Cavitation, "sigma"
);

impl Reynolds {
    /// Builds a Reynolds number from velocity, characteristic length, and
    /// dynamic viscosity via `Re = rho V L / mu`.
    pub fn from_viscosity(
        velocity: Velocity,
        length: Length,
        mu: DynamicViscosity,
        rho: Density,
    ) -> Self {
        let denom = mu.value();
        if denom == 0.0 {
            return Self(f64::INFINITY);
        }
        Self(rho.value() * velocity.value() * length.value() / denom)
    }

    /// Builds a Reynolds number from velocity, characteristic length, and
    /// kinematic viscosity via `Re = V L / nu`.
    pub fn from_kinematic(velocity: Velocity, length: Length, nu: KinematicViscosity) -> Self {
        let denom = nu.value();
        if denom == 0.0 {
            return Self(f64::INFINITY);
        }
        Self(velocity.value() * length.value() / denom)
    }
}

impl Froude {
    /// Builds a Froude number from velocity and a characteristic length
    /// using standard gravity: `Fr = V / sqrt(g L)`.
    pub fn new(velocity: Velocity, length: Length) -> Self {
        Self::with_gravity(velocity, length, consts::STANDARD_GRAVITY)
    }

    /// Builds a Froude number with an explicitly supplied gravitational
    /// acceleration, as required for reduced-gravity manoeuvring studies.
    pub fn with_gravity(velocity: Velocity, length: Length, gravity: f64) -> Self {
        let denom = crate::math::sqrt(gravity * length.value());
        if denom == 0.0 {
            return Self(f64::INFINITY);
        }
        Self(velocity.value() / denom)
    }
}

impl Weber {
    /// Builds a Weber number: `We = rho V^2 L / sigma`.
    pub fn new(
        rho: Density,
        velocity: Velocity,
        length: Length,
        surface_tension: SurfaceTension,
    ) -> Self {
        let denom = surface_tension.value();
        if denom == 0.0 {
            return Self(f64::INFINITY);
        }
        Self(rho.value() * velocity.value() * velocity.value() * length.value() / denom)
    }
}

impl Mach {
    /// Builds a Mach number from flow speed and local speed of sound.
    pub fn new(velocity: Velocity, speed_of_sound: Velocity) -> Self {
        let denom = speed_of_sound.value();
        if denom == 0.0 {
            return Self(f64::INFINITY);
        }
        Self(velocity.value() / denom)
    }

    /// Builds a Mach number from a stagnation-to-static pressure ratio via
    /// the isentropic relation `Ma = sqrt((2 / gamma) (p0/p - 1))`, with
    /// `gamma = 1.4` for air.
    pub fn from_stagnation_pressure_ratio(stagnation: Pressure, static_p: Pressure) -> Self {
        let ratio = stagnation.value() / static_p.value();
        Self(crate::math::sqrt(2.0 / 1.4 * (ratio - 1.0).max(0.0)))
    }
}

impl Cavitation {
    /// Builds a cavitation number from the local absolute pressure, the
    /// fluid's vapour pressure, density, and velocity.
    pub fn new(
        pressure: Pressure,
        vapour_pressure: Pressure,
        rho: Density,
        velocity: Velocity,
    ) -> Self {
        let dynamic = 0.5 * rho.value() * velocity.value() * velocity.value();
        if dynamic == 0.0 {
            return Self(f64::INFINITY);
        }
        Self((pressure.value() - vapour_pressure.value()) / dynamic)
    }
}

/// The characteristic velocity implied by a power budget and a density,
/// `c = sqrt(P / rho)`. Returns zero for a non-positive input rather than a
/// NaN, so it is safe to feed from a clamped solver.
pub fn wave_velocity(power: Power, rho: Density) -> Velocity {
    let v2 = power.value() / rho.value();
    if v2 <= 0.0 {
        return Velocity::new(0.0);
    }
    Velocity::new(crate::math::sqrt(v2))
}

/// The time for a `volume` of fluid to pass at a given `flow` rate.
pub fn transit_time(volume: Volume, flow: VolumetricFlow) -> Time {
    if flow.value() == 0.0 {
        return Time::new(f64::INFINITY);
    }
    Time::new(volume.value() / flow.value())
}

/// An absolute temperature in kelvin, from degrees Celsius.
pub fn celsius_to_kelvin(celsius: f64) -> AbsoluteTemperature {
    AbsoluteTemperature::new(celsius + 273.15)
}

/// Degrees Celsius for an absolute temperature in kelvin.
pub fn kelvin_to_celsius(kelvin: AbsoluteTemperature) -> f64 {
    kelvin.value() - 273.15
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quantity::{Area, Force, MassFlow, VelocityPerArea};

    #[test]
    fn reynolds_from_kinematic_matches_reference() {
        // Water at 20 C: nu ~ 1.004e-6 m^2/s, V = 1 m/s, L = 1 m -> Re ~ 9.96e5.
        let re = Reynolds::from_kinematic(
            Velocity::new(1.0),
            Length::new(1.0),
            KinematicViscosity::new(1.004e-6),
        );
        assert!(
            (re.value() - 996_016.0).abs() / 996_016.0 < 1e-3,
            "{}",
            re.value()
        );
    }

    #[test]
    fn reynolds_from_dynamic_matches_kinematic_form() {
        let rho = Density::new(998.2);
        let nu = KinematicViscosity::new(1.004e-6);
        let mu = rho * nu;
        let a = Reynolds::from_kinematic(Velocity::new(2.5), Length::new(0.1), nu);
        let b = Reynolds::from_viscosity(Velocity::new(2.5), Length::new(0.1), mu, rho);
        assert!((a.value() - b.value()).abs() < 1e-6);
    }

    #[test]
    fn froude_of_one_is_the_critical_speed() {
        // Fr = 1 at V = sqrt(g L) = sqrt(9.80665) ~ 3.1312 m/s for L = 1 m.
        let fr = Froude::new(Velocity::new(3.1312), Length::new(1.0));
        assert!((fr.value() - 1.0).abs() < 1e-3, "{}", fr.value());
    }

    #[test]
    fn weber_and_cavitation_evaluate_correctly() {
        // We = rho V^2 L / sigma = 1000 * 1 * 1 / 0.072 = 13888.9
        let we = Weber::new(
            Density::new(1000.0),
            Velocity::new(1.0),
            Length::new(1.0),
            SurfaceTension::new(0.072),
        );
        assert!((we.value() - 13_888.9).abs() < 1.0, "{}", we.value());

        // sigma = (p - pv) / (0.5 rho V^2) = (200000 - 2000) / 500 = 396
        let sigma = Cavitation::new(
            Pressure::new(200_000.0),
            Pressure::new(2_000.0),
            Density::new(1000.0),
            Velocity::new(1.0),
        );
        assert!((sigma.value() - 396.0).abs() < 1e-9, "{}", sigma.value());
    }

    #[test]
    fn mach_number_is_velocity_over_sound_speed() {
        let m = Mach::new(Velocity::new(340.0), Velocity::new(1700.0));
        assert!((m.value() - 0.2).abs() < 1e-12);
        // The isentropic relation p0/p = 1 + gamma Ma^2 / 2 gives p0/p =
        // 1.028 for Ma = 0.2 in air.
        let s = Mach::from_stagnation_pressure_ratio(Pressure::new(1.028), Pressure::new(1.0));
        assert!((s.value() - 0.2).abs() < 1e-3, "{}", s.value());
    }

    #[test]
    fn zero_denominators_yield_infinity_not_nan() {
        assert!(Reynolds::from_kinematic(
            Velocity::new(1.0),
            Length::new(1.0),
            KinematicViscosity::new(0.0)
        )
        .value()
        .is_infinite());
        assert!(Cavitation::new(
            Pressure::new(1.0),
            Pressure::new(0.0),
            Density::new(1.0),
            Velocity::new(0.0)
        )
        .value()
        .is_infinite());
    }

    #[test]
    fn helpers_behave() {
        let c = celsius_to_kelvin(20.0);
        assert!((c.value() - 293.15).abs() < 1e-9);
        assert!((kelvin_to_celsius(c) - 20.0).abs() < 1e-9);

        let t = transit_time(Volume::new(10.0), VolumetricFlow::new(2.0));
        assert!((t.value() - 5.0).abs() < 1e-12);

        let v = wave_velocity(Power::new(2.0), Density::new(2.0));
        assert!((v.value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn conversions_reuse_quantity_algebra() {
        let q = VolumetricFlow::new(3.0);
        let v: VelocityPerArea = q / Area::new(4.0);
        assert!((v.value() - 0.75).abs() < 1e-12);
        let f: Force = Power::new(10.0) / Velocity::new(2.0);
        assert!((f.value() - 5.0).abs() < 1e-12);
        let g: Force = Velocity::new(2.0) * Power::new(10.0);
        assert!((g.value() - 20.0).abs() < 1e-12);
        let _: MassFlow = Density::new(1000.0) * q;
    }
}

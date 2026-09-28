//! Compile-time unit-safe physical quantities.
//!
//! Every quantity is a distinct Rust type carrying a `Dimension` — the SI
//! exponent vector `(L, M, T, Θ)`. Because the types are distinct, the
//! compiler rejects dimensionally meaningless operations:
//!
//! ```compile_fail
//! use tpt_fluids_core::quantity::{Head, Pressure};
//! let h = Head::new(12.0);
//! let p = Pressure::new(101_325.0);
//! let _ = h + p; // error: cannot add `Head` to `Pressure`
//! ```
//!
//! Arithmetic is only *implemented* where it is physically meaningful, so the
//! type system also encodes the algebra of the quantities: a `VolumetricFlow`
//! divided by a `Length` is a [`VolumetricFlowPerLength`], not an `f64`.

use core::fmt;
use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// The SI dimension of a quantity, as a vector of exponents over
/// `(metre, kilogram, second, kelvin)`.
///
/// A `Pressure` is `(-1, 1, -2, 0)`; a `Head` is `(1, 0, 0, 0)`. Exponents may
/// be negative and are not required to sum to zero.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Dimension {
    /// Power of the metre.
    pub length: i32,
    /// Power of the kilogram.
    pub mass: i32,
    /// Power of the second.
    pub time: i32,
    /// Power of the kelvin.
    pub temperature: i32,
}

impl Dimension {
    /// No dimension at all — the basis of every non-dimensional number.
    pub const DIMENSIONLESS: Self = Self::new(0, 0, 0, 0);
    /// Metre.
    pub const LENGTH: Self = Self::new(1, 0, 0, 0);
    /// Square metre.
    pub const AREA: Self = Self::new(2, 0, 0, 0);
    /// Cubic metre.
    pub const VOLUME: Self = Self::new(3, 0, 0, 0);
    /// Kilogram.
    pub const MASS: Self = Self::new(0, 1, 0, 0);
    /// Second.
    pub const TIME: Self = Self::new(0, 0, 1, 0);
    /// Kelvin (an absolute temperature).
    pub const ABSOLUTE_TEMPERATURE: Self = Self::new(0, 0, 0, 1);
    /// Kelvin (a temperature *difference*).
    pub const TEMPERATURE: Self = Self::new(0, 0, 1, 1);
    /// Metres per second.
    pub const VELOCITY: Self = Self::new(1, 0, -1, 0);
    /// Metres per second squared.
    pub const ACCELERATION: Self = Self::new(1, 0, -2, 0);
    /// Kilograms per cubic metre.
    pub const DENSITY: Self = Self::new(-3, 1, 0, 0);
    /// Radians per second. The radian is dimensionless, so this carries the
    /// reciprocal-second exponent in the time slot.
    pub const ANGULAR_RATE: Self = Self::new(0, 0, 0, -1);
    /// Newton-metres.
    pub const TORQUE: Self = Self::new(2, 1, -2, 0);
    /// Pascals.
    pub const PRESSURE: Self = Self::new(-1, 1, -2, 0);
    /// Newtons.
    pub const FORCE: Self = Self::new(1, 1, -2, 0);
    /// Watts.
    pub const POWER: Self = Self::new(2, 1, -3, 0);
    /// Joules.
    pub const ENERGY: Self = Self::new(2, 1, -2, 0);
    /// Cubic metres per second.
    pub const VOLUMETRIC_FLOW: Self = Self::new(3, 0, -1, 0);
    /// Kilograms per second.
    pub const MASS_FLOW: Self = Self::new(0, 1, -1, 0);
    /// Pascal-seconds.
    pub const DYNAMIC_VISCOSITY: Self = Self::new(-1, 1, -1, 0);
    /// Square metres per second.
    pub const KINEMATIC_VISCOSITY: Self = Self::new(2, 0, -1, 0);
    /// Newtons per metre (surface tension).
    pub const SURFACE_TENSION: Self = Self::new(0, 1, -2, 0);
    /// Radians per second.
    pub const ANGULAR_VELOCITY: Self = Self::new(0, 0, -1, 0);
    /// Kilograms per square metre.
    pub const MASS_PER_AREA: Self = Self::new(-2, 1, 0, 0);
    /// Pascals per metre.
    pub const PRESSURE_GRADIENT: Self = Self::new(-2, 1, -2, 0);
    /// Metres to the fourth power (second moment of area).
    pub const SECOND_MOMENT_OF_AREA: Self = Self::new(4, 0, 0, 0);

    /// Builds a dimension from its four exponents.
    pub const fn new(length: i32, mass: i32, time: i32, temperature: i32) -> Self {
        Self {
            length,
            mass,
            time,
            temperature,
        }
    }

    /// The dimension of `self * other`.
    pub const fn mul(self, other: Self) -> Self {
        Self::new(
            self.length + other.length,
            self.mass + other.mass,
            self.time + other.time,
            self.temperature + other.temperature,
        )
    }

    /// The dimension of `self / other`.
    pub const fn div(self, other: Self) -> Self {
        Self::new(
            self.length - other.length,
            self.mass - other.mass,
            self.time - other.time,
            self.temperature - other.temperature,
        )
    }

    /// The dimension of `self^n`, for integer `n`.
    pub const fn powi(self, n: i32) -> Self {
        Self::new(
            self.length * n,
            self.mass * n,
            self.time * n,
            self.temperature * n,
        )
    }

    /// Whether this dimension is the empty vector `(0, 0, 0, 0)`.
    pub const fn is_dimensionless(self) -> bool {
        self.length == 0 && self.mass == 0 && self.time == 0 && self.temperature == 0
    }
}

impl fmt::Display for Dimension {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "L^{} M^{} T^{} K^{}",
            self.length, self.mass, self.time, self.temperature
        )
    }
}

/// A physical quantity expressed in SI base units.
///
/// Implementors are thin `f64` newtypes. The dimension is carried as an
/// associated constant rather than a type parameter so that each physical
/// quantity is one concrete, `Copy` type with its own distinct identity.
pub trait Quantity: Copy + PartialOrd {
    /// The SI dimension of this quantity.
    const DIMENSION: Dimension;
    /// The quantity's magnitude in SI base units.
    fn value(self) -> f64;
}

/// Declares a dimensioned newtype with its constant and the unary/additive
/// operator set that every quantity shares.
macro_rules! dimensioned {
    (
        $(#[$meta:meta])*
        $name:ident, $unit:literal, $dim:expr
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, PartialOrd, Default)]
        pub struct $name(f64);

        impl $name {
            #[doc = concat!("Wraps a magnitude expressed in ", $unit, ".")]
            #[inline]
            pub const fn new(value: f64) -> Self {
                Self(value)
            }

            /// The additive identity.
            #[inline]
            pub const fn zero() -> Self {
                Self(0.0)
            }

            #[doc = concat!("This quantity's SI dimension, measured in ", $unit, ".")]
            pub const DIMENSION: Dimension = $dim;

            #[doc = concat!("The magnitude expressed in ", $unit, ".")]
            #[inline]
            pub const fn value(self) -> f64 {
                self.0
            }
        }

        impl Quantity for $name {
            const DIMENSION: Dimension = $dim;
            #[inline]
            fn value(self) -> f64 {
                self.0
            }
        }

        impl Add for $name {
            type Output = Self;
            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self(self.0 + rhs.0)
            }
        }

        impl Sub for $name {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self(self.0 - rhs.0)
            }
        }

        impl AddAssign for $name {
            #[inline]
            fn add_assign(&mut self, rhs: Self) {
                self.0 += rhs.0;
            }
        }

        impl SubAssign for $name {
            #[inline]
            fn sub_assign(&mut self, rhs: Self) {
                self.0 -= rhs.0;
            }
        }

        impl Neg for $name {
            type Output = Self;
            #[inline]
            fn neg(self) -> Self {
                Self(-self.0)
            }
        }

        impl Mul<f64> for $name {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: f64) -> Self {
                Self(self.0 * rhs)
            }
        }

        impl Mul<$name> for f64 {
            type Output = $name;
            #[inline]
            fn mul(self, rhs: $name) -> $name {
                $name(self * rhs.0)
            }
        }

        impl Div<f64> for $name {
            type Output = Self;
            #[inline]
            fn div(self, rhs: f64) -> Self {
                Self(self.0 / rhs)
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, concat!(stringify!($name), "({} ", $unit, ")"), self.0)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{} {}", self.0, $unit)
            }
        }
    };
}

/// Implements `A * B = C` and `A / B = C` (plus the commuted `B * A = C`) for
/// three dimensioned quantities whose product is dimensionally meaningful.
///
/// Self-products (`Length * Length`) would otherwise expand both the `A * B`
/// and `B * A` arms into the same impl, so those cases are written out by
/// hand just below rather than through this macro.
macro_rules! dimensioned_ops {
    ($a:ident, $b:ident, $c:ident) => {
        impl Mul<$b> for $a {
            type Output = $c;
            #[inline]
            fn mul(self, rhs: $b) -> $c {
                $c::new(self.value() * rhs.value())
            }
        }

        impl Mul<$a> for $b {
            type Output = $c;
            #[inline]
            fn mul(self, rhs: $a) -> $c {
                $c::new(self.value() * rhs.value())
            }
        }

        impl Div<$b> for $a {
            type Output = $c;
            #[inline]
            fn div(self, rhs: $b) -> $c {
                $c::new(self.value() / rhs.value())
            }
        }
    };
}

dimensioned!(
    /// A linear length, in metres.
    Length, "m", Dimension::LENGTH
);
dimensioned!(
    /// An area, in square metres.
    Area, "m^2", Dimension::AREA
);
dimensioned!(
    /// A volume, in cubic metres.
    Volume, "m^3", Dimension::VOLUME
);
dimensioned!(
    /// A mass, in kilograms.
    Mass, "kg", Dimension::MASS
);
dimensioned!(
    /// A duration, in seconds.
    Time, "s", Dimension::TIME
);
dimensioned!(
    /// An absolute temperature, in kelvin.
    AbsoluteTemperature, "K", Dimension::ABSOLUTE_TEMPERATURE
);
dimensioned!(
    /// A temperature difference, in kelvin.
    Temperature, "K", Dimension::TEMPERATURE
);
dimensioned!(
    /// A velocity, in metres per second.
    Velocity, "m/s", Dimension::VELOCITY
);
dimensioned!(
    /// An acceleration, in metres per second squared.
    Acceleration, "m/s^2", Dimension::ACCELERATION
);
dimensioned!(
    /// A mass density, in kilograms per cubic metre.
    Density, "kg/m^3", Dimension::DENSITY
);
dimensioned!(
    /// A pressure, in pascals.
    Pressure, "Pa", Dimension::PRESSURE
);
dimensioned!(
    /// An angular rate, in radians per second.
    ///
    /// Note this is radians per second, *not* revolutions per minute. The
    /// two are used interchangeably in propeller work, and confusing them is a
    /// factor of about 9.5 in any computed thrust, so the unit is stated
    /// rather than left to the reader.
    AngularRate, "rad/s", Dimension::ANGULAR_RATE
);
dimensioned!(
    /// A torque, in newton-metres.
    Torque, "N*m", Dimension::TORQUE
);
dimensioned!(
    /// A force, in newtons.
    Force, "N", Dimension::FORCE
);
dimensioned!(
    /// A power, in watts.
    Power, "W", Dimension::POWER
);
dimensioned!(
    /// An energy, in joules.
    Energy, "J", Dimension::ENERGY
);
dimensioned!(
    /// A volumetric flow rate, in cubic metres per second.
    VolumetricFlow, "m^3/s", Dimension::VOLUMETRIC_FLOW
);
dimensioned!(
    /// A mass flow rate, in kilograms per second.
    MassFlow, "kg/s", Dimension::MASS_FLOW
);
dimensioned!(
    /// Dynamic (Newtonian) viscosity, in pascal-seconds.
    DynamicViscosity, "Pa*s", Dimension::DYNAMIC_VISCOSITY
);
dimensioned!(
    /// Kinematic viscosity, in square metres per second.
    KinematicViscosity, "m^2/s", Dimension::KINEMATIC_VISCOSITY
);
dimensioned!(
    /// Surface tension, in newtons per metre.
    SurfaceTension, "N/m", Dimension::SURFACE_TENSION
);
dimensioned!(
    /// An angular velocity, in radians per second.
    AngularVelocity, "rad/s", Dimension::ANGULAR_VELOCITY
);
dimensioned!(
    /// An areal mass density, in kilograms per square metre.
    MassPerArea, "kg/m^2", Dimension::MASS_PER_AREA
);
dimensioned!(
    /// A pressure gradient, in pascals per metre.
    PressureGradient, "Pa/m", Dimension::PRESSURE_GRADIENT
);
dimensioned!(
    /// A second moment of area, in metres to the fourth power.
    SecondMomentOfArea, "m^4", Dimension::SECOND_MOMENT_OF_AREA
);
dimensioned!(
    /// A pressure expressed as an equivalent column of fluid, in metres.
    ///
    /// This is `Head` in the hydraulic sense: the height of fluid producing a
    /// given pressure. It shares its dimension with [`Length`] but is a
    /// distinct type, so a geometrical length can never be silently used
    /// where a hydraulic head is meant.
    Head, "m", Dimension::LENGTH
);
dimensioned!(
    /// A volumetric flow rate per unit length, in square metres per second.
    VolumetricFlowPerLength, "m^2/s", Dimension::VOLUMETRIC_FLOW.div(Dimension::LENGTH)
);
dimensioned!(
    /// A volumetric flow rate per unit area, in metres per second.
    VelocityPerArea, "m/s", Dimension::VOLUMETRIC_FLOW.div(Dimension::AREA)
);

dimensioned_ops!(Area, Length, Volume);
dimensioned_ops!(Length, Time, Velocity);
dimensioned_ops!(Velocity, Time, Length);
dimensioned_ops!(Area, Velocity, VolumetricFlow);
dimensioned_ops!(Length, Velocity, VolumetricFlow);
dimensioned_ops!(Mass, Volume, Density);
dimensioned_ops!(Density, VolumetricFlow, MassFlow);
dimensioned_ops!(Density, KinematicViscosity, DynamicViscosity);
dimensioned_ops!(VolumetricFlow, Length, VolumetricFlowPerLength);
dimensioned_ops!(VolumetricFlow, Area, VelocityPerArea);
dimensioned_ops!(Energy, Volume, Pressure);

/// `Power / Velocity = Force`, the work rate per unit of a body's speed.
///
/// Hand-written rather than generated: `Velocity * Time` is already `Length`,
/// and routing `Power * Velocity` through the macro would collide on the
/// commuted arm.
impl Div<Velocity> for Power {
    type Output = Force;
    #[inline]
    fn div(self, rhs: Velocity) -> Force {
        Force::new(self.value() / rhs.value())
    }
}

/// The commuted product, `Velocity * Power = Force`.
impl Mul<Power> for Velocity {
    type Output = Force;
    #[inline]
    fn mul(self, rhs: Power) -> Force {
        Force::new(self.value() * rhs.value())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimension_algebra_is_correct() {
        assert_eq!(Dimension::PRESSURE.mul(Dimension::AREA), Dimension::FORCE);
        assert_eq!(Dimension::FORCE.div(Dimension::AREA), Dimension::PRESSURE);
        assert_eq!(Dimension::LENGTH.powi(3), Dimension::VOLUME);
        assert!(Dimension::DIMENSIONLESS.is_dimensionless());
        assert!(!Dimension::LENGTH.is_dimensionless());
    }

    #[test]
    fn products_produce_the_right_types() {
        // Q = A * V -> a cubic metre per second.
        let a = Area::new(2.0);
        let v = Velocity::new(3.0);
        let q: VolumetricFlow = a * v;
        assert!((q.value() - 6.0).abs() < 1e-12);

        // m_dot = rho * Q.
        let rho = Density::new(1000.0);
        let mdot: MassFlow = rho * q;
        assert!((mdot.value() - 6000.0).abs() < 1e-9);

        // mu = rho * nu.
        let nu = KinematicViscosity::new(1.0e-6);
        let mu: DynamicViscosity = rho * nu;
        assert!((mu.value() - 1.0e-3).abs() < 1e-15);
    }

    #[test]
    fn head_and_length_are_distinct_types() {
        // These two share a dimension but must not be interchangeable.
        let head = Head::new(30.0);
        let length = Length::new(30.0);
        assert_eq!(head.value(), length.value());
        assert_eq!(Head::DIMENSION, Length::DIMENSION);
    }

    #[test]
    fn additive_operators_behave() {
        let mut p = Pressure::new(1.0);
        p += Pressure::new(2.0);
        p -= Pressure::new(0.5);
        assert!((p.value() - 2.5).abs() < 1e-12);
        assert!(((-p).value() + 2.5).abs() < 1e-12);
        assert_eq!(Pressure::default().value(), 0.0);
        assert_eq!(Pressure::zero().value(), 0.0);
    }

    #[test]
    fn scalar_scaling_behaves() {
        let v = Velocity::new(4.0) * 2.5;
        assert!((v.value() - 10.0).abs() < 1e-12);
        let w: Velocity = 2.0 * Velocity::new(4.0);
        assert!((w.value() - 8.0).abs() < 1e-12);
        let half = Pressure::new(10.0) / 4.0;
        assert!((half.value() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn derived_flow_quantities_behave() {
        let q = VolumetricFlow::new(3.0);
        let per_len: VolumetricFlowPerLength = q / Length::new(2.0);
        assert!((per_len.value() - 1.5).abs() < 1e-12);
        let vel: VelocityPerArea = q / Area::new(4.0);
        assert!((vel.value() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn display_reports_unit() {
        let s = Pressure::new(101_325.0).to_string();
        assert!(s.contains("Pa"), "{s}");
        let d = format!("{:?}", Length::new(2.0));
        assert!(d.contains("Length"), "{d}");
    }
}

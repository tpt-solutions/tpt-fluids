//! Viscosity models: Newtonian, power-law, Bingham plastic, and Sutherland's
//! law for gases.
//!
//! The first three are *shear-rate* relations, so they are asked for a
//! viscosity at a given `dV/dy`; Sutherland's is a pure temperature
//! correlation. [`ViscosityModel`] unifies them so a solver can carry one
//! field type.

use crate::math;
use crate::quantity::{
    AbsoluteTemperature, Density, DynamicViscosity, KinematicViscosity, Temperature,
};

/// A rheological model relating shear stress to shear rate and temperature.
///
/// The first three variants are shear-rate relations, so they are asked for a
/// viscosity at a given `dV/dy`; [`Sutherland`](Self::Sutherland) is a pure
/// temperature correlation. One type covers all of them so a solver can carry
/// a single viscosity field.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ViscosityModel {
    /// A constant-viscosity Newtonian fluid: `tau = mu dV/dy`.
    Newtonian {
        /// The constant dynamic viscosity, in pascal-seconds.
        mu: f64,
    },
    /// A power-law (Ostwald-de Waele) fluid: `tau = K (dV/dy)^n`, so the
    /// apparent viscosity is `mu = K (dV/dy)^(n - 1)`.
    PowerLaw {
        /// The consistency index `K`, in pascal-seconds to the power `n`.
        consistency: f64,
        /// The flow-behaviour index `n`. `n < 1` shear-thins, `n > 1`
        /// shear-thickens.
        exponent: f64,
    },
    /// A Bingham plastic: `tau = tau_y + mu_p dV/dy` once the yield stress
    /// is exceeded, and rigid below it.
    BinghamPlastic {
        /// The yield stress, in pascals.
        yield_stress: f64,
        /// The plastic viscosity, in pascal-seconds.
        plastic_viscosity: f64,
    },
    /// Sutherland's law for a gas:
    /// `mu = mu0 (T0/S) (T/T0)^(3/2) (T0 + S) / (T + S)`.
    Sutherland {
        /// The reference viscosity, in pascal-seconds.
        mu0: f64,
        /// The reference temperature, in kelvin.
        t0: f64,
        /// Sutherland's constant, in kelvin.
        s: f64,
    },
}

impl ViscosityModel {
    /// A constant-viscosity Newtonian model.
    pub const fn newtonian(mu: f64) -> Self {
        Self::Newtonian { mu }
    }

    /// A power-law model.
    pub const fn power_law(consistency: f64, exponent: f64) -> Self {
        Self::PowerLaw {
            consistency,
            exponent,
        }
    }

    /// A Bingham plastic model.
    pub const fn bingham_plastic(yield_stress: f64, plastic_viscosity: f64) -> Self {
        Self::BinghamPlastic {
            yield_stress,
            plastic_viscosity,
        }
    }

    /// Sutherland's law for dry air at standard conditions.
    pub const fn air() -> Self {
        Self::Sutherland {
            mu0: crate::consts::SUTHERLAND_REFERENCE_VISCOSITY_AIR,
            t0: 273.15,
            s: crate::consts::SUTHERLAND_CONSTANT_AIR,
        }
    }

    /// The dynamic viscosity this model reports at the given conditions.
    ///
    /// `shear_rate` is `dV/dy` in reciprocal seconds and is ignored by
    /// Sutherland's law. A non-positive shear rate is treated as the
    /// un-sheared limit, which is what makes the power-law and Bingham
    /// models well defined at rest.
    pub fn dynamic_viscosity(
        &self,
        temperature: AbsoluteTemperature,
        shear_rate: f64,
    ) -> DynamicViscosity {
        let rate = if shear_rate > 0.0 { shear_rate } else { 0.0 };
        match *self {
            Self::Newtonian { mu } => DynamicViscosity::new(mu),
            Self::PowerLaw {
                consistency,
                exponent,
            } => {
                // mu = K * rate^(n - 1)
                DynamicViscosity::new(consistency * math::powf(rate, exponent - 1.0))
            }
            Self::BinghamPlastic {
                yield_stress,
                plastic_viscosity,
            } => {
                // tau = tau_y + mu_p rate, so mu = tau / rate.
                if rate == 0.0 {
                    // Unyielded: the material behaves as a solid of infinite
                    // apparent viscosity rather than reporting a finite one.
                    return DynamicViscosity::new(f64::INFINITY);
                }
                DynamicViscosity::new(yield_stress / rate + plastic_viscosity)
            }
            Self::Sutherland { mu0, t0, s } => {
                if temperature.value() <= 0.0 {
                    return DynamicViscosity::new(f64::NAN);
                }
                let t = temperature.value();
                // mu / mu0 = (T/T0)^(3/2) (T0 + S) / (T + S)
                let ratio = math::powf(t / t0, 1.5) * (t0 + s) / (t + s);
                DynamicViscosity::new(mu0 * ratio)
            }
        }
    }

    /// The shear stress `tau` implied by a shear rate, in pascals.
    pub fn shear_stress(&self, temperature: AbsoluteTemperature, shear_rate: f64) -> f64 {
        let rate = if shear_rate > 0.0 { shear_rate } else { 0.0 };
        match *self {
            Self::Newtonian { mu } => mu * rate,
            Self::PowerLaw {
                consistency,
                exponent,
            } => consistency * math::powf(rate, exponent),
            Self::BinghamPlastic {
                yield_stress,
                plastic_viscosity,
            } => {
                if rate == 0.0 {
                    0.0
                } else {
                    yield_stress + plastic_viscosity * rate
                }
            }
            Self::Sutherland { .. } => {
                // Sutherland's law has no shear-rate dependence, so this is
                // the Newtonian stress of the temperature-dependent viscosity.
                self.dynamic_viscosity(temperature, rate).value() * rate
            }
        }
    }

    /// The kinematic viscosity `nu = mu / rho`, in square metres per second.
    pub fn kinematic_viscosity(
        &self,
        temperature: AbsoluteTemperature,
        shear_rate: f64,
        density: Density,
    ) -> KinematicViscosity {
        let mu = self.dynamic_viscosity(temperature, shear_rate).value();
        if density.value() == 0.0 {
            return KinematicViscosity::new(f64::NAN);
        }
        KinematicViscosity::new(mu / density.value())
    }

    /// The viscous-heating temperature rise across a `gap` at the given shear
    /// rate for a fluid of thermal conductivity `k`:
    /// `dT = mu (dV/dy)^2 gap / k`.
    pub fn viscous_heating(
        &self,
        temperature: AbsoluteTemperature,
        shear_rate: f64,
        gap: f64,
        conductivity: f64,
    ) -> Temperature {
        let mu = self.dynamic_viscosity(temperature, shear_rate).value();
        if conductivity == 0.0 {
            return Temperature::new(0.0);
        }
        Temperature::new(mu * shear_rate * shear_rate * gap / conductivity)
    }
}

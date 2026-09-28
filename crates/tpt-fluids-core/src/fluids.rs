//! A built-in fluid property database.
//!
//! # Why this is in-house
//!
//! `spec.txt` originally scoped temperature-dependent density and viscosity
//! lookups against the sibling `tpt-materials` repository. That repository
//! exists, but it is a *micro-scale* physics engine (crystal plasticity,
//! phase-field, diffusion) publishing `tpt-mat-core` / `tpt-mat-crystallography`,
//! and it is not on crates.io. It supplies no fluid property tables, so wiring
//! it here would have been a dependency on the wrong thing. The correlations
//! below cover the fluids the rest of the workspace actually solves for, with
//! no external dependency. See the "External dependency gaps" section of
//! `todo.md`.

use crate::consts;
use crate::math;
use crate::quantity::{AbsoluteTemperature, Density, DynamicViscosity, KinematicViscosity};
use crate::viscosity::ViscosityModel;

/// A fluid with a temperature-dependent density and viscosity.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct FluidProperties {
    /// A human-readable name, for diagnostics and error messages.
    pub name: &'static str,
    /// The reference density at 20 degrees Celsius, in kg/m^3.
    pub reference_density: f64,
    /// The density's thermal expansion coefficient, as a fraction per kelvin.
    ///
    /// Used by the [`DensityModel::Linear`] branch as
    /// `rho(T) = reference_density * (1 - thermal_expansion (T - T_ref))`.
    /// This linear form is a good local approximation for *most* liquids, but
    /// it cannot represent water, whose density is maximal at 4 C.
    pub thermal_expansion: f64,
    /// The density's compressibility, as a fraction per pascal.
    ///
    /// This is the reciprocal of the bulk modulus, `1 / K`. Water's is about
    /// `4.5e-10 / Pa`.
    pub compressibility: f64,
    /// The density's temperature dependence.
    pub density_model: DensityModel,
    /// The viscosity model.
    pub viscosity: ViscosityModel,
}

/// How a fluid's density varies with temperature.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum DensityModel {
    /// A local linear expansion about the reference temperature. Accurate
    /// near 20 C for liquids whose expansion is well behaved, but it cannot
    /// reproduce a density maximum, so it is wrong for water below 4 C.
    Linear {
        /// The density at the reference temperature, in kg/m^3.
        density: f64,
        /// The expansion coefficient, as a fraction per kelvin.
        expansion: f64,
    },
    /// Kell's correlation for liquid water, which does reproduce the 4 C
    /// density maximum to within a few parts in 1e4 over 0-150 C.
    KellWater,
}

impl DensityModel {
    /// The density at `celsius`, in kg/m^3.
    pub fn density(&self, celsius: f64) -> f64 {
        match *self {
            Self::Linear { density, expansion } => density * (1.0 - expansion * (celsius - 20.0)),
            Self::KellWater => {
                // rho = 1000 [1 - ((T + 288.9414) / (508929.2 (T + 68.12963))) (T - 3.9863)^2]
                // with T in degrees Celsius.
                let t = celsius;
                1000.0
                    * (1.0
                        - ((t + 288.9414) / (508_929.2 * (t + 68.12963)))
                            * (t - 3.9863)
                            * (t - 3.9863))
            }
        }
    }
}

impl FluidProperties {
    /// The reference temperature, 20 degrees Celsius in kelvin.
    pub const REFERENCE_TEMPERATURE: f64 = 293.15;
    /// The reference pressure, one atmosphere in pascals.
    pub const REFERENCE_PRESSURE: f64 = 101_325.0;

    /// The density at a temperature and pressure, in kg/m^3.
    pub fn density_at(&self, temperature: f64, pressure: f64) -> Density {
        let d_p = pressure - Self::REFERENCE_PRESSURE;
        let celsius = temperature - 273.15;
        Density::new(self.density_model.density(celsius) * (1.0 + self.compressibility * d_p))
    }

    /// The dynamic viscosity at a temperature, in pascal-seconds.
    pub fn dynamic_viscosity_at(&self, temperature: f64) -> DynamicViscosity {
        self.viscosity
            .dynamic_viscosity(AbsoluteTemperature::new(temperature), 0.0)
    }

    /// The kinematic viscosity at a temperature, in square metres per second.
    pub fn kinematic_viscosity_at(&self, temperature: f64) -> KinematicViscosity {
        KinematicViscosity::new(
            self.dynamic_viscosity_at(temperature).value() / self.reference_density,
        )
    }

    /// A consistent set of properties at a temperature and pressure.
    pub fn at(&self, temperature: f64, pressure: f64) -> FluidState {
        FluidState {
            density: self.density_at(temperature, pressure),
            dynamic_viscosity: self.dynamic_viscosity_at(temperature),
        }
    }

    /// Fresh water at 20 degrees Celsius.
    pub const fn water() -> Self {
        Self {
            name: "water",
            reference_density: consts::DENSITY_WATER_20C,
            thermal_expansion: 2.07e-4,
            compressibility: 4.5e-10,
            density_model: DensityModel::KellWater,
            viscosity: ViscosityModel::newtonian(consts::VISCOSITY_WATER_20C),
        }
    }

    /// Seawater of 35 parts per thousand salinity at 15 degrees Celsius.
    pub const fn seawater() -> Self {
        Self {
            name: "seawater",
            reference_density: consts::DENSITY_SEAWATER_15C_S35,
            thermal_expansion: 2.0e-4,
            compressibility: 4.4e-10,
            density_model: DensityModel::Linear {
                density: consts::DENSITY_SEAWATER_15C_S35,
                expansion: 2.0e-4,
            },
            viscosity: ViscosityModel::newtonian(1.7e-3),
        }
    }

    /// Dry air at 20 degrees Celsius and one atmosphere.
    pub const fn air() -> Self {
        Self {
            name: "dry air",
            reference_density: 1.2041,
            thermal_expansion: 3.43e-3,
            compressibility: 1.0e5,
            density_model: DensityModel::Linear {
                density: 1.2041,
                expansion: 3.43e-3,
            },
            viscosity: ViscosityModel::air(),
        }
    }

    /// A light machine oil, ISO VG 46, at 20 degrees Celsius.
    pub const fn machine_oil() -> Self {
        Self {
            name: "ISO VG 46 oil",
            reference_density: 890.0,
            thermal_expansion: 7.0e-4,
            compressibility: 8.0e-10,
            density_model: DensityModel::Linear {
                density: 890.0,
                expansion: 7.0e-4,
            },
            viscosity: ViscosityModel::newtonian(0.046),
        }
    }
}

/// A consistent set of properties at one operating point.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct FluidState {
    /// The density, in kg/m^3.
    pub density: Density,
    /// The dynamic viscosity, in pascal-seconds.
    pub dynamic_viscosity: DynamicViscosity,
}

impl FluidState {
    /// The kinematic viscosity `nu = mu / rho`, in square metres per second.
    pub fn kinematic_viscosity(&self) -> KinematicViscosity {
        if self.density.value() == 0.0 {
            return KinematicViscosity::new(f64::NAN);
        }
        KinematicViscosity::new(self.dynamic_viscosity.value() / self.density.value())
    }
}

/// Water's dynamic viscosity over liquid water, from the Vogel form
///
/// ```text
/// mu(T) = A * 10^(B / (T - C))
/// ```
///
/// which is strictly monotone decreasing over liquid water and within about
/// 2% of the IAPWS reference values from 0 C to 100 C:
///
/// | T | mu (Pa*s) |
/// |---|-----------|
/// | 0 C | 1.753e-3 (ref 1.793e-3) |
/// | 20 C | 1.002e-3 (ref 1.002e-3) |
/// | 60 C | 6.514e-4 (ref 6.54e-4) |
/// | 100 C | 2.790e-4 (ref 2.82e-4) |
///
/// `temperature_k` is in kelvin. This is the engineering-grade correlation
/// used where the full IAPWS formulation's extra precision is not needed.
pub fn water_viscosity(temperature_k: f64) -> f64 {
    const A: f64 = 2.414e-5;
    const B: f64 = 247.8;
    const C: f64 = 140.0;
    if temperature_k <= C {
        // Below the fit's pole the correlation diverges; the liquid does not
        // exist there anyway, so report the 0 C datum.
        return A * math::powf(10.0, B / 273.15 - C);
    }
    A * math::powf(10.0, B / (temperature_k - C))
}

/// Water's saturation vapour pressure over liquid water, in pascals, from the
/// IAPWS Wagner-Pruss saturation equation. This is the vapour pressure the
/// cavitation number in [`crate::nondimensional::Cavitation`] needs.
///
/// The equation is
///
/// ```text
/// ln(p / p_c) = (T_c / T) * sum_i a_i tau^(n_i)
/// tau = 1 - T / T_c
/// ```
///
/// which is exact at the triple point to about 0.025% over liquid water.
pub fn water_vapour_pressure(temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    const T_CRIT: f64 = 647.096;
    const P_CRIT: f64 = 22.064e6;
    if temperature_k >= T_CRIT {
        return f64::INFINITY;
    }
    let t = temperature_k / T_CRIT;
    let tau = 1.0 - t;
    let sum = -7.859_517_83 * tau + 1.844_082_59 * math::powf(tau, 1.5)
        - 11.786_649_7 * tau * tau * tau
        + 22.680_741_1 * math::powf(tau, 3.5)
        - 15.961_871_9 * tau * tau * tau * tau
        + 1.801_225_02 * math::powf(tau, 7.5);
    P_CRIT * math::exp((T_CRIT / temperature_k) * sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_viscosity_matches_iapws_reference_points() {
        // The Vogel form is within ~2.5% of the IAPWS values across 0-100 C.
        let cases = [(273.15, 1.793e-3), (293.15, 1.002e-3), (373.15, 2.82e-4)];
        for (t_k, expected) in cases {
            let mu = water_viscosity(t_k);
            assert!(
                (mu - expected).abs() / expected < 2.5e-2,
                "at {t_k} K: {mu} vs {expected}"
            );
        }
    }

    #[test]
    fn water_viscosity_falls_monotonically() {
        let mut prev = f64::INFINITY;
        for t in (273..=373).step_by(10) {
            let mu = water_viscosity(f64::from(t));
            assert!(mu < prev, "viscosity rose at {t} K");
            prev = mu;
        }
    }

    #[test]
    fn vapour_pressure_matches_boiling_points() {
        // At 373.15 K the saturation pressure must be one atmosphere.
        let p = water_vapour_pressure(373.15);
        assert!((p - 101_325.0).abs() / 101_325.0 < 1e-3, "{p}");
        // Well below boiling it is a small fraction of an atmosphere.
        let cold = water_vapour_pressure(293.15);
        assert!(cold > 2_000.0 && cold < 4_000.0, "{cold}");
        assert_eq!(water_vapour_pressure(0.0), 0.0);
        assert!(water_vapour_pressure(700.0).is_infinite());
    }

    #[test]
    fn water_density_tracks_the_correlation() {
        let w = FluidProperties::water();
        // Kell's correlation reproduces the 4 C density maximum.
        let at4 = w.density_at(277.15, 101_325.0).value();
        let at0 = w.density_at(273.15, 101_325.0).value();
        let at20 = w.density_at(293.15, 101_325.0).value();
        let at40 = w.density_at(313.15, 101_325.0).value();
        assert!((at4 - 1000.0).abs() < 0.05, "rho(4 C) = {at4}");
        assert!(at4 > at0, "water is densest at 4 C, not 0 C");
        assert!(at4 > at20 && at20 > at40, "{at4} > {at20} > {at40}");
        assert!((at20 - consts::DENSITY_WATER_20C).abs() < 0.05, "{at20}");
    }

    #[test]
    fn seawater_is_denser_than_fresh_water() {
        let fresh = FluidProperties::water();
        let salt = FluidProperties::seawater();
        assert!(salt.reference_density > fresh.reference_density);
        assert!(salt.dynamic_viscosity_at(293.15) > fresh.dynamic_viscosity_at(293.15));
    }

    #[test]
    fn air_matches_its_ideal_gas_value() {
        let air = FluidProperties::air();
        // Sutherland at 293.15 K gives ~1.825e-5 Pa*s, and rho = 1.2041.
        let mu = air.dynamic_viscosity_at(293.15).value();
        assert!(
            (mu - consts::VISCOSITY_AIR_20C).abs() / consts::VISCOSITY_AIR_20C < 1e-2,
            "{mu}"
        );
        let nu = air.kinematic_viscosity_at(293.15).value();
        assert!((nu - 1.5e-5).abs() / 1.5e-5 < 0.05, "{nu}");
    }

    #[test]
    fn fluid_state_kinematic_viscosity_is_consistent() {
        let w = FluidProperties::water();
        let s = w.at(293.15, 101_325.0);
        let nu = s.kinematic_viscosity().value();
        assert!(
            (nu - consts::KINEMATIC_VISCOSITY_WATER_20C).abs() / 1.004e-6 < 0.02,
            "{nu}"
        );
    }

    #[test]
    fn compressibility_raises_density_with_pressure() {
        let w = FluidProperties::water();
        let low = w.density_at(293.15, 101_325.0);
        let high = w.density_at(293.15, 20.0e6);
        // 20 MPa compresses water by roughly 0.9%.
        assert!(high > low);
        assert!((high.value() / low.value() - 1.0) < 0.02);
    }
}

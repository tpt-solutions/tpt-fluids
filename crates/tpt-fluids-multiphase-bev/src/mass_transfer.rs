//! CO2 diffusivity, liquid-side mass-transfer coefficient, and diffusion-limited
//! bubble growth.
//!
//! Diffusivity in the formulation follows Stokes-Einstein scaling from water,
//! `D = D_w (T / T_w) (mu_w / mu)`, which is the right trend for a viscous
//! syrup even though the absolute value needs the measured viscosity ratio.
//!
//! The liquid-side coefficient is the Higbie penetration result for a
//! stagnant surface, with the Frossling correlation for a moving bubble:
//!
//! ```text
//! Sh = 2 + 0.6 Re^(1/2) Sc^(1/3)      (bubble, Re and Sc from its slip)
//! k_L = Sh D / d                      (bubble)
//! k_L = 2 sqrt(D / (pi t))            (penetration, contact time t)
//! ```
//!
//! A growing bubble in a supersaturated liquid is limited by diffusion of CO2
//! into it; the Epstein-Plesset result gives its radius over time:
//!
//! ```text
//! dR/dt = D dc / (rho_g R) * (1 + R / sqrt(pi D t)) ...
//! ```
//!
//! Only the first-order (no-surface-tension) growth rate is used here, which is
//! the leading term for the millimetre-scale bubbles that matter for foam.

use tpt_fluids_core::math;

use crate::consts::{BOLTZMANN, MOLAR_MASS_CO2};
use crate::error::{positive, Result};
use crate::salting_out::Formulation;

/// Reference water diffusivity of CO2 at 25 C, in m^2/s.
const D_CO2_WATER_25C: f64 = 1.92e-9;
/// Reference water viscosity at 25 C, in Pa s.
const MU_WATER_25C: f64 = 8.9e-4;
/// Reference temperature, in kelvin.
const T_25C: f64 = 298.15;

/// The CO2 diffusivity in the formulation at `t_k`, in m^2/s.
pub fn diffusivity(formulation: &Formulation, t_k: f64) -> Result<f64> {
    positive(t_k, "temperature")?;
    let mu = tpt_fluids_core::fluids::water_viscosity(t_k) * formulation.viscosity_ratio;
    positive(mu, "solution viscosity")?;
    Ok(D_CO2_WATER_25C * (t_k / T_25C) * (MU_WATER_25C / mu))
}

/// The Schmidt number `nu / D` of CO2 in the formulation.
pub fn schmidt(formulation: &Formulation, t_k: f64) -> Result<f64> {
    let d = diffusivity(formulation, t_k)?;
    let mu = tpt_fluids_core::fluids::water_viscosity(t_k) * formulation.viscosity_ratio;
    Ok(mu / (1000.0 * d))
}

/// The Higbie penetration coefficient for a stagnant surface contacted for
/// `t_contact` seconds, in m/s.
pub fn penetration_coefficient(d: f64, t_contact: f64) -> Result<f64> {
    positive(d, "diffusivity")?;
    positive(t_contact, "contact time")?;
    Ok(2.0 * math::sqrt(d / (core::f64::consts::PI * t_contact)))
}

/// The liquid-side coefficient for a bubble of diameter `d_bubble` rising at
/// slip velocity `slip`, from the Frossling correlation, in m/s.
pub fn bubble_mass_transfer(
    formulation: &Formulation,
    t_k: f64,
    d_bubble: f64,
    slip: f64,
) -> Result<f64> {
    positive(d_bubble, "bubble diameter")?;
    positive(slip, "slip velocity")?;
    let d = diffusivity(formulation, t_k)?;
    let nu = tpt_fluids_core::fluids::water_viscosity(t_k) * formulation.viscosity_ratio / 1000.0;
    let re = slip * d_bubble / nu;
    let sc = nu / d;
    let sh = 2.0 + 0.6 * math::sqrt(re) * math::powf(sc, 1.0 / 3.0);
    Ok(sh * d / d_bubble)
}

/// The rate of radius growth of a bubble of radius `r` (metres) in a liquid
/// with dissolved-gas excess `dc_mol_per_m3` (the concentration above the
/// bubble's equilibrium), in m/s. Diffusion-limited, first-order in the
/// excess; `rho_gas` is the gas density in kg/m^3.
pub fn growth_rate(d: f64, dc_mol_per_m3: f64, r: f64, rho_gas: f64) -> Result<f64> {
    positive(d, "diffusivity")?;
    positive(r, "bubble radius")?;
    positive(rho_gas, "gas density")?;
    if dc_mol_per_m3 <= 0.0 {
        return Ok(0.0);
    }
    // Mass flux into the bubble, mol/(m^2 s), divided by the molar gas
    // density at the bubble surface, gives the surface regression rate.
    let molar_density = rho_gas / MOLAR_MASS_CO2;
    Ok(d * dc_mol_per_m3 / (molar_density * r))
}

/// The thermal energy `k_B T`, in joules, for convenience in rate expressions.
pub fn thermal_energy(t_k: f64) -> f64 {
    BOLTZMANN * t_k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_diffusivity_is_recovered_at_reference() {
        let d = diffusivity(&Formulation::water(), T_25C).unwrap();
        assert!((d - D_CO2_WATER_25C).abs() / D_CO2_WATER_25C < 0.02, "{d}");
    }

    #[test]
    fn syrup_slows_diffusion() {
        let syrup = Formulation::water().with_viscosity_ratio(3.0).unwrap();
        let d_w = diffusivity(&Formulation::water(), 277.15).unwrap();
        let d_s = diffusivity(&syrup, 277.15).unwrap();
        assert!(d_s < d_w / 2.5);
    }

    #[test]
    fn penetration_coefficient_falls_with_contact_time() {
        let early = penetration_coefficient(D_CO2_WATER_25C, 0.01).unwrap();
        let late = penetration_coefficient(D_CO2_WATER_25C, 1.0).unwrap();
        assert!(late < early);
    }

    #[test]
    fn bubble_transfer_is_at_least_the_stagnant_limit() {
        // Sh >= 2 means k_L >= 2 D / d.
        let k = bubble_mass_transfer(&Formulation::water(), T_25C, 1.0e-3, 0.2).unwrap();
        assert!(k >= 2.0 * D_CO2_WATER_25C / 1.0e-3);
    }

    #[test]
    fn growth_is_zero_without_excess_and_faster_for_small_bubbles() {
        assert_eq!(growth_rate(D_CO2_WATER_25C, 0.0, 1e-5, 1.8).unwrap(), 0.0);
        let small = growth_rate(D_CO2_WATER_25C, 50.0, 1e-5, 1.8).unwrap();
        let big = growth_rate(D_CO2_WATER_25C, 50.0, 1e-4, 1.8).unwrap();
        assert!(small > big);
    }
}

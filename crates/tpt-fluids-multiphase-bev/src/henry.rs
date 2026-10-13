//! Henry's law for CO2 in pure water.
//!
//! The solubility of CO2 is taken from the Weiss (1974) correlation, which is
//! fitted to measured solubilities between 0 and 40 C and is the standard
//! reference for carbonation work:
//!
//! ```text
//! ln K0 = -60.2409 + 93.4517 (100/T) + 23.3585 ln(T/100)
//! ```
//!
//! with `K0` in mol/(kg atm) of *fugacity*. The fugacity coefficient comes
//! from the Weiss-Price virial form, so the model stays accurate at bottling
//! pressures of several atmospheres rather than assuming an ideal gas.
//!
//! | T (C) | K0 (mol/(kg atm)) | g CO2 per kg water at 1 atm |
//! |-------|-------------------|-----------------------------|
//! | 0     | 0.0777            | 3.42                        |
//! | 25    | 0.0345            | 1.52                        |
//! | 40    | 0.0207            | 0.91                        |
//!
//! Salts and sugars reduce these values; see [`crate::salting_out`].

use tpt_fluids_core::fluids::water_vapour_pressure;
use tpt_fluids_core::math;

use crate::consts::{ATM, MOLAR_GAS_CONSTANT, MOLAR_MASS_CO2, MOLAR_VOLUME_STP};
use crate::error::{positive, BevError, Result};

/// The temperature range, in kelvin, the correlations here are fitted over.
pub const VALID_RANGE_K: (f64, f64) = (273.15, 323.15);

fn check_temperature(t_k: f64) -> Result<f64> {
    positive(t_k, "temperature")?;
    if t_k < VALID_RANGE_K.0 - 5.0 || t_k > VALID_RANGE_K.1 + 5.0 {
        return Err(BevError::OutsideValidRange("0-50 C"));
    }
    Ok(t_k)
}

/// The mass density of pure liquid water, in kg/m^3, at a temperature in kelvin.
pub fn water_density(t_k: f64) -> f64 {
    tpt_fluids_core::fluids::FluidProperties::water()
        .density_at(t_k, ATM)
        .value()
}

/// The Henry solubility `K0` of CO2 in pure water, in mol per kg of water per
/// pascal of CO2 fugacity, from the Weiss (1974) correlation.
pub fn henry_pure_mol_per_kg_pa(t_k: f64) -> Result<f64> {
    let t = check_temperature(t_k)?;
    let ln_k0 = -60.2409 + 93.4517 * (100.0 / t) + 23.3585 * math::ln(t / 100.0);
    Ok(math::exp(ln_k0) / ATM)
}

/// The Henry solubility of CO2 in pure water, in mol/(m^3 Pa), using the
/// density of water to convert from molality.
pub fn henry_pure_mol_per_m3_pa(t_k: f64) -> Result<f64> {
    Ok(henry_pure_mol_per_kg_pa(t_k)? * water_density(t_k))
}

/// The fugacity coefficient of pure CO2 at `p_pa` and `t_k`, from the
/// Weiss-Price virial form `ln phi = (B + 2 delta) p / (R T)`.
///
/// `B(T)` is the second virial coefficient of CO2 and `delta` the cross term
/// with water vapour. The form is good to about 10 atm, which covers every
/// bottled soft drink.
pub fn fugacity_coefficient(t_k: f64, p_pa: f64) -> Result<f64> {
    let t = check_temperature(t_k)?;
    positive(p_pa, "pressure")?;
    // cm^3/mol.
    let b = -1636.75 + 12.0408 * t - 0.032_795_7 * t * t + 3.165_28e-5 * t * t * t;
    let delta = 57.7 - 0.118 * t;
    let b_eff = (b + 2.0 * delta) * 1.0e-6;
    Ok(math::exp(b_eff * p_pa / (MOLAR_GAS_CONSTANT * t)))
}

/// The CO2 partial pressure, in pascals, of a headspace at total pressure
/// `p_total_pa` that is saturated with water vapour. Dissolved gas is in
/// equilibrium with this, not with the total.
pub fn co2_partial_pressure(t_k: f64, p_total_pa: f64) -> Result<f64> {
    let t = check_temperature(t_k)?;
    positive(p_total_pa, "total pressure")?;
    let p = p_total_pa - water_vapour_pressure(t);
    positive(p, "CO2 partial pressure")
}

/// The equilibrium dissolved CO2 in pure water at `p_co2_pa`, in mol per kg of
/// water, including the fugacity correction.
pub fn saturation_molality_pure(t_k: f64, p_co2_pa: f64) -> Result<f64> {
    let k0 = henry_pure_mol_per_kg_pa(t_k)?;
    let phi = fugacity_coefficient(t_k, p_co2_pa)?;
    Ok(k0 * phi * p_co2_pa)
}

/// Converts a dissolved concentration in mol/m^3 to grams per litre.
pub fn mol_per_m3_to_g_per_l(c: f64) -> f64 {
    c * MOLAR_MASS_CO2
}

/// Converts a dissolved concentration in mol/m^3 to "volumes of CO2": litres
/// of gas at 0 C and 1 atm per litre of liquid, the unit bottlers specify.
pub fn mol_per_m3_to_volumes(c: f64) -> f64 {
    c * MOLAR_VOLUME_STP
}

/// Converts volumes of CO2 to a dissolved concentration in mol/m^3.
pub fn volumes_to_mol_per_m3(volumes: f64) -> f64 {
    volumes / MOLAR_VOLUME_STP
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solubility_matches_measured_values() {
        // Mol/kg at 1 atm of CO2, from the standard tabulation.
        let cases = [(0.0, 0.0777), (25.0, 0.0345), (40.0, 0.0239)];
        for (c, expect) in cases {
            let k = henry_pure_mol_per_kg_pa(c + 273.15).unwrap() * ATM;
            assert!((k - expect).abs() / expect < 0.03, "{c} C: {k} vs {expect}");
        }
    }

    #[test]
    fn solubility_falls_with_temperature() {
        let mut last = f64::INFINITY;
        for t in (273..=323).step_by(2) {
            let k = henry_pure_mol_per_kg_pa(f64::from(t)).unwrap();
            assert!(k < last);
            last = k;
        }
    }

    #[test]
    fn one_atmosphere_at_25c_is_about_1_5_grams_per_kilogram() {
        let m = saturation_molality_pure(298.15, ATM).unwrap();
        let g = m * MOLAR_MASS_CO2 * 1000.0;
        assert!((g - 1.5).abs() < 0.1, "{g}");
    }

    #[test]
    fn fugacity_is_ideal_at_low_pressure_and_below_one_when_compressed() {
        let low = fugacity_coefficient(298.15, 100.0).unwrap();
        assert!((low - 1.0).abs() < 1e-5);
        let phi = fugacity_coefficient(298.15, 5.0 * ATM).unwrap();
        assert!(phi < 1.0 && phi > 0.97, "{phi}");
    }

    #[test]
    fn volume_conversions_round_trip() {
        let c = volumes_to_mol_per_m3(3.5);
        assert!((mol_per_m3_to_volumes(c) - 3.5).abs() < 1e-12);
        // One volume is about 1.96 g/L.
        assert!((mol_per_m3_to_g_per_l(volumes_to_mol_per_m3(1.0)) - 1.964).abs() < 0.01);
    }

    #[test]
    fn rejects_unphysical_inputs() {
        assert!(henry_pure_mol_per_kg_pa(-1.0).is_err());
        assert!(henry_pure_mol_per_kg_pa(500.0).is_err());
        assert!(co2_partial_pressure(298.15, 100.0).is_err());
    }
}

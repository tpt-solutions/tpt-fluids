//! Salting-out of CO2 by dissolved ions and polyols (Sechenov relation).
//!
//! Dissolved solutes lower gas solubility. The Sechenov (Setschenow) relation
//! gives the ratio of the Henry solubility in the solution to that in pure
//! water:
//!
//! ```text
//! log10(S0 / S) = sum_i (h_i + h_G(T)) c_i          (ions, Weisenberger-Schumpe)
//! log10(S0 / S) = sum_j Ks_j c_j                    (neutral solutes)
//! h_G(T)        = h_G0 + h_T (T - 298.15)
//! ```
//!
//! with `c` in mol/L and the coefficients in L/mol. The ion coefficients and
//! the CO2 gas parameters (`h_G0 = -0.0172`, `h_T = -3.38e-4`) are from the
//! Weisenberger-Schumpe (1996) tabulation as commonly reproduced; check them
//! against the source before relying on a result to better than a few percent.
//!
//! # Calibrate the polyol coefficients
//!
//! There is no well-established CO2 Sechenov constant for xylitol, nor for a
//! microencapsulating polymer. [`Solute::xylitol_estimate`] therefore returns
//! a *placeholder* tagged [`Provenance::Estimate`], of the same order as other
//! sugars and polyols, so that formulations can be compared relatively. Replace
//! it with a measured value before using absolute numbers.

use alloc::string::String;
use alloc::vec::Vec;

use tpt_fluids_core::math;

use crate::consts::T_REF;
use crate::error::{non_negative, positive, Result};

/// The Weisenberger-Schumpe CO2 gas parameter at 298.15 K, in L/mol.
const H_G0_CO2: f64 = -0.0172;
/// The temperature slope of the CO2 gas parameter, in L/(mol K).
const H_T_CO2: f64 = -3.38e-4;

/// Where a coefficient came from, so a screening result can say how far to
/// trust it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Provenance {
    /// A published, tabulated value.
    Literature,
    /// A fitted value from the user's own measurements.
    Measured,
    /// A placeholder of the right order of magnitude; not a measurement.
    Estimate,
}

/// How a solute salts CO2 out.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SaltingModel {
    /// An ion with a Weisenberger-Schumpe parameter `h_i`, in L/mol. The gas
    /// term `h_G(T)` is added.
    Ion {
        /// The ion-specific parameter, in L/mol.
        h_i: f64,
    },
    /// A neutral solute with a Sechenov constant `Ks`, in L/mol.
    Neutral {
        /// The Sechenov constant, in L/mol.
        ks: f64,
    },
}

/// A dissolved species that affects CO2 solubility and surface tension.
#[derive(Clone, PartialEq, Debug)]
pub struct Solute {
    /// A human-readable name.
    pub name: String,
    /// The salting-out model.
    pub model: SaltingModel,
    /// The surface-tension increment `d(sigma)/dc`, in N/m per mol/L. Ions
    /// raise surface tension; many polyols raise it slightly.
    pub dsigma_dc: f64,
    /// Where `model` and `dsigma_dc` came from.
    pub provenance: Provenance,
}

impl Solute {
    fn ion(name: &str, h_i: f64, dsigma_dc: f64) -> Self {
        Self {
            name: name.into(),
            model: SaltingModel::Ion { h_i },
            dsigma_dc,
            provenance: Provenance::Literature,
        }
    }

    /// The sodium ion.
    pub fn sodium() -> Self {
        Self::ion("Na+", 0.1143, 1.0e-3)
    }

    /// The potassium ion.
    pub fn potassium() -> Self {
        Self::ion("K+", 0.0922, 0.8e-3)
    }

    /// The calcium ion.
    pub fn calcium() -> Self {
        Self::ion("Ca2+", 0.1762, 2.0e-3)
    }

    /// The magnesium ion.
    pub fn magnesium() -> Self {
        Self::ion("Mg2+", 0.1694, 2.0e-3)
    }

    /// The chloride ion.
    pub fn chloride() -> Self {
        Self::ion("Cl-", 0.0318, 0.7e-3)
    }

    /// The bicarbonate ion.
    pub fn bicarbonate() -> Self {
        Self::ion("HCO3-", 0.0967, 0.5e-3)
    }

    /// The sulfate ion.
    pub fn sulfate() -> Self {
        Self::ion("SO4 2-", 0.1117, 1.5e-3)
    }

    /// A neutral solute with a caller-supplied Sechenov constant, in L/mol.
    pub fn neutral(name: &str, ks: f64, dsigma_dc: f64, provenance: Provenance) -> Self {
        Self {
            name: name.into(),
            model: SaltingModel::Neutral { ks },
            dsigma_dc,
            provenance,
        }
    }

    /// A **placeholder** for xylitol: `Ks = 0.10` L/mol and a small positive
    /// surface-tension increment. Tagged [`Provenance::Estimate`]; calibrate it.
    pub fn xylitol_estimate() -> Self {
        Self::neutral("xylitol", 0.10, 0.3e-3, Provenance::Estimate)
    }
}

/// The CO2 gas parameter `h_G(T)`, in L/mol.
pub fn gas_parameter(t_k: f64) -> f64 {
    H_G0_CO2 + H_T_CO2 * (t_k - T_REF)
}

impl Solute {
    /// This solute's contribution to `log10(S0/S)` per mol/L at `t_k`.
    pub fn sechenov_per_molar(&self, t_k: f64) -> f64 {
        match self.model {
            SaltingModel::Ion { h_i } => h_i + gas_parameter(t_k),
            SaltingModel::Neutral { ks } => ks,
        }
    }
}

/// A dissolved-solute recipe: species with their molar concentrations.
#[derive(Clone, PartialEq, Debug, Default)]
pub struct Formulation {
    solutes: Vec<(Solute, f64)>,
    /// The viscosity of the solution relative to water at the same
    /// temperature. Polyol syrups slow CO2 diffusion through this.
    pub viscosity_ratio: f64,
}

impl Formulation {
    /// Plain water: no solutes, relative viscosity 1.
    pub fn water() -> Self {
        Self {
            solutes: Vec::new(),
            viscosity_ratio: 1.0,
        }
    }

    /// Adds `molar` mol/L of `solute`.
    pub fn with(mut self, solute: Solute, molar: f64) -> Result<Self> {
        non_negative(molar, "solute concentration")?;
        self.solutes.push((solute, molar));
        Ok(self)
    }

    /// Sets the viscosity of the solution relative to water.
    pub fn with_viscosity_ratio(mut self, ratio: f64) -> Result<Self> {
        positive(ratio, "relative viscosity")?;
        self.viscosity_ratio = ratio;
        Ok(self)
    }

    /// The listed solutes and their concentrations.
    pub fn solutes(&self) -> &[(Solute, f64)] {
        &self.solutes
    }

    /// `log10(S0/S)` for this formulation at `t_k`: zero for pure water,
    /// growing with every solute that salts CO2 out.
    pub fn sechenov_log10(&self, t_k: f64) -> f64 {
        self.solutes
            .iter()
            .map(|(s, c)| s.sechenov_per_molar(t_k) * c)
            .sum()
    }

    /// The ratio `S/S0` of CO2 solubility in this solution to that in water.
    pub fn solubility_ratio(&self, t_k: f64) -> f64 {
        math::powf(10.0, -self.sechenov_log10(t_k))
    }

    /// The fractional loss of CO2 solubility versus pure water, in `[0, 1)`
    /// for a salting-out solution.
    pub fn solubility_loss(&self, t_k: f64) -> f64 {
        1.0 - self.solubility_ratio(t_k)
    }

    /// The Henry solubility of CO2 in this solution, in mol/(m^3 Pa).
    pub fn henry_mol_per_m3_pa(&self, t_k: f64) -> Result<f64> {
        Ok(crate::henry::henry_pure_mol_per_m3_pa(t_k)? * self.solubility_ratio(t_k))
    }

    /// The equilibrium dissolved CO2, in mol/m^3, at a given CO2 partial
    /// pressure in pascals.
    pub fn saturation_mol_per_m3(&self, t_k: f64, p_co2_pa: f64) -> Result<f64> {
        let phi = crate::henry::fugacity_coefficient(t_k, p_co2_pa)?;
        Ok(self.henry_mol_per_m3_pa(t_k)? * phi * p_co2_pa)
    }

    /// The CO2 partial pressure, in pascals, in equilibrium with a dissolved
    /// concentration `c` in mol/m^3. This is the pressure a sealed container
    /// holds at that carbonation.
    pub fn equilibrium_pressure(&self, t_k: f64, c_mol_per_m3: f64) -> Result<f64> {
        non_negative(c_mol_per_m3, "dissolved CO2")?;
        let kh = self.henry_mol_per_m3_pa(t_k)?;
        // Fixed-point on the fugacity coefficient, which varies slowly.
        let mut p = c_mol_per_m3 / kh;
        for _ in 0..6 {
            if p <= 0.0 {
                return Ok(0.0);
            }
            let phi = crate::henry::fugacity_coefficient(t_k, p)?;
            p = c_mol_per_m3 / (kh * phi);
        }
        Ok(p)
    }

    /// The surface tension of the solution, in N/m, from the water value at
    /// `t_k` plus the linear solute increments.
    pub fn surface_tension(&self, t_k: f64) -> f64 {
        use tpt_fluids_core::surface::SurfaceTensionModel;
        let base = SurfaceTensionModel::water().surface_tension(t_k).value();
        base + self
            .solutes
            .iter()
            .map(|(s, c)| s.dsigma_dc * c)
            .sum::<f64>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_water_is_unchanged() {
        let f = Formulation::water();
        assert_eq!(f.solubility_ratio(298.15), 1.0);
        assert_eq!(f.solubility_loss(298.15), 0.0);
    }

    #[test]
    fn one_molar_nacl_cuts_solubility_by_about_a_fifth() {
        let f = Formulation::water()
            .with(Solute::sodium(), 1.0)
            .unwrap()
            .with(Solute::chloride(), 1.0)
            .unwrap();
        // Weisenberger-Schumpe at 25 C: log10 = 0.1143 + 0.0318 - 2 * 0.0172.
        assert!((f.sechenov_log10(298.15) - 0.1117).abs() < 1e-9);
        let loss = f.solubility_loss(298.15);
        assert!((0.20..0.26).contains(&loss), "{loss}");
    }

    #[test]
    fn every_added_solute_lowers_solubility() {
        let mut f = Formulation::water();
        let mut last = 1.0;
        for s in [
            Solute::magnesium(),
            Solute::calcium(),
            Solute::xylitol_estimate(),
        ] {
            f = f.with(s, 0.05).unwrap();
            let r = f.solubility_ratio(277.15);
            assert!(r < last);
            last = r;
        }
    }

    #[test]
    fn xylitol_is_flagged_as_an_estimate() {
        assert_eq!(Solute::xylitol_estimate().provenance, Provenance::Estimate);
    }

    #[test]
    fn equilibrium_pressure_inverts_saturation() {
        let f = Formulation::water()
            .with(Solute::xylitol_estimate(), 0.5)
            .unwrap();
        let c = f.saturation_mol_per_m3(277.15, 2.0 * 101_325.0).unwrap();
        let p = f.equilibrium_pressure(277.15, c).unwrap();
        assert!(
            (p - 2.0 * 101_325.0).abs() / (2.0 * 101_325.0) < 1e-6,
            "{p}"
        );
    }

    #[test]
    fn salts_raise_surface_tension() {
        let w = Formulation::water().surface_tension(293.15);
        let s = Formulation::water()
            .with(Solute::sodium(), 0.5)
            .unwrap()
            .surface_tension(293.15);
        assert!(s > w);
    }

    #[test]
    fn rejects_negative_concentration() {
        assert!(Formulation::water().with(Solute::sodium(), -1.0).is_err());
    }
}

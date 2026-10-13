//! Carbonate speciation, pH buffers, and the rate CO2 is released from them.
//!
//! Dissolved CO2 is not all molecular. At the pH of a buffered drink part of
//! it sits as bicarbonate, which cannot nucleate a bubble and can only leave
//! after it has been protonated and dehydrated. This module computes how the
//! total splits for a given buffer and how fast the reservoir drains.
//!
//! ```text
//! [CO2*] = kH p_CO2
//! [H+][HCO3-]   = K1 [CO2*]
//! [H+][CO3 2-]  = K2 [HCO3-]
//! [H+] + C_net = [OH-] + [HCO3-] + 2[CO3 2-] + sum_buffers z_bar(pH) C_b
//! ```
//!
//! `C_net` is the net excess of strong-base cations (for example from mineral
//! carbonates or citrates) over strong-acid anions. Activity coefficients are
//! taken as one, which is adequate below roughly 0.1 mol/L ionic strength and a
//! screening approximation above it.
//!
//! `K1` and `K2` are the Millero (1979) freshwater expressions; the hydration
//! kinetics use `k_h = 0.039 /s` at 25 C with a 70 kJ/mol activation energy.

use tpt_fluids_core::math;

use crate::consts::{MOLAR_GAS_CONSTANT, T_REF};
use crate::error::{non_negative, positive, BevError, Result};

/// A polyprotic buffer acid `H_nA` with up to three dissociation steps.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Buffer {
    /// The total analytical concentration, in mol/L.
    pub total_molar: f64,
    /// The acid dissociation constants `pKa`, in order. Only the first
    /// `steps` are used.
    pub pka: [f64; 3],
    /// The number of protons the acid can lose, 1 to 3.
    pub steps: usize,
}

impl Buffer {
    /// Builds a buffer from its total concentration and `pKa` values.
    pub fn new(total_molar: f64, pka: &[f64]) -> Result<Self> {
        non_negative(total_molar, "buffer concentration")?;
        if pka.is_empty() || pka.len() > 3 {
            return Err(BevError::OutOfRange("number of dissociation steps"));
        }
        let mut arr = [0.0; 3];
        arr[..pka.len()].copy_from_slice(pka);
        Ok(Self {
            total_molar,
            pka: arr,
            steps: pka.len(),
        })
    }

    /// Citric acid (`pKa` 3.13, 4.76, 6.40).
    pub fn citric(total_molar: f64) -> Result<Self> {
        Self::new(total_molar, &[3.13, 4.76, 6.40])
    }

    /// Phosphoric acid (`pKa` 2.15, 7.20, 12.35).
    pub fn phosphoric(total_molar: f64) -> Result<Self> {
        Self::new(total_molar, &[2.15, 7.20, 12.35])
    }

    /// Malic acid (`pKa` 3.40, 5.20).
    pub fn malic(total_molar: f64) -> Result<Self> {
        Self::new(total_molar, &[3.40, 5.20])
    }

    /// The mean negative charge per acid molecule at hydrogen-ion
    /// concentration `h` (mol/L).
    fn mean_charge(&self, h: f64) -> f64 {
        // Relative weights of H_nA, H_(n-1)A-, ...: w_0 = 1, w_k = w_(k-1) Ka_k / h.
        let mut w = 1.0;
        let mut sum_w = 1.0;
        let mut sum_zw = 0.0;
        for k in 0..self.steps {
            let ka = math::powf(10.0, -self.pka[k]);
            w *= ka / h;
            sum_w += w;
            sum_zw += (k as f64 + 1.0) * w;
        }
        sum_zw / sum_w
    }
}

/// The first dissociation constant of carbonic acid, `K1`, in mol/L.
pub fn k1(t_k: f64) -> f64 {
    let pk = -126.343_21 + 6320.813 / t_k + 19.568_224 * math::ln(t_k);
    math::powf(10.0, -pk)
}

/// The second dissociation constant of carbonic acid, `K2`, in mol/L.
pub fn k2(t_k: f64) -> f64 {
    let pk = -90.183_33 + 5143.692 / t_k + 14.613_358 * math::ln(t_k);
    math::powf(10.0, -pk)
}

/// The ion product of water, `Kw`, in (mol/L)^2.
pub fn kw(t_k: f64) -> f64 {
    let pk = 4471.33 / t_k + 0.017_060 * t_k - 6.0846;
    math::powf(10.0, -pk)
}

/// The equilibrium state of a carbonated, buffered solution.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CarbonateState {
    /// The pH.
    pub ph: f64,
    /// Molecular dissolved CO2, in mol/L.
    pub co2_molar: f64,
    /// Bicarbonate, in mol/L.
    pub hco3_molar: f64,
    /// Carbonate, in mol/L.
    pub co3_molar: f64,
}

impl CarbonateState {
    /// Total dissolved inorganic carbon, in mol/L.
    pub fn dic_molar(&self) -> f64 {
        self.co2_molar + self.hco3_molar + self.co3_molar
    }

    /// The fraction of dissolved inorganic carbon that is molecular CO2 and so
    /// able to supersaturate and nucleate bubbles.
    pub fn molecular_fraction(&self) -> f64 {
        self.co2_molar / self.dic_molar()
    }

    /// The bicarbonate-to-CO2 ratio `[HCO3-]/[CO2*]`.
    pub fn bicarbonate_ratio(&self) -> f64 {
        self.hco3_molar / self.co2_molar
    }
}

/// Solves the speciation at fixed molecular CO2.
///
/// - `co2_molar`: molecular dissolved CO2 `[CO2*]`, in mol/L.
/// - `net_strong_cation_molar`: net strong-base equivalents (cations minus
///   strong-acid anions), in mol/L; zero for an unbuffered drink.
/// - `buffers`: the buffer acids present.
pub fn speciate(
    t_k: f64,
    co2_molar: f64,
    net_strong_cation_molar: f64,
    buffers: &[Buffer],
) -> Result<CarbonateState> {
    positive(t_k, "temperature")?;
    positive(co2_molar, "dissolved CO2")?;
    let (k1, k2, kw) = (k1(t_k), k2(t_k), kw(t_k));
    let charge_balance = |ph: f64| -> f64 {
        let h = math::powf(10.0, -ph);
        let hco3 = k1 * co2_molar / h;
        let co3 = k2 * hco3 / h;
        let buffer: f64 = buffers
            .iter()
            .map(|b| b.mean_charge(h) * b.total_molar)
            .sum();
        // Positive minus negative charge; falls monotonically as pH rises.
        h + net_strong_cation_molar - kw / h - hco3 - 2.0 * co3 - buffer
    };
    let (mut lo, mut hi) = (0.0_f64, 14.0_f64);
    if charge_balance(lo) < 0.0 || charge_balance(hi) > 0.0 {
        return Err(BevError::NoConvergence("pH bracket"));
    }
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if charge_balance(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let ph = 0.5 * (lo + hi);
    let h = math::powf(10.0, -ph);
    let hco3 = k1 * co2_molar / h;
    Ok(CarbonateState {
        ph,
        co2_molar,
        hco3_molar: hco3,
        co3_molar: k2 * hco3 / h,
    })
}

/// The rate constant of CO2 hydration, `CO2 + H2O -> H2CO3`, in 1/s.
pub fn hydration_rate(t_k: f64) -> f64 {
    const K25: f64 = 0.039;
    const EA: f64 = 70_000.0;
    K25 * math::exp(-EA / MOLAR_GAS_CONSTANT * (1.0 / t_k - 1.0 / T_REF))
}

/// The pseudo-first-order rate, in 1/s, at which the bicarbonate reservoir is
/// converted back to molecular CO2 at a given pH.
///
/// Bicarbonate is protonated to `H2CO3` (true `pKa` 3.45), which dehydrates at
/// `k_d = k_h / K_h` with `K_h = [H2CO3]/[CO2] = 1.7e-3`:
///
/// ```text
/// k = k_d [H+] / Ka(H2CO3)
/// ```
///
/// This is why a drink buffered near neutral pH holds its CO2 as bicarbonate
/// and releases it slowly, while an acidic one releases it within a second or
/// two. It is an engineering estimate.
pub fn bicarbonate_release_rate(t_k: f64, ph: f64) -> f64 {
    const K_HYDRATION_EQ: f64 = 1.7e-3;
    const KA_H2CO3: f64 = 3.55e-4;
    let k_dehydration = hydration_rate(t_k) / K_HYDRATION_EQ;
    k_dehydration * math::powf(10.0, -ph) / KA_H2CO3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dissociation_constants_match_reference_at_25c() {
        let pk1 = -math::log10(k1(298.15));
        let pk2 = -math::log10(k2(298.15));
        let pkw = -math::log10(kw(298.15));
        assert!((pk1 - 6.35).abs() < 0.03, "{pk1}");
        assert!((pk2 - 10.33).abs() < 0.03, "{pk2}");
        assert!((pkw - 14.0).abs() < 0.03, "{pkw}");
    }

    #[test]
    fn saturated_pure_water_is_acidic_at_about_ph_3_9() {
        // 1 atm CO2 at 25 C is about 33 mmol/L.
        let s = speciate(298.15, 0.033, 0.0, &[]).unwrap();
        assert!((3.8..4.0).contains(&s.ph), "{}", s.ph);
        assert!(s.molecular_fraction() > 0.99);
    }

    #[test]
    fn a_buffer_with_base_raises_ph_and_bicarbonate() {
        let bare = speciate(277.15, 0.15, 0.0, &[]).unwrap();
        let citrate = Buffer::citric(0.01).unwrap();
        // Partly neutralised: 20 mmol/L of strong-base cation against 10 mmol/L citrate.
        let buffered = speciate(277.15, 0.15, 0.02, &[citrate]).unwrap();
        assert!(buffered.ph > bare.ph + 1.0);
        assert!(buffered.hco3_molar > bare.hco3_molar);
        assert!(buffered.molecular_fraction() < bare.molecular_fraction());
    }

    #[test]
    fn charge_balance_closes() {
        let citrate = Buffer::citric(0.005).unwrap();
        let s = speciate(293.15, 0.1, 0.012, &[citrate]).unwrap();
        let h = math::powf(10.0, -s.ph);
        let neg = kw(293.15) / h
            + s.hco3_molar
            + 2.0 * s.co3_molar
            + citrate.mean_charge(h) * citrate.total_molar;
        assert!((h + 0.012 - neg).abs() < 1e-9, "{}", h + 0.012 - neg);
    }

    #[test]
    fn acidic_drinks_release_fast_and_neutral_ones_slowly() {
        let acidic = bicarbonate_release_rate(277.15, 3.0);
        let neutral = bicarbonate_release_rate(277.15, 7.0);
        assert!(acidic > 1.0, "{acidic}");
        assert!(neutral < 0.1, "{neutral}");
        assert!(acidic / neutral > 1.0e3);
    }

    #[test]
    fn hydration_is_faster_when_warm() {
        assert!(hydration_rate(310.0) > hydration_rate(280.0));
        assert!((hydration_rate(298.15) - 0.039).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_buffers() {
        assert!(Buffer::new(0.01, &[]).is_err());
        assert!(Buffer::new(-1.0, &[4.0]).is_err());
    }
}

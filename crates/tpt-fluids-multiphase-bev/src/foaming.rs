//! Screening for foam-over versus a flat drink.
//!
//! The two failure modes are opposite ends of one dial, the number of active
//! nucleation sites:
//!
//! - **Foam-over.** Too many sites, or sites that grow bubbles fast, release
//!   gas faster than the headspace and pour can contain it. The gas-release
//!   rate on opening exceeds the foam-stability threshold.
//! - **Flat.** Too few sites or no bicarbonate reservoir, so the drink loses
//!   its carbonation within the shelf life before it is poured.
//!
//! The regime is read from the release trace, not from one instant. Both
//! limits are inputs to [`Thresholds`] and the defaults are **uncalibrated**:
//! they are placeholders chosen to make the ordering of regimes visible. Set
//! them from bench data before trusting an absolute verdict.

use crate::error::{positive, Result};
use crate::kinetics::Trace;

/// The screening verdict.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Regime {
    /// Carbonation is lost within the shelf life.
    Flat,
    /// Carbonation survives the shelf life and the pour is controlled.
    Stable,
    /// The initial release rate would overflow the glass.
    FoamOver,
}

/// The limits that separate the regimes. All in SI units.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Thresholds {
    /// The peak volumetric release rate above which the pour foams over, in
    /// mol/(m^3 s). **Uncalibrated.**
    pub foam_over_rate: f64,
    /// The dissolved concentration that counts as "still carbonated" at the
    /// end of shelf life, in mol/m^3. **Uncalibrated.**
    pub flat_below: f64,
    /// The shelf life, in seconds.
    pub shelf_life: f64,
}

impl Default for Thresholds {
    /// Placeholder limits: foam over at a peak release of 0.5 mol/(m^3 s),
    /// flat once below 1 volume of CO2 after a 90-day shelf life.
    fn default() -> Self {
        Self {
            foam_over_rate: 0.5,
            flat_below: crate::henry::volumes_to_mol_per_m3(1.0),
            shelf_life: 90.0 * 86_400.0,
        }
    }
}

/// The screening result, with margins so a near miss is visible.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Assessment {
    /// The verdict.
    pub regime: Regime,
    /// The peak volumetric release rate on opening, in mol/(m^3 s).
    pub peak_release: f64,
    /// The dissolved concentration at the end of shelf life, in mol/m^3.
    pub dissolved_at_shelf_life: f64,
}

impl Assessment {
    /// Reads the verdict from a release trace that covers the shelf life.
    pub fn from_trace(trace: &Trace, th: &Thresholds) -> Result<Self> {
        positive(th.foam_over_rate, "foam-over threshold")?;
        positive(th.shelf_life, "shelf life")?;
        let peak = trace.release_rate.iter().copied().fold(0.0, f64::max);
        let end = trace.dissolved.last().copied().unwrap_or(0.0);
        let regime = if peak > th.foam_over_rate {
            Regime::FoamOver
        } else if end < th.flat_below {
            Regime::Flat
        } else {
            Regime::Stable
        };
        Ok(Self {
            regime,
            peak_release: peak,
            dissolved_at_shelf_life: end,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kinetics::{integrate, Release};

    fn run(site_rate: f64, k_bicarb: f64) -> Trace {
        let rel = Release {
            c_eq: 0.0,
            k_surface: 1.0e-8,
            k_bicarbonate: k_bicarb,
            bicarbonate_reservoir: 0.0,
        };
        integrate(100.0, rel, 90.0 * 86_400.0, 4000, |c| {
            if c > 0.0 {
                site_rate
            } else {
                0.0
            }
        })
        .unwrap()
    }

    #[test]
    fn many_sites_foam_over() {
        let t = run(5.0, 0.0);
        let a = Assessment::from_trace(&t, &Thresholds::default()).unwrap();
        assert_eq!(a.regime, Regime::FoamOver);
    }

    #[test]
    fn a_quiet_drink_with_no_reservoir_goes_flat() {
        // Surface loss alone empties the drink well below one volume.
        let rel = Release {
            c_eq: 0.0,
            k_surface: 1.0e-6,
            k_bicarbonate: 0.0,
            bicarbonate_reservoir: 0.0,
        };
        let t = integrate(100.0, rel, 90.0 * 86_400.0, 4000, |_| 0.0).unwrap();
        let a = Assessment::from_trace(&t, &Thresholds::default()).unwrap();
        assert_eq!(a.regime, Regime::Flat);
    }

    #[test]
    fn slow_surface_loss_alone_stays_stable_over_shelf_life() {
        // Barely any surface loss: still above one volume at 90 days.
        let rel = Release {
            c_eq: 0.0,
            k_surface: 1.0e-8,
            k_bicarbonate: 0.0,
            bicarbonate_reservoir: 0.0,
        };
        let t = integrate(100.0, rel, 90.0 * 86_400.0, 4000, |_| 0.0).unwrap();
        let a = Assessment::from_trace(&t, &Thresholds::default()).unwrap();
        assert_eq!(a.regime, Regime::Stable);
    }

    #[test]
    fn a_reservoir_holds_a_quiet_drink_stable() {
        // Bicarbonate refill balances surface loss at c = 50 mol/m^3 > 1 volume.
        let rel = Release {
            c_eq: 0.0,
            k_surface: 1.0e-6,
            k_bicarbonate: 1.0e-6,
            bicarbonate_reservoir: 50.0,
        };
        let t = integrate(100.0, rel, 90.0 * 86_400.0, 4000, |_| 0.0).unwrap();
        let a = Assessment::from_trace(&t, &Thresholds::default()).unwrap();
        assert_eq!(a.regime, Regime::Stable);
    }

    #[test]
    fn peak_release_rises_with_site_activity() {
        // Foam-over is decided by the instantaneous peak, which is monotone in
        // site activity. The shelf-life verdict is not: a site-rich drink
        // depletes itself and can end up flat, so only the peak is asserted.
        let mut last = 0.0;
        for rate in [0.0, 0.1, 1.0, 10.0] {
            let a = Assessment::from_trace(&run(rate, 0.0), &Thresholds::default()).unwrap();
            assert!(a.peak_release >= last);
            last = a.peak_release;
        }
    }

    #[test]
    fn extreme_site_activity_foams_over_and_moderate_goes_flat() {
        let th = Thresholds::default();
        assert_eq!(
            Assessment::from_trace(&run(10.0, 0.0), &th).unwrap().regime,
            Regime::FoamOver
        );
        assert_eq!(
            Assessment::from_trace(&run(0.1, 0.0), &th).unwrap().regime,
            Regime::Flat
        );
    }

    #[test]
    fn rejects_nonpositive_thresholds() {
        let t = run(0.0, 0.0);
        let th = Thresholds {
            foam_over_rate: 0.0,
            ..Thresholds::default()
        };
        assert!(Assessment::from_trace(&t, &th).is_err());
    }
}

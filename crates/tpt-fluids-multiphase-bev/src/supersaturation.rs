//! Supersaturation of dissolved CO2 and the driving pressure for bubbles.
//!
//! A drink at dissolved concentration `c` is in equilibrium with a CO2
//! partial pressure `p_eq = c / kH`. Once the container is opened the liquid
//! sits at ambient pressure, and the excess is the thermodynamic driving force
//! for every bubble process:
//!
//! ```text
//! S  = p_eq / p_ambient
//! dP = p_eq + p_vap - p_ambient
//! ```
//!
//! `dP` is the overpressure that a gas embryo of critical radius
//! `r* = 2 sigma / dP` has to balance, and it enters the nucleation barrier
//! squared, so the rate is extremely sensitive to it.

use tpt_fluids_core::fluids::water_vapour_pressure;

use crate::error::{non_negative, positive, Result};
use crate::salting_out::Formulation;

/// The supersaturation state of an opened drink.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Supersaturation {
    /// The dissolved CO2, in mol/m^3.
    pub dissolved: f64,
    /// The CO2 partial pressure the liquid is in equilibrium with, in Pa.
    pub equilibrium_pressure: f64,
    /// The ambient (liquid) pressure, in Pa.
    pub ambient_pressure: f64,
    /// The water vapour pressure at the drink temperature, in Pa.
    pub vapour_pressure: f64,
}

impl Supersaturation {
    /// Evaluates the state of `formulation` at `t_k` holding `dissolved`
    /// mol/m^3 of CO2 at ambient pressure `p_ambient_pa`.
    pub fn new(
        formulation: &Formulation,
        t_k: f64,
        dissolved_mol_per_m3: f64,
        p_ambient_pa: f64,
    ) -> Result<Self> {
        non_negative(dissolved_mol_per_m3, "dissolved CO2")?;
        positive(p_ambient_pa, "ambient pressure")?;
        Ok(Self {
            dissolved: dissolved_mol_per_m3,
            equilibrium_pressure: formulation.equilibrium_pressure(t_k, dissolved_mol_per_m3)?,
            ambient_pressure: p_ambient_pa,
            vapour_pressure: water_vapour_pressure(t_k),
        })
    }

    /// The supersaturation ratio `p_eq / p_ambient`; above one the liquid is
    /// supersaturated.
    pub fn ratio(&self) -> f64 {
        self.equilibrium_pressure / self.ambient_pressure
    }

    /// The bubble driving overpressure `dP`, in pascals, floored at zero.
    pub fn driving_pressure(&self) -> f64 {
        (self.equilibrium_pressure + self.vapour_pressure - self.ambient_pressure).max(0.0)
    }

    /// The critical embryo radius `2 sigma / dP`, in metres, for a surface
    /// tension `sigma_n_per_m`. Infinite when there is no driving pressure.
    pub fn critical_radius(&self, sigma_n_per_m: f64) -> f64 {
        let dp = self.driving_pressure();
        if dp <= 0.0 {
            f64::INFINITY
        } else {
            2.0 * sigma_n_per_m / dp
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::ATM;
    use crate::henry::volumes_to_mol_per_m3;

    fn state(volumes: f64, t_k: f64) -> Supersaturation {
        Supersaturation::new(
            &Formulation::water(),
            t_k,
            volumes_to_mol_per_m3(volumes),
            ATM,
        )
        .unwrap()
    }

    #[test]
    fn a_cold_carbonated_drink_is_supersaturated_by_over_one_atmosphere() {
        let s = state(3.5, 277.15);
        assert!((2.0..2.8).contains(&s.ratio()), "{}", s.ratio());
        assert!(s.driving_pressure() > 0.9 * ATM);
    }

    #[test]
    fn flat_water_has_no_driving_pressure() {
        let s = state(0.0, 277.15);
        assert_eq!(s.driving_pressure(), 0.0);
        assert!(s.critical_radius(0.073).is_infinite());
    }

    #[test]
    fn critical_radius_shrinks_with_carbonation() {
        let low = state(2.0, 277.15).critical_radius(0.074);
        let high = state(4.0, 277.15).critical_radius(0.074);
        assert!(high < low);
        // Micron scale at a bottling level.
        assert!((0.1e-6..5.0e-6).contains(&high), "{high}");
    }

    #[test]
    fn warmer_drinks_are_more_supersaturated_at_equal_content() {
        assert!(state(3.0, 295.15).ratio() > state(3.0, 277.15).ratio());
    }
}

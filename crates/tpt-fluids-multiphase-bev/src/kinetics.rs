//! Time evolution of dissolved CO2 in an opened drink.
//!
//! Dissolved concentration `c` (mol/m^3) falls as gas leaves through two
//! routes: bubbles nucleated at the active sites, and the free surface. The
//! bicarbonate reservoir refills `c` at a rate set by its release constant:
//!
//! ```text
//! dc/dt = -Q_sites(c) - k_surf (c - c_eq) + k_HCO3 * c_HCO3
//! ```
//!
//! `Q_sites` is the volumetric gas loss from active sites, each bubble taking
//! its mass out of the liquid. It is held in the callback so the caller can
//! supply the nucleation model for its own formulation. The integrator is a
//! fixed-step classical Runge-Kutta, with an explicit step count so results
//! are reproducible.

use alloc::vec::Vec;

use crate::error::{non_negative, positive, Result};

/// A time series of the state, sampled at every step.
#[derive(Clone, PartialEq, Debug)]
pub struct Trace {
    /// Sample times, in seconds.
    pub time: Vec<f64>,
    /// Dissolved CO2 at each sample, in mol/m^3.
    pub dissolved: Vec<f64>,
    /// Volumetric gas-release rate at each sample, in mol/(m^3 s).
    pub release_rate: Vec<f64>,
}

impl Trace {
    /// The first time the dissolved concentration falls to `target`, in
    /// seconds, found by linear interpolation. `None` if it never does.
    pub fn time_to_reach(&self, target: f64) -> Option<f64> {
        for i in 1..self.dissolved.len() {
            let (c0, c1) = (self.dissolved[i - 1], self.dissolved[i]);
            if c1 <= target && c0 > target {
                let f = (c0 - target) / (c0 - c1);
                return Some(self.time[i - 1] + f * (self.time[i] - self.time[i - 1]));
            }
        }
        None
    }
}

/// The parameters of the release model.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Release {
    /// The equilibrium concentration at ambient pressure, in mol/m^3.
    pub c_eq: f64,
    /// The free-surface mass-transfer rate `k_L a / V`, in 1/s.
    pub k_surface: f64,
    /// The bicarbonate release rate, in 1/s.
    pub k_bicarbonate: f64,
    /// The bicarbonate reservoir, in mol/m^3.
    pub bicarbonate_reservoir: f64,
}

/// Integrates the release from `c0` for `t_end` seconds in `steps` steps.
///
/// `sites` returns the volumetric gas-release rate from bubbles, in
/// mol/(m^3 s), for a given dissolved concentration. It must be non-negative.
pub fn integrate<F>(
    c0: f64,
    release: Release,
    t_end: f64,
    steps: usize,
    mut sites: F,
) -> Result<Trace>
where
    F: FnMut(f64) -> f64,
{
    non_negative(c0, "initial dissolved CO2")?;
    positive(t_end, "end time")?;
    positive(steps as f64, "step count")?;
    let h = t_end / steps as f64;
    let mut trace = Trace {
        time: Vec::with_capacity(steps + 1),
        dissolved: Vec::with_capacity(steps + 1),
        release_rate: Vec::with_capacity(steps + 1),
    };
    let mut c = c0;
    let mut t = 0.0;
    let rhs = |c: f64, sites: &mut F| -> (f64, f64) {
        let q = sites(c).max(0.0);
        let surf = release.k_surface * (c - release.c_eq);
        let bic = release.k_bicarbonate * release.bicarbonate_reservoir;
        (-(q + surf) + bic, q)
    };
    for _ in 0..=steps {
        let (_, q) = rhs(c, &mut sites);
        trace.time.push(t);
        trace.dissolved.push(c);
        trace.release_rate.push(q);
        if t >= t_end {
            break;
        }
        let (k1, _) = rhs(c, &mut sites);
        let (k2, _) = rhs(c + 0.5 * h * k1, &mut sites);
        let (k3, _) = rhs(c + 0.5 * h * k2, &mut sites);
        let (k4, _) = rhs(c + h * k3, &mut sites);
        c = (c + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)).max(0.0);
        t += h;
    }
    Ok(trace)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quiet() -> Release {
        Release {
            c_eq: 0.0,
            k_surface: 0.0,
            k_bicarbonate: 0.0,
            bicarbonate_reservoir: 0.0,
        }
    }

    #[test]
    fn with_no_loss_terms_concentration_is_constant() {
        let t = integrate(50.0, quiet(), 100.0, 100, |_| 0.0).unwrap();
        assert!(t.dissolved.iter().all(|&c| (c - 50.0).abs() < 1e-12));
    }

    #[test]
    fn surface_loss_decays_exponentially() {
        let rel = Release {
            k_surface: 0.01,
            ..quiet()
        };
        let t = integrate(50.0, rel, 100.0, 1000, |_| 0.0).unwrap();
        let last = *t.dissolved.last().unwrap();
        let expect = 50.0 * tpt_fluids_core::math::exp(-1.0);
        assert!((last - expect).abs() / expect < 1e-6, "{last} vs {expect}");
    }

    #[test]
    fn gas_removed_from_sites_equals_loss_from_liquid() {
        // With a constant site rate Q, dc/dt = -Q so c(T) = c0 - Q T.
        let t = integrate(50.0, quiet(), 10.0, 100, |_| 2.0).unwrap();
        let last = *t.dissolved.last().unwrap();
        assert!((last - (50.0 - 20.0)).abs() < 1e-9, "{last}");
        assert!(t.release_rate.iter().all(|&q| (q - 2.0).abs() < 1e-12));
    }

    #[test]
    fn refinement_converges() {
        let rel = Release {
            k_surface: 0.05,
            c_eq: 5.0,
            ..quiet()
        };
        let coarse = integrate(50.0, rel, 60.0, 60, |c| 0.1 * c).unwrap();
        let fine = integrate(50.0, rel, 60.0, 600, |c| 0.1 * c).unwrap();
        let a = *coarse.dissolved.last().unwrap();
        let b = *fine.dissolved.last().unwrap();
        assert!((a - b).abs() < 1e-3, "{a} vs {b}");
    }

    #[test]
    fn half_life_is_found_by_interpolation() {
        let rel = Release {
            k_surface: 0.1,
            ..quiet()
        };
        let t = integrate(100.0, rel, 50.0, 5000, |_| 0.0).unwrap();
        let half = t.time_to_reach(50.0).unwrap();
        assert!(
            (half - core::f64::consts::LN_2 / 0.1).abs() < 1e-3,
            "{half}"
        );
    }

    #[test]
    fn rejects_negative_start() {
        assert!(integrate(-1.0, quiet(), 1.0, 10, |_| 0.0).is_err());
    }
}

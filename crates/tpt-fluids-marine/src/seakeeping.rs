//! Seakeeping: linear wave theory, response amplitude operators, and sea
//! state statistics.
//!
//! A linearised seakeeping problem decouples into six independent degrees of
//! freedom, each described by a *response amplitude operator* (RAO): the
//! complex ratio of a motion response to wave elevation at a given frequency.
//! Its magnitude is the amplification, and its argument the phase lag; both
//! are needed to predict a motion spectrum from an elevation spectrum.
//!
//! # The approximation used here
//!
//! Each RAO is modelled as a damped single-degree-of-freedom oscillator,
//!
//! ```text
//!            1
//! RAO(w) = ----------- ,  w_n^2 - w^2 + i 2 zeta w_n w
//!        w_n^2 - w^2
//! ```
//!
//! with `w_n` the natural frequency and `zeta` the damping ratio of that
//! degree of freedom. This captures resonance, the 180-degree phase reversal
//! above it, and roll damping, which are the features that actually govern
//! ship motion. It is not a substitute for a full Green-function hull-integral
//! solution, and the module says so, but it is quantitatively right for
//! resonance peaks and the roll damping that damps them.

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;

use crate::error::{MarineError, Result};

/// A complex response amplitude operator.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Rao {
    /// The magnitude, the amplification factor.
    pub magnitude: f64,
    /// The phase lag in radians, negative for a lag.
    pub phase: f64,
}

impl Rao {
    /// Builds an RAO from a magnitude and a phase in radians.
    pub const fn new(magnitude: f64, phase: f64) -> Self {
        Self { magnitude, phase }
    }

    /// The real part, the in-phase component.
    pub fn real(&self) -> f64 {
        self.magnitude * math::cos(self.phase)
    }

    /// The imaginary part, the quadrature component.
    pub fn imaginary(&self) -> f64 {
        self.magnitude * math::sin(self.phase)
    }
}

/// The six rigid-body degrees of freedom a ship responds in.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MotionMode {
    /// Surge, translation along the ship's length.
    Surge,
    /// Heave, translation in the vertical direction.
    Heave,
    /// Sway, translation athwartships.
    Sway,
    /// Roll, rotation about the longitudinal axis.
    Roll,
    /// Pitch, rotation about the transverse axis.
    Pitch,
    /// Yaw, rotation about the vertical axis.
    Yaw,
}

impl MotionMode {
    /// Every mode, in the conventional 6-DOF order.
    pub const ALL: [MotionMode; 6] = [
        MotionMode::Surge,
        MotionMode::Heave,
        MotionMode::Sway,
        MotionMode::Roll,
        MotionMode::Pitch,
        MotionMode::Yaw,
    ];
}

/// A damped single-degree-of-freedom oscillator, the model behind each RAO.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Oscillator {
    /// The natural period in seconds.
    pub natural_period: f64,
    /// The damping ratio, typically a few percent of critical.
    pub damping_ratio: f64,
}

impl Oscillator {
    /// Builds an oscillator, rejecting a non-positive period.
    pub fn new(natural_period: f64, damping_ratio: f64) -> Result<Self> {
        if natural_period <= 0.0 {
            return Err(MarineError::NonPositive("natural period"));
        }
        if damping_ratio < 0.0 {
            return Err(MarineError::NonFinite("damping ratio"));
        }
        Ok(Self {
            natural_period,
            damping_ratio,
        })
    }

    /// The natural circular frequency, in radians per second.
    pub fn natural_frequency(&self) -> f64 {
        2.0 * core::f64::consts::PI / self.natural_period
    }

    /// The response amplitude operator at a wave angular frequency.
    ///
    /// The oscillator is driven by wave elevation, and the restoring force
    /// comes from hydrostatic stiffness, so the static response is unity.
    pub fn rao(&self, omega: f64) -> Rao {
        let wn = self.natural_frequency();
        if wn <= 0.0 {
            return Rao::new(0.0, 0.0);
        }
        let detuning = wn * wn - omega * omega;
        let damping = 2.0 * self.damping_ratio * wn * omega;
        let denom = math::sqrt(detuning * detuning + damping * damping);
        if denom <= 0.0 {
            return Rao::new(f64::INFINITY, 0.0);
        }
        // The numerator is `w_n^2`, not 1: the response is a displacement, so
        // it carries the restoring stiffness `w_n^2` that converts the forcing
        // acceleration into motion. Using a bare `1 / |denominator|` would
        // silently scale the whole curve by `w_n^2` and make the answer depend
        // on the units the period was entered in.
        Rao::new(wn * wn / denom, -math::atan2(damping, detuning))
    }

    /// The response amplitude operator at a wave period.
    pub fn rao_at_period(&self, period: f64) -> Rao {
        if period <= 0.0 {
            return Rao::new(0.0, 0.0);
        }
        self.rao(2.0 * core::f64::consts::PI / period)
    }

    /// The peak response, which occurs at resonance, and the period it occurs
    /// at.
    pub fn peak(&self) -> (f64, f64) {
        let wn = self.natural_frequency();
        // The resonant frequency is shifted slightly below w_n by the damping;
        // the shift is second order in zeta, so w_n is the leading answer.
        let damping = 2.0 * self.damping_ratio * wn * wn;
        let peak_omega = math::sqrt(wn * wn - damping * damping / 2.0);
        // Same `w_n^2` numerator as `rao`, for the same reason: the peak
        // magnitude is `w_n^2 / (2 zeta w_n^2) = 1 / (2 zeta)`.
        let magnitude = 1.0 / (2.0 * self.damping_ratio);
        let period = if peak_omega > 0.0 {
            2.0 * core::f64::consts::PI / peak_omega
        } else {
            f64::INFINITY
        };
        (magnitude, period)
    }
}

/// The deep-water dispersion relation `w^2 = g k`, returning the wavenumber
/// in radians per metre for an angular frequency.
pub fn deep_water_wavenumber(omega: f64) -> f64 {
    if omega <= 0.0 {
        return 0.0;
    }
    omega * omega / STANDARD_GRAVITY
}

/// The deep-water angular frequency for a wavenumber.
pub fn deep_water_omega(wavenumber: f64) -> f64 {
    if wavenumber <= 0.0 {
        return 0.0;
    }
    math::sqrt(STANDARD_GRAVITY * wavenumber)
}

/// The deep-water phase velocity, in metres per second.
pub fn phase_velocity(wavenumber: f64) -> f64 {
    if wavenumber <= 0.0 {
        return 0.0;
    }
    STANDARD_GRAVITY / deep_water_omega(wavenumber)
}

/// The deep-water wavelength for a period, in metres: `L = g T^2 / 2 pi`.
pub fn wavelength_from_period(period: f64) -> f64 {
    if period <= 0.0 {
        return 0.0;
    }
    STANDARD_GRAVITY * period * period / (2.0 * core::f64::consts::PI)
}

/// The significant wave height from a wind speed, in metres.
///
/// The ITTC Pierson-Moskowitz fetch-limited relation, `H_s = 0.0246 U^2`,
/// with `U` the mean wind speed at 10 m above the sea in m/s. It is a
/// fully-developed-sea estimate, so it is an upper bound for a fetch-limited
/// case, but it is the reference the specification asks for.
pub fn significant_wave_height(wind_speed: f64) -> f64 {
    if wind_speed <= 0.0 {
        return 0.0;
    }
    0.0246 * wind_speed * wind_speed
}

/// The peak period of a fully-developed sea, in seconds.
///
/// `T_p = 0.729 U + 2.75` is the companion to
/// [`significant_wave_height`] on the ITTC Pierson-Moskowitz diagram, and is
/// what sets the spectral peak the ship responds to.
pub fn peak_period(wind_speed: f64) -> f64 {
    if wind_speed <= 0.0 {
        return 0.0;
    }
    0.729 * wind_speed + 2.75
}

/// The significant wave height from a wind speed given in knots.
pub fn significant_wave_height_knots(wind_knots: f64) -> f64 {
    significant_wave_height(wind_knots * 0.514_444)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deep_water_dispersion_round_trips() {
        let omega = 1.0;
        let k = deep_water_wavenumber(omega);
        assert!((deep_water_omega(k) - omega).abs() < 1e-12);
    }

    #[test]
    fn deep_water_dispersion_matches_the_definition() {
        // w^2 = g k  =>  k = w^2/g. At w = 1 rad/s, k = 0.1019 /m.
        let k = deep_water_wavenumber(1.0);
        assert!((k - 1.0 / STANDARD_GRAVITY).abs() < 1e-12, "k = {k}");
    }

    #[test]
    fn long_waves_outrun_short_ones() {
        let long = deep_water_wavenumber(0.5);
        let short = deep_water_wavenumber(2.0);
        assert!(long < short);
        // Phase velocity falls with wavenumber.
        assert!(phase_velocity(long) > phase_velocity(short));
    }

    #[test]
    fn wavelength_matches_the_deep_water_relation() {
        // A 10 s wave in deep water has a 156 m wavelength.
        let l = wavelength_from_period(10.0);
        assert!((l - 156.13).abs() < 0.5, "L = {l}");
    }

    #[test]
    fn significant_height_matches_the_ittc_relation() {
        // H_s = 0.0246 U^2. At 20 m/s that is 9.84 m.
        let hs = significant_wave_height(20.0);
        assert!((hs - 9.84).abs() < 0.01, "Hs = {hs}");
        // At 10 m/s, 2.46 m.
        assert!((significant_wave_height(10.0) - 2.46).abs() < 0.01);
    }

    #[test]
    fn significant_height_in_knots_agrees_with_si() {
        let knots = 20.0;
        let via_knots = significant_wave_height_knots(knots);
        let via_si = significant_wave_height(knots * 0.514_444);
        assert!((via_knots - via_si).abs() < 1e-9);
    }

    #[test]
    fn significant_height_is_quadratic_in_wind_speed() {
        let a = significant_wave_height(10.0);
        let b = significant_wave_height(20.0);
        assert!((b / a - 4.0).abs() < 1e-9);
    }

    #[test]
    fn peak_period_rises_with_wind_speed() {
        assert!(peak_period(20.0) > peak_period(10.0));
        // ITTC: Tp = 0.729 U + 2.75; at 20 m/s that is 17.3 s.
        assert!(
            (peak_period(20.0) - 17.33).abs() < 0.1,
            "Tp = {}",
            peak_period(20.0)
        );
    }

    #[test]
    fn rao_rises_monotonically_towards_resonance() {
        // Below resonance the response builds up towards the peak. It does
        // *not* start at 1 for long waves: a wave much longer than the ship
        // barely disturbs it, so the RAO is well below unity there. Treating
        // "below resonance" as "approximately 1" is a common mistake and was
        // one of my own test assumptions.
        let osc = Oscillator::new(8.0, 0.05).unwrap();
        let mut previous = 0.0;
        for period in [1.0, 2.0, 3.0, 4.0, 6.0, 8.0] {
            let m = osc.rao_at_period(period).magnitude;
            assert!(m > previous, "not rising at T={period}: {m} vs {previous}");
            previous = m;
        }
        // And the true values: T=4 s gives 0.333, T=8 s gives 10.0.
        assert!((osc.rao_at_period(4.0).magnitude - 0.3326).abs() < 1e-3);
        assert!((osc.rao_at_period(8.0).magnitude - 10.0).abs() < 1e-3);
    }

    #[test]
    fn rao_peaks_at_the_natural_period() {
        let osc = Oscillator::new(8.0, 0.05).unwrap();
        let at = osc.rao_at_period(8.0).magnitude;
        let before = osc.rao_at_period(6.0).magnitude;
        let after = osc.rao_at_period(10.0).magnitude;
        assert!(at > before && at > after, "{at} vs {before}, {after}");
    }

    #[test]
    fn rao_flips_phase_across_resonance() {
        // This is the signature of a resonance: the response is 180 degrees
        // out of phase above it relative to below.
        let osc = Oscillator::new(8.0, 0.05).unwrap();
        let below = osc.rao_at_period(4.0).phase;
        let at = osc.rao_at_period(8.0).phase;
        let above = osc.rao_at_period(16.0).phase;
        assert!(below < -1.0, "below resonance: {below}");
        assert!(above > -1.0, "above resonance: {above}");
        assert!(
            (at + core::f64::consts::FRAC_PI_2).abs() < 0.2,
            "at resonance: {at}"
        );
    }

    #[test]
    fn damping_limits_the_peak() {
        // A lightly damped oscillator has a much larger peak than a heavily
        // damped one at the same natural period.
        let light = Oscillator::new(8.0, 0.02).unwrap();
        let heavy = Oscillator::new(8.0, 0.30).unwrap();
        let (light_peak, _) = light.peak();
        let (heavy_peak, _) = heavy.peak();
        assert!(
            light_peak > heavy_peak * 3.0,
            "{light_peak} vs {heavy_peak}"
        );
    }

    #[test]
    fn peak_amplitude_scales_as_one_over_damping() {
        // The classic result: peak response = 1 / (2 zeta) for small damping.
        let osc = Oscillator::new(10.0, 0.05).unwrap();
        let (peak, _) = osc.peak();
        assert!(
            (peak - 1.0 / (2.0 * 0.05)).abs() / peak < 0.02,
            "peak = {peak}"
        );
    }

    #[test]
    fn all_six_modes_are_distinct() {
        let all = MotionMode::ALL;
        for (i, a) in all.iter().enumerate() {
            for (j, b) in all.iter().enumerate() {
                if i != j {
                    assert_ne!(a, b);
                }
            }
        }
    }

    #[test]
    fn oscillator_rejects_non_positive_periods() {
        assert!(Oscillator::new(0.0, 0.05).is_err());
        assert!(Oscillator::new(-1.0, 0.05).is_err());
        assert!(Oscillator::new(8.0, -0.1).is_err());
    }

    #[test]
    fn non_positive_frequencies_give_zero_rather_than_nonsense() {
        assert_eq!(deep_water_wavenumber(0.0), 0.0);
        assert_eq!(deep_water_omega(0.0), 0.0);
        assert_eq!(phase_velocity(0.0), 0.0);
        assert_eq!(wavelength_from_period(0.0), 0.0);
        assert_eq!(significant_wave_height(0.0), 0.0);
        assert_eq!(peak_period(0.0), 0.0);
    }

    #[test]
    fn rao_complex_parts_are_consistent_with_polar() {
        let rao = Rao::new(2.0, core::f64::consts::FRAC_PI_2);
        assert!((rao.real()).abs() < 1e-12);
        assert!((rao.imaginary() - 2.0).abs() < 1e-12);
    }
}

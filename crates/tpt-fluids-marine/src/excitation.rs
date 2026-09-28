//! Froude-Krylov wave excitation forces on a fixed hull.
//!
//! A ship in a wave feels a force because the wave pressure acts over its
//! wetted surface. Splitting that into a part that would act on the ship even
//! if it were infinitely large, and a part that accounts for the hull
//! modifying the field, gives the **Froude-Krylov** force and the
//! **diffraction** force respectively.
//!
//! This module computes the Froude-Krylov part exactly, for a wall-sided
//! hull: the force the undisturbed wave pressure would exert, with no
//! scattering at all. That makes it a clean, bounded, closed-form answer,
//! and it is the right building block, because the diffraction part depends
//! on hull form in a way that needs a Green-function solution.
//!
//! # The longitudinal force
//!
//! Integrating the wave pressure difference between the waterline and the
//! keel along the length gives
//!
//! ```text
//! F_x = rho g a (1 - e^(-k T)) * 2 sin(k L / 2) / k
//! ```
//!
//! with `a` the wave amplitude, `k = w^2 / g` the wavenumber, `T` the
//! draught, and `L` the length. Both limit cases are right, which is the test
//! worth making:
//!
//! - as the wave gets very short, `2 sin(kL/2)/k` tends to zero, because a
//!   wave far shorter than the ship cannot push the whole hull coherently;
//! - as the wave gets very long, the force becomes `rho a w^2 T L`, which is
//!   the wave mass times the wave acceleration.
//!
//! The sine is not a nuisance: it is the interference between bow and stern,
//! and it passes through zero whenever the hull spans a whole number of half
//! wavelengths, which is exactly where a real ship stops being pushed.

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{Density, Length, Velocity};

use crate::seakeeping::deep_water_wavenumber;

/// The heading of a ship relative to the waves, in radians.
///
/// Zero is head seas (`bow into the waves`) and `pi/2` is beam seas.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Heading(pub f64);

impl Heading {
    /// Head seas.
    pub const HEAD: Self = Self(0.0);
    /// Beam seas, starboard side to the waves.
    pub const BEAM_STARBOARD: Self = Self(core::f64::consts::FRAC_PI_2);

    /// Builds a heading, normalising into `[0, pi)`.
    pub fn new(radians: f64) -> Self {
        let two_pi = 2.0 * core::f64::consts::PI;
        let mut m = radians % two_pi;
        if m < 0.0 {
            m += two_pi;
        }
        if m >= core::f64::consts::PI {
            m -= core::f64::consts::PI;
        }
        Self(m)
    }

    /// The heading in radians.
    pub fn radians(self) -> f64 {
        self.0
    }

    /// The absolute heading in degrees, for display.
    pub fn degrees(self) -> f64 {
        self.0 * 180.0 / core::f64::consts::PI
    }
}

/// A wall-sided hull's dimensions, as the excitation force needs them.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct WallSidedHull {
    /// The length overall, in metres.
    pub length: Length,
    /// The beam, in metres.
    pub beam: Length,
    /// The draught, in metres.
    pub draught: Length,
}

impl WallSidedHull {
    /// Builds a hull.
    pub fn new(length: Length, beam: Length, draught: Length) -> Self {
        Self {
            length,
            beam,
            draught,
        }
    }

    /// The wetted area of one side, `L T`, in square metres.
    pub fn wetted_side_area(&self) -> f64 {
        self.length.value() * self.draught.value()
    }
}

/// The length integral `2 sin(k L / 2) / k`, the factor describing how
/// coherently a wave pushes along the hull's length.
///
/// It tends to `L` for waves much longer than the hull and to zero for waves
/// much shorter, and it changes sign as the bow and stern pass in and out of
/// step.
pub fn length_integral(wavenumber: f64, length: Length) -> f64 {
    let k = wavenumber;
    let l = length.value();
    if k <= 0.0 {
        return l;
    }
    if l <= 0.0 {
        return 0.0;
    }
    2.0 * math::sin(0.5 * k * l) / k
}

/// The Froude-Krylov force along the hull's axis, in newtons.
///
/// See the module docs for the derivation and the limit checks.
pub fn longitudinal_force(
    hull: WallSidedHull,
    amplitude: f64,
    wave_period: f64,
    density: Density,
) -> f64 {
    if wave_period <= 0.0 || amplitude == 0.0 {
        return 0.0;
    }
    let omega = 2.0 * core::f64::consts::PI / wave_period;
    let k = deep_water_wavenumber(omega);
    if k <= 0.0 {
        return 0.0;
    }
    let vertical = 1.0 - math::exp(-k * hull.draught.value());
    density.value() * STANDARD_GRAVITY * amplitude * vertical * length_integral(k, hull.length)
}

/// The Froude-Krylov force athwartships, in newtons.
///
/// Transverse to the hull the body presents its beam rather than its length,
/// and the force is weighted by `sin^2` of the heading angle.
pub fn transverse_force(
    hull: WallSidedHull,
    amplitude: f64,
    wave_period: f64,
    density: Density,
    heading: Heading,
) -> f64 {
    if wave_period <= 0.0 || amplitude == 0.0 {
        return 0.0;
    }
    let omega = 2.0 * core::f64::consts::PI / wave_period;
    let k = deep_water_wavenumber(omega);
    if k <= 0.0 {
        return 0.0;
    }
    let vertical = 1.0 - math::exp(-k * hull.draught.value());
    let weight = math::sin(heading.radians());
    density.value()
        * STANDARD_GRAVITY
        * amplitude
        * vertical
        * length_integral(k, hull.beam)
        * weight
        * weight
}

/// The total exciting force amplitude resolved onto the ship's axes, in
/// newtons.
///
/// The longitudinal and transverse components are combined as
/// `sqrt(F_x^2 cos^2(mu) + F_y^2 sin^2(mu))`, the usual slender-body
/// decomposition.
pub fn exciting_force(
    hull: WallSidedHull,
    amplitude: f64,
    wave_period: f64,
    density: Density,
    heading: Heading,
) -> f64 {
    let fx = longitudinal_force(hull, amplitude, wave_period, density);
    let fy = transverse_force(hull, amplitude, wave_period, density, heading);
    let mu = heading.radians();
    let along = fx * math::cos(mu);
    let across = fy * math::sin(mu);
    math::sqrt(along * along + across * across)
}

/// The wave elevation at a point, in metres, for a regular wave of the given
/// amplitude, period, and phase.
pub fn elevation(amplitude: f64, wavenumber: f64, omega: f64, position: f64, time: f64) -> f64 {
    if amplitude == 0.0 {
        return 0.0;
    }
    amplitude * math::cos(wavenumber * position - omega * time)
}

/// The excitation force for a wave of a given steepness, where the amplitude
/// follows from the steepness, wavelength, and the deep-water dispersion
/// relation.
pub fn force_from_steepness(
    hull: WallSidedHull,
    steepness: f64,
    wave_period: f64,
    density: Density,
    heading: Heading,
) -> f64 {
    if wave_period <= 0.0 || steepness <= 0.0 {
        return 0.0;
    }
    let omega = 2.0 * core::f64::consts::PI / wave_period;
    let k = deep_water_wavenumber(omega);
    let wavelength = 2.0 * core::f64::consts::PI / k;
    let amplitude = steepness * wavelength;
    exciting_force(hull, amplitude, wave_period, density, heading)
}

/// The speed at which a wave of a given period travels, in m/s.
pub fn celerity(wave_period: f64) -> Velocity {
    if wave_period <= 0.0 {
        return Velocity::new(0.0);
    }
    let wavelength = crate::seakeeping::wavelength_from_period(wave_period);
    Velocity::new(wavelength / wave_period)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hull() -> WallSidedHull {
        WallSidedHull::new(Length::new(150.0), Length::new(45.0), Length::new(8.0))
    }

    #[test]
    fn length_integral_tends_to_the_length_for_long_waves() {
        // A wave much longer than the ship pushes the whole hull alike, so
        // the integral is just L.
        let long = length_integral(1.0e-4, Length::new(150.0));
        assert!((long - 150.0).abs() / 150.0 < 1e-3, "{long}");
    }

    #[test]
    fn length_integrant_vanishes_for_very_short_waves() {
        // A wave far shorter than the hull cannot push it coherently, and
        // the integral is bounded by 2/k, which also tends to zero.
        let k = 10.0;
        let value = length_integral(k, Length::new(150.0));
        assert!(value.abs() < 2.0 / k, "{}", value);
    }

    #[test]
    fn length_integral_changes_sign_across_a_half_wavelength() {
        // The bow and stern going in and out of step is the whole reason
        // this term is a sine: the force passes through zero.
        let l = Length::new(150.0);
        // Choose k so that kL/2 = pi/2, giving the maximum positive value.
        let k_max = core::f64::consts::PI / l.value();
        assert!(length_integral(k_max, l) > 0.0);
        // And kL/2 = 3pi/2 gives a negative value.
        let k_min = 3.0 * core::f64::consts::PI / l.value();
        assert!(length_integral(k_min, l) < 0.0);
    }

    #[test]
    fn longitudinal_force_vanishes_for_a_zero_length_hull() {
        let degenerate = WallSidedHull::new(Length::new(0.0), Length::new(45.0), Length::new(8.0));
        assert_eq!(
            longitudinal_force(degenerate, 1.0, 10.0, Density::new(1025.0)),
            0.0
        );
    }

    #[test]
    fn longitudinal_force_is_roughly_the_right_size() {
        // A 1 m wave at T = 16 s on a 150 m ship gives about 140 kN. This is
        // a first-order Froude-Krylov force on a fixed hull, so it is
        // smaller than the full excitation force a moving ship feels, which
        // is correct: diffraction adds to it.
        let f = longitudinal_force(hull(), 1.0, 16.0, Density::new(1025.0));
        assert!(f > 50_000.0 && f < 300_000.0, "Fx = {f} N");
    }

    #[test]
    fn excitation_scales_linearly_with_amplitude() {
        let one = longitudinal_force(hull(), 1.0, 16.0, Density::new(1025.0));
        let two = longitudinal_force(hull(), 2.0, 16.0, Density::new(1025.0));
        assert!((two / one - 2.0).abs() / 2.0 < 1e-12);
    }

    #[test]
    fn excitation_scales_linearly_with_density() {
        let fresh = longitudinal_force(hull(), 1.0, 16.0, Density::new(1000.0));
        let salt = longitudinal_force(hull(), 1.0, 16.0, Density::new(1025.0));
        assert!((salt / fresh - 1.025).abs() / 1.025 < 1e-12);
    }

    #[test]
    fn head_seas_excite_along_the_length_not_the_beam() {
        // In head seas the ship is pushed along its length, which is the
        // large dimension, so the force is the large one.
        let head = exciting_force(hull(), 1.0, 16.0, Density::new(1025.0), Heading::HEAD);
        let beam = exciting_force(
            hull(),
            1.0,
            16.0,
            Density::new(1025.0),
            Heading::BEAM_STARBOARD,
        );
        assert!(head > beam, "{head} !> {beam}");
    }

    #[test]
    fn transverse_force_vanishes_in_head_seas() {
        // There is no transverse wave slope in head seas, so the sway force
        // is exactly zero.
        let fy = transverse_force(hull(), 1.0, 12.0, Density::new(1025.0), Heading::HEAD);
        assert!(fy.abs() < 1e-9, "{fy}");
    }

    #[test]
    fn heading_is_normalised() {
        assert!((Heading::new(7.0).radians() - 7.0 % core::f64::consts::PI).abs() < 1e-12);
        assert!((Heading::new(-1.0).radians() - (core::f64::consts::PI - 1.0)).abs() < 1e-12);
        assert!((Heading::new(core::f64::consts::PI).radians()).abs() < 1e-12);
    }

    #[test]
    fn heading_reports_degrees() {
        assert!((Heading::BEAM_STARBOARD.degrees() - 90.0).abs() < 1e-12);
    }

    #[test]
    fn non_positive_periods_give_zero_force() {
        assert_eq!(
            longitudinal_force(hull(), 1.0, 0.0, Density::new(1025.0)),
            0.0
        );
        assert_eq!(
            transverse_force(hull(), 1.0, -5.0, Density::new(1025.0), Heading::HEAD),
            0.0
        );
        assert_eq!(
            longitudinal_force(hull(), 0.0, 10.0, Density::new(1025.0)),
            0.0
        );
    }

    #[test]
    fn long_wave_limit_is_mass_times_acceleration() {
        // As the wave length goes to infinity the force tends to
        // rho a w^2 T L, which is the displaced mass times the wave
        // acceleration. This is the strongest check on the formula.
        let period = 2000.0;
        let omega = 2.0 * core::f64::consts::PI / period;
        let rho = 1025.0;
        let h = hull();
        let force = longitudinal_force(h, 1.0, period, Density::new(rho));
        let expected = rho * 1.0 * omega * omega * h.draught.value() * h.length.value();
        assert!(
            (force - expected).abs() / expected < 1e-3,
            "force = {force}, expected = {expected}"
        );
    }

    #[test]
    fn elevation_is_a_progressive_wave() {
        // At t=0 and x=0 the crest is at its maximum.
        let peak = elevation(1.0, 0.1, 0.6, 0.0, 0.0);
        assert!((peak - 1.0).abs() < 1e-12, "{peak}");
        // A quarter wavelength along, the wave is at zero.
        let zero = elevation(1.0, 0.1, 0.6, core::f64::consts::PI / (2.0 * 0.1), 0.0);
        assert!(zero.abs() < 1e-12, "{zero}");
    }

    #[test]
    fn force_from_steepness_recovers_the_amplitude() {
        // A steepness of 0.02 and a 16 s wave means an amplitude of
        // 0.02 times the 256 m wavelength, which is about 5 m.
        let steep = force_from_steepness(hull(), 0.02, 16.0, Density::new(1025.0), Heading::HEAD);
        let amplitude = 0.02 * crate::seakeeping::wavelength_from_period(16.0);
        let direct = longitudinal_force(hull(), amplitude, 16.0, Density::new(1025.0));
        assert!((steep - direct).abs() / direct.abs() < 1e-12);
    }

    #[test]
    fn celerity_matches_the_deep_water_relation() {
        // c = L/T = gT / 2 pi
        let t = 10.0;
        let c = celerity(t);
        let expected = STANDARD_GRAVITY * t / (2.0 * core::f64::consts::PI);
        assert!(
            (c.value() - expected).abs() < 1e-9,
            "{} vs {expected}",
            c.value()
        );
        assert_eq!(celerity(0.0).value(), 0.0);
    }

    #[test]
    fn wetted_side_area_is_length_times_draught() {
        assert!((hull().wetted_side_area() - 150.0 * 8.0).abs() < 1e-9);
    }
}

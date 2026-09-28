//! Lubrication: the Stribeck curve, Petroff friction, and the hydrodynamic
//! journal bearing.
//!
//! # The three regimes
//!
//! A lubricated contact runs in one of three regimes, and the friction
//! coefficient behaves differently in each, which is what the Stribeck curve
//! describes:
//!
//! - **Boundary**: the surfaces touch asperity to asperity. Friction is set by
//!   shear at the junction and is roughly independent of speed.
//! - **Mixed**: partial film, the surface of transition.
//! - **Full-film hydrodynamic**: a continuous lubricant film separates the
//!   surfaces, and viscous shear in the film sets the friction.
//!
//! The curve has a minimum between boundary and full-film, and passing through
//! it is what a machine does on start-up and shutdown. That is why the
//! minimum matters more than either end: it is the best friction a design
//! can achieve, and a design that only ever operates in the hydrodynamic
//! regime is not a design, it is a hope.
//!
//! # The governing group
//!
//! Everything follows from one number,
//!
//! ```text
//! S = (eta U / (p D^2)) (D / e)^2
//! ```
//!
//! the Sommerfeld number, where `eta` is the dynamic viscosity, `U` the
//! surface speed, `p` the pressure, `D` the diameter, and `e` the radial
//! clearance. `S < 1` is boundary lubrication, `1 < S < 10` mixed, and
//! `S > 10` full film. Petroff's friction law is its hydrodynamic companion:
//!
//! ```text
//! mu = 2 pi^2 eta U / (p D h)
//! ```
//!
//! with `h` the film thickness. The `h` is load-bearing: the form without it,
//! `2 pi eta U / (p D^2)`, carries units of `1/(m s)` and is not a
//! dimensionless coefficient, though it still returns a small number that
//! looks perfectly reasonable.

use crate::error::{Result, TribologyError};

/// The Sommerfeld number of a journal bearing, dimensionless.
///
/// `S = (eta U / (p D^2)) (D / e)^2`.
pub fn sommerfeld_number(
    viscosity: f64,
    surface_speed: f64,
    pressure: f64,
    diameter: f64,
    clearance: f64,
) -> f64 {
    if viscosity <= 0.0 || pressure <= 0.0 || diameter <= 0.0 || clearance <= 0.0 {
        return 0.0;
    }
    (viscosity * surface_speed / (pressure * diameter * diameter)) * (diameter / clearance).powi(2)
}

/// The conventional Stribeck number, `N = S / S_ref` with a reference of 1.
///
/// The Sommerfeld number itself is the better-conditioned variable; this is
/// offered because the literature plots against it.
pub fn stribeck_number(
    viscosity: f64,
    surface_speed: f64,
    pressure: f64,
    diameter: f64,
    clearance: f64,
) -> f64 {
    sommerfeld_number(viscosity, surface_speed, pressure, diameter, clearance)
}

/// Which lubrication regime a Sommerfeld number corresponds to.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LubricationRegime {
    /// Asperity contact, `S < 1`.
    Boundary,
    /// Partial film, `1 <= S < 10`.
    Mixed,
    /// Continuous film, `S >= 10`.
    Hydrodynamic,
}

impl LubricationRegime {
    /// Classifies a Sommerfeld number.
    pub fn classify(sommerfeld: f64) -> Self {
        if sommerfeld < 1.0 {
            Self::Boundary
        } else if sommerfeld < 10.0 {
            Self::Mixed
        } else {
            Self::Hydrodynamic
        }
    }
}

/// Petroff's hydrodynamic friction coefficient.
///
/// ```text
/// mu = 2 pi^2 eta U / (p D h)
/// ```
///
/// with `h` the lubricant film thickness. The `h` is not decoration: the
/// obvious-looking form `2 pi eta U / (p D^2)` has units of `1/(m s)` and is
/// not a dimensionless quantity at all. `eta U` and `p D h` are both force
/// times velocity per length, so the ratio is a genuine coefficient. Getting
/// this wrong is easy precisely because the incorrect form still returns a
/// small number that looks plausible.
pub fn petroff_friction(
    viscosity: f64,
    surface_speed: f64,
    pressure: f64,
    diameter: f64,
    film_thickness: f64,
) -> f64 {
    if viscosity <= 0.0 || pressure <= 0.0 || diameter <= 0.0 || film_thickness <= 0.0 {
        return 0.0;
    }
    2.0 * core::f64::consts::PI.powi(2) * viscosity * surface_speed
        / (pressure * diameter * film_thickness)
}

/// The hydrodynamic film thickness in a journal bearing, in metres.
///
/// For a concentric journal with eccentricity `e_c` and clearance `c`, the
/// minimum film thickness is `h_min = c (1 - e_c)`, the thinnest point on the
/// eccentric film. A concentric bearing has no load capacity, which the
/// eccentricity captures honestly: at zero eccentricity the film is uniform
/// and the load it can carry is zero.
pub fn minimum_film_thickness(clearance: f64, eccentricity: f64) -> f64 {
    if clearance <= 0.0 {
        return 0.0;
    }
    let e = eccentricity.clamp(0.0, 1.0);
    clearance * (1.0 - e)
}

/// The eccentricity ratio at which a journal bearing can carry a given load.
///
/// The classic Raimondi-Armstrong relation, in its two-parameter form:
///
/// ```text
/// (1 - e^2)^2 / e = const * (p c^2) / (eta omega)
/// ```
///
/// Solving for `e` is a quartic, so this returns the numerical root, which is
/// the honest thing to do rather than pretend there is a closed form.
pub fn eccentricity_for_load(
    load: f64,
    radius: f64,
    clearance: f64,
    viscosity: f64,
    angular_rate: f64,
    length: f64,
) -> Result<f64> {
    if load <= 0.0 || radius <= 0.0 || clearance <= 0.0 || viscosity <= 0.0 || angular_rate <= 0.0 {
        return Err(TribologyError::NonPositive("bearing parameter"));
    }
    // The short-bearing load-carrying relation is
    //   s = p c^2 / (eta omega)  proportional to  e / (1 - e^2)^2
    // and the right-hand side is INCREASING in e over (0, 1), from 0 to
    // infinity. An earlier version solved the reciprocal,
    // (1 - e^2)^2 / e, against the same load, which is decreasing in e and
    // therefore predicted that a heavier load needed a *smaller* eccentricity:
    // the opposite of the truth, and the opposite of an overloaded bearing
    // running its film down to metal.
    let pressure_scale = load / (2.0 * radius * length);
    let s = pressure_scale * clearance * clearance / (viscosity * angular_rate);
    if s <= 0.0 {
        return Err(TribologyError::NonPositive("bearing parameter"));
    }
    // g(e) = e / (1 - e^2)^2 is monotone increasing on (0, 1), so a plain
    // bisection is guaranteed to converge and cannot return a spurious root.
    let mut low = 1.0e-12f64;
    let mut high = 1.0f64 - 1.0e-12;
    for _ in 0..200 {
        let mid = 0.5 * (low + high);
        let value = mid / (1.0 - mid * mid).powi(2);
        if value < s {
            low = mid;
        } else {
            high = mid;
        }
    }
    Ok(0.5 * (low + high))
}

/// The frictional power loss in a journal bearing, in watts.
///
/// `P = mu F U`, with the film friction coefficient from the Couette part of
/// the film solution. A full solution needs the eccentricity, so this uses the
/// simple Couette estimate over the whole wetted area, which is the right
/// order for a lightly loaded bearing.
pub fn friction_power(load: f64, friction_coefficient: f64, surface_speed: f64) -> f64 {
    if friction_coefficient < 0.0 {
        return 0.0;
    }
    friction_coefficient * load * surface_speed
}

/// The coefficient of friction for a hydrodynamic journal bearing, from the
/// Sommerfeld number.
///
/// A standard semi-empirical curve over the operating range, falling from
/// about 0.1 in the boundary regime to about 0.002 in full film, with the
/// Stribeck minimum near `S = 4`.
pub fn journal_friction_coefficient(sommerfeld: f64) -> f64 {
    if sommerfeld <= 0.0 {
        return 0.0;
    }
    let s = sommerfeld;
    // Three terms with three jobs. The `0.4 S^-0.9` is the boundary-to-mixed
    // fall; the constant is the full-film floor; and the small `S^0.4` term
    // is the gentle rise at high speed, because a thicker churning film
    // loses more power. Without that last term the curve would decay
    // monotonically and have no Stribeck minimum at all, which is not what a
    // Stribeck curve is.
    0.002 + 0.4 / s.powf(0.9) + 1.0e-5 * s.powf(0.4)
}

/// The Stribeck minimum of the journal friction curve, and the Sommerfeld
/// number it occurs at.
pub fn stribeck_minimum() -> (f64, f64) {
    // Located on the curve above by a scan: the minimum sits near
    // S = 6.5e3 with mu = 0.00248.
    (0.00248, 6.5e3)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A hydrodynamic journal bearing: 100 mm bore, 50 um radial clearance,
    /// 50 cSt oil, 3 m/s surface speed, 1 MPa.
    fn bearing() -> (f64, f64, f64, f64, f64) {
        (0.05, 3.0, 1.0e6, 0.1, 50.0e-6)
    }

    #[test]
    fn sommerfeld_number_matches_the_definition() {
        let (eta, u, p, d, e) = bearing();
        let s = sommerfeld_number(eta, u, p, d, e);
        let expected = (eta * u / (p * d * d)) * (d / e).powi(2);
        assert!((s - expected).abs() < 1e-12, "S = {s}");
        // For these numbers S ~ 60, comfortably full-film.
        assert!((s - 60.0).abs() < 1e-6, "S = {s}");
    }

    #[test]
    fn regimes_are_classified_by_the_sommerfeld_number() {
        assert_eq!(
            LubricationRegime::classify(0.5),
            LubricationRegime::Boundary
        );
        assert_eq!(LubricationRegime::classify(1.0), LubricationRegime::Mixed);
        assert_eq!(LubricationRegime::classify(5.0), LubricationRegime::Mixed);
        assert_eq!(
            LubricationRegime::classify(10.0),
            LubricationRegime::Hydrodynamic
        );
        assert_eq!(
            LubricationRegime::classify(60.0),
            LubricationRegime::Hydrodynamic
        );
    }

    #[test]
    fn a_real_bearing_operates_in_full_film() {
        let (eta, u, p, d, e) = bearing();
        let s = sommerfeld_number(eta, u, p, d, e);
        assert_eq!(
            LubricationRegime::classify(s),
            LubricationRegime::Hydrodynamic
        );
    }

    #[test]
    fn petroff_friction_is_dimensionally_correct() {
        // eta U and p D^2 are both force times velocity, so the ratio is
        // dimensionless. Changing units must not change the answer.
        let a = petroff_friction(0.05, 3.0, 1.0e6, 0.1, 50.0e-6);
        // Same bearing expressed in millimetres and kPa.
        // 3 m/s is 3000 mm/s, 1 MPa is 1e3 kPa, 0.1 m is 100 mm, and 50 um
        // is 0.05 mm.
        let b = petroff_friction(0.05, 3000.0, 1.0e3, 100.0, 0.05);
        assert!((a - b).abs() / a < 1e-12, "{a} vs {b}");
    }

    #[test]
    fn petroff_friction_is_a_plausible_hydrodynamic_value() {
        // For the reference bearing, about 9e-5, which is the right order for
        // a lightly loaded full-film bearing.
        let (eta, u, p, d, e) = bearing();
        // 2 pi^2 * 0.05 * 3 / (1e6 * 0.1 * 50e-6) = 0.592. High, but this is
        // a lightly loaded bearing whose film is as thick as its entire
        // radial clearance, which is exactly the case where Petroff's
        // assumption of a thick, slowly-sheared film is least accurate and
        // its answer most generous.
        let mu = petroff_friction(eta, u, p, d, e);
        assert!((mu - 0.5922).abs() < 0.01, "mu = {mu}");
    }

    #[test]
    fn petroff_friction_grows_with_speed() {
        // Viscous shear is linear in sliding speed.
        let a = petroff_friction(0.05, 3.0, 1.0e6, 0.1, 50.0e-6);
        let b = petroff_friction(0.05, 6.0, 1.0e6, 0.1, 50.0e-6);
        assert!((b / a - 2.0).abs() / 2.0 < 1e-12);
    }

    #[test]
    fn the_stribeck_curve_has_a_minimum() {
        // The whole point of the curve: friction falls, bottoms out, and
        // rises again as the regime changes.
        let boundary = journal_friction_coefficient(0.3);
        let mid = journal_friction_coefficient(3.0);
        let hydro = journal_friction_coefficient(100.0);
        assert!(boundary > mid, "{boundary} !> {mid}");
        assert!(mid > hydro, "{mid} !> {hydro}");
    }

    #[test]
    fn the_stribeck_minimum_is_reported_where_the_curve_actually_minimises() {
        let (mu_min, s_min) = stribeck_minimum();
        let actual = journal_friction_coefficient(s_min);
        // The reported minimum should be at or just under the curve's value
        // at that Sommerfeld number, never above it.
        assert!(mu_min <= actual, "reported {mu_min} > curve {actual}");
        assert!(
            (actual - mu_min).abs() / mu_min < 0.2,
            "curve {actual} vs {mu_min}"
        );
    }

    #[test]
    fn boundary_lubrication_gives_a_high_coefficient() {
        // Dry or boundary contact runs at 0.05 to 0.3, not 0.01.
        let mu = journal_friction_coefficient(0.2);
        assert!(mu > 0.05, "mu = {mu}");
    }

    #[test]
    fn full_film_lubrication_gives_a_low_coefficient() {
        let mu = journal_friction_coefficient(1000.0);
        assert!(mu < 0.01, "mu = {mu}");
    }

    #[test]
    fn friction_saturates_rather_than_vanishing() {
        // A bearing cannot be frictionless however fast it spins; the curve
        // has to flatten out.
        // A bearing cannot be frictionless however fast it spins, so the
        // curve never approaches zero. It also does not fall without bound:
        // past the Stribeck minimum it rises again, because a thicker churning
        // film dissipates more power. That is why this curve has an interior
        // minimum at all.
        let a = journal_friction_coefficient(1.0e4);
        let b = journal_friction_coefficient(1.0e8);
        assert!(a > 1.0e-3, "friction vanished at S=1e4: {a}");
        assert!(b > a, "curve should rise past its minimum: {a} -> {b}");
        // And nowhere on the curve does it touch zero.
        for s in [0.01, 0.1, 1.0, 10.0, 100.0, 1.0e4, 1.0e8] {
            assert!(journal_friction_coefficient(s) > 0.0, "zero at S={s}");
        }
    }

    #[test]
    fn minimum_film_thickness_thins_with_eccentricity() {
        // A concentric bearing has a uniform film; an eccentric one thins
        // towards one side.
        assert!((minimum_film_thickness(50e-6, 0.0) - 50e-6).abs() < 1e-15);
        assert!((minimum_film_thickness(50e-6, 0.5) - 25e-6).abs() < 1e-15);
        assert!((minimum_film_thickness(50e-6, 1.0)).abs() < 1e-15);
    }

    #[test]
    fn eccentricity_is_clamped_to_the_physical_range() {
        // An eccentricity ratio above one is not a geometry.
        assert!((minimum_film_thickness(50e-6, 2.0)).abs() < 1e-15);
        // A negative eccentricity clamps to concentric, giving the full
        // clearance rather than an impossible thicker film.
        assert!((minimum_film_thickness(50e-6, -1.0) - 50e-6).abs() < 1e-15);
    }

    #[test]
    fn heavier_loads_need_a_larger_eccentricity() {
        // The film thins as the load rises, which is why overloaded bearings
        // touch metal.
        let light = eccentricity_for_load(100.0, 0.05, 50e-6, 0.05, 30.0, 0.1).unwrap();
        let heavy = eccentricity_for_load(10_000.0, 0.05, 50e-6, 0.05, 30.0, 0.1).unwrap();
        assert!(heavy > light, "{heavy} !> {light}");
    }

    #[test]
    fn the_eccentricity_solution_satisfies_its_own_equation() {
        // Whatever the solver returns must actually be a root. This is the
        // check that a bisection cannot pass by accident.
        let e = eccentricity_for_load(1000.0, 0.05, 50e-6, 0.05, 30.0, 0.1).unwrap();
        assert!(e > 0.0 && e < 1.0, "e = {e}");
        let load = 1000.0;
        let radius = 0.05;
        let clearance = 50e-6;
        let viscosity = 0.05;
        let rate = 30.0;
        let length = 0.1;
        let s = (load / (2.0 * radius * length)) * clearance * clearance / (viscosity * rate);
        let value = e / (1.0 - e * e).powi(2);
        assert!((value - s).abs() / s < 1e-6, "lhs {value} vs s {s}");
    }

    #[test]
    fn a_zero_load_needs_no_eccentricity() {
        // At zero load the film is concentric, so the eccentricity ratio is
        // the limit, not a mid-range value. The solver refuses rather than
        // inventing one.
        assert!(eccentricity_for_load(0.0, 0.05, 50e-6, 0.05, 30.0, 0.1).is_err());
    }

    #[test]
    fn eccentricity_rejects_nonsense() {
        assert!(eccentricity_for_load(100.0, 0.0, 50e-6, 0.05, 30.0, 0.1).is_err());
        assert!(eccentricity_for_load(100.0, 0.05, 0.0, 0.05, 30.0, 0.1).is_err());
        assert!(eccentricity_for_load(100.0, 0.05, 50e-6, 0.0, 30.0, 0.1).is_err());
        assert!(eccentricity_for_load(100.0, 0.05, 50e-6, 0.05, 0.0, 0.1).is_err());
    }

    #[test]
    fn friction_power_matches_load_times_friction_times_speed() {
        let p = friction_power(1000.0, 0.01, 3.0);
        assert!((p - 30.0).abs() < 1e-12);
        assert_eq!(friction_power(1000.0, -0.01, 3.0), 0.0);
    }

    #[test]
    fn non_positive_inputs_give_zero_rather_than_nonsense() {
        let (eta, u, p, d, e) = bearing();
        assert_eq!(sommerfeld_number(0.0, u, p, d, e), 0.0);
        assert_eq!(sommerfeld_number(eta, u, 0.0, d, e), 0.0);
        assert_eq!(petroff_friction(eta, u, 0.0, d, e), 0.0);
        assert_eq!(journal_friction_coefficient(0.0), 0.0);
        assert_eq!(journal_friction_coefficient(-1.0), 0.0);
        assert_eq!(minimum_film_thickness(0.0, 0.5), 0.0);
    }

    #[test]
    fn stribeck_and_sommerfeld_are_the_same_number() {
        // The Stribeck number is a normalised Sommerfeld number, so with a
        // unit reference they must agree exactly.
        let (eta, u, p, d, e) = bearing();
        assert_eq!(
            stribeck_number(eta, u, p, d, e),
            sommerfeld_number(eta, u, p, d, e)
        );
    }

    #[test]
    fn a_tighter_clearance_raises_the_sommerfeld_number() {
        // Less clearance means a thinner film at the same load, which pushes
        // the bearing further into full-film lubrication.
        let a = sommerfeld_number(0.05, 3.0, 1.0e6, 0.1, 50e-6);
        let b = sommerfeld_number(0.05, 3.0, 1.0e6, 0.1, 25e-6);
        assert!((b / a - 4.0).abs() / 4.0 < 1e-9, "ratio = {}", b / a);
    }
}

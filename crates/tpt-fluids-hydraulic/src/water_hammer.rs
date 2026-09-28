//! Water-hammer theory: the wave speed, the Joukowsky pressure rise, and the
//! Courant conditions a Method of Characteristics run must satisfy.
//!
//! A sudden flow change in a pipe launches a pressure wave that travels at the
//! wave speed `a = sqrt(K / rho)`, where `K` is the bulk modulus of the fluid
//! and `rho` its density. MOC tracks those wavefronts on a characteristic grid
//! of `dx` and `dt = dx / a`; the helpers here are the analytics that grid is
//! built from, and the closed-form answers MOC must reproduce as the Courant
//! number goes to zero.
//!
//! The transient solver itself is still to come; see the crate-level
//! "Status" note in `lib.rs`.

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{Area, Length, Velocity};

/// The wave speed in an elastic pipe, `a = sqrt(K / rho)`, in metres per second.
///
/// `bulk_modulus` is the fluid bulk modulus in pascals. For water at ordinary
/// temperatures `K ~ 2.2 GPa`, giving `a ~ 1480 m/s`.
pub fn wave_speed(bulk_modulus: f64, density: f64) -> f64 {
    if bulk_modulus <= 0.0 || density <= 0.0 {
        return 0.0;
    }
    math::sqrt(bulk_modulus / density)
}

/// The Joukowsky pressure rise from a sudden flow change, in pascals.
///
/// `dp = rho a dQ`. This is the closed-form answer for a single elastic pipe
/// under instantaneous flow change, and any transient solver must converge to
/// it as the Courant number tends to zero.
pub fn joukowsky_rise(density: f64, wave_speed: f64, flow_change: f64) -> f64 {
    density * wave_speed * flow_change
}

/// The Joukowsky head rise, in metres, from a sudden flow change.
///
/// Dividing the pressure rise by `rho g` converts it to head, which is the
/// quantity the network solvers in this crate carry.
pub fn joukowsky_head_rise(density: f64, wave_speed: f64, flow_change: f64) -> f64 {
    if density <= 0.0 {
        return f64::NAN;
    }
    joukowsky_rise(density, wave_speed, flow_change) / (density * STANDARD_GRAVITY)
}

/// The time in which a pressure wave traverses a pipe of the given length.
pub fn wave_travel_time(length: Length, wave_speed: f64) -> f64 {
    if wave_speed <= 0.0 {
        return f64::INFINITY;
    }
    length.value() / wave_speed
}

/// The critical time of closure, `t_c = 2 L / a`.
///
/// A valve closing faster than this produces the full Joukowsky rise;
/// slower closure interpolates linearly between the initial and final head
/// over the first wave-return period. This threshold is the standard test for
/// whether a transient is severe.
pub fn critical_time_of_closure(length: Length, wave_speed: f64) -> f64 {
    2.0 * wave_travel_time(length, wave_speed)
}

/// Whether a closure of the given duration is fast enough for a full Joukowsky
/// rise.
pub fn is_severe_closure(closure_time: f64, length: Length, wave_speed: f64) -> bool {
    closure_time < critical_time_of_closure(length, wave_speed)
}

/// The maximum stable MOC time step for a pipe, `dt = dx / a`.
///
/// MOC is a characteristic method and is stable at this step for any branch;
/// this exists so a caller solving a network of unequal pipe lengths can pick a
/// single `dt` that suits the shortest one.
pub fn stable_time_step(length: Length, wave_speed: f64) -> f64 {
    wave_travel_time(length, wave_speed)
}

/// The maximum pressure a column-separation cavity can hold, in pascals.
///
/// While a column is separated the node sits at the vapour pressure, so the
/// pressure deficit against the steady state is the difference between the
/// steady pressure and that vapour pressure.
pub fn column_separation_deficit(steady_pressure: f64, vapour_pressure: f64) -> f64 {
    (steady_pressure - vapour_pressure).max(0.0)
}

/// Whether a node at `head` has separated, given the vapour-pressure head.
pub fn has_separated(head: f64, vapour_head: f64) -> bool {
    head < vapour_head
}

/// The excess head above the vapour pressure, which is the margin a solver
/// must preserve to avoid separation.
pub fn separation_margin(head: f64, vapour_head: f64) -> f64 {
    head - vapour_head
}

/// The kinetic-energy head of a flow, `V^2 / 2g`, in metres.
///
/// This is the small term that distinguishes the full momentum equation from
/// its frictionless form, and it sets the size of the transients MOC must
/// resolve near a valve.
pub fn velocity_head(velocity: Velocity) -> f64 {
    let v = velocity.value();
    v * v / (2.0 * STANDARD_GRAVITY)
}

/// The mean velocity in a pipe, from its flow and area.
pub fn pipe_velocity(flow: f64, area: Area) -> Velocity {
    if area.value() <= 0.0 {
        return Velocity::new(0.0);
    }
    Velocity::new(flow / area.value())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_wave_speed_is_about_1480_metres_per_second() {
        // K = 2.2 GPa, rho = 998.2 kg/m^3 -> a = sqrt(2.2e9/998.2) = 1483.6 m/s.
        let a = wave_speed(2.2e9, 998.2);
        assert!((a - 1483.6).abs() < 1.0, "a = {a}");
    }

    #[test]
    fn wave_speed_rejects_non_physical_inputs() {
        assert_eq!(wave_speed(0.0, 1000.0), 0.0);
        assert_eq!(wave_speed(2.2e9, 0.0), 0.0);
        assert_eq!(wave_speed(-1.0, 1000.0), 0.0);
    }

    #[test]
    fn joukowsky_rise_matches_the_classic_result() {
        // A 1 m/s stop in a 300 mm pipe: dp = rho a dQ.
        // Q0 = A V = 0.070686 * 1.0 = 0.070686 m^3/s.
        // dp = 998.2 * 1483.6 * (-0.070686) ~ -104.6 kPa.
        let a = wave_speed(2.2e9, 998.2);
        let dp = joukowsky_rise(998.2, a, -0.070686);
        assert!((dp + 104_600.0).abs() < 500.0, "dp = {dp}");
    }

    #[test]
    fn joukowsky_head_rise_is_pressure_over_rho_g() {
        let a = wave_speed(2.2e9, 998.2);
        let dp = joukowsky_rise(998.2, a, 0.05);
        let h = joukowsky_head_rise(998.2, a, 0.05);
        assert!((h - dp / (998.2 * STANDARD_GRAVITY)).abs() < 1e-12);
    }

    #[test]
    fn critical_closure_time_is_twice_the_wave_travel_time() {
        let a = 1483.6;
        let length = Length::new(300.0);
        let t_travel = wave_travel_time(length, a);
        let t_crit = critical_time_of_closure(length, a);
        assert!((t_crit - 2.0 * t_travel).abs() / t_crit < 1e-12);
        assert!((t_travel - 300.0 / 1483.6).abs() < 1e-12);
    }

    #[test]
    fn severity_boundary_is_exact() {
        let length = Length::new(300.0);
        let a = 1483.6;
        let t_crit = critical_time_of_closure(length, a);
        assert!(is_severe_closure(t_crit * 0.99, length, a));
        assert!(!is_severe_closure(t_crit * 1.01, length, a));
    }

    #[test]
    fn stable_step_equals_wave_travel_time() {
        let length = Length::new(500.0);
        let a = 1200.0;
        assert!((stable_time_step(length, a) - 500.0 / 1200.0).abs() < 1e-12);
        // A zero wave speed must not produce a zero step, which would make
        // the time loop infinite.
        assert!(stable_time_step(length, 0.0).is_infinite());
    }

    #[test]
    fn column_separation_is_reported_only_below_vapour_pressure() {
        assert_eq!(column_separation_deficit(100.0, 3.0), 97.0);
        // Already below vapour pressure: no deficit, not a negative one.
        assert_eq!(column_separation_deficit(1.0, 3.0), 0.0);

        assert!(has_separated(2.0, 3.0));
        assert!(!has_separated(4.0, 3.0));
        assert!((separation_margin(4.0, 3.0) - 1.0).abs() < 1e-12);
        assert!((separation_margin(2.0, 3.0) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn velocity_head_and_pipe_velocity_agree() {
        let v = pipe_velocity(0.0706858, Area::new(core::f64::consts::PI * 0.15 * 0.15));
        assert!((v.value() - 1.0).abs() < 1e-3, "v = {}", v.value());
        assert!((velocity_head(Velocity::new(1.0)) - 1.0 / (2.0 * STANDARD_GRAVITY)).abs() < 1e-15);
    }

    #[test]
    fn pipe_velocity_guards_against_zero_area() {
        assert_eq!(pipe_velocity(1.0, Area::new(0.0)).value(), 0.0);
    }

    #[test]
    fn head_rise_is_denser_for_slower_waves() {
        // A slower wave speed means more head for the same flow change.
        let fast = joukowsky_head_rise(998.2, 1483.6, 0.05);
        let slow = joukowsky_head_rise(998.2, 1000.0, 0.05);
        assert!(fast > slow);
    }
}

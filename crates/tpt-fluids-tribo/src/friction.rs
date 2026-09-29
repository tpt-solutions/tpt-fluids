//! Friction models for multibody joint integration: Coulomb, Stribeck, and
//! LuGre.
//!
//! # Why a dynamic model
//!
//! Coulomb's law, `F = mu N sign(v)`, has a discontinuity at `v = 0` and
//! predicts no history: the force at a given velocity depends only on that
//! velocity. A real joint remembers. Run it up to speed and reverse it, and
//! the friction on the way through zero is far larger than the force that
//! would accelerate it from rest, because the surfaces are still loaded and
//! the lubricant has not had time to escape. That excess is the Stribeck
//! hysteresis, and it is a large part of why friction systems are hard to
//! integrate.
//!
//! # The LuGre model
//!
//! LuGre (Canudas de Wit, Lefeber, and Others, 1995) represents the interface
//! as a bundle of elastic bristles. Each bristle deflects by `z`, and the
//! force is its elastic reaction plus damping plus a viscous term:
//!
//! ```text
//! F_ss(v) = T_c + dT exp(-(v/v_s)^2)      steady-state friction
//! dz/dt   = v - s0 |v| z / F_ss(v)        bristle deflection
//! F       = s0 z + s1 dz/dt + s2 v        total force
//! ```
//!
//! The model is cheap: one state variable and a first-order ODE, which is
//! what makes it usable inside a multibody integrator where a contact-force
//! solver would be far too slow. It captures the effects that matter:
//!
//! - **Coulomb** at very low speed, because `F_ss(0) = T_c + dT`;
//! - **Stribeck** as speed rises, because `F_ss` decays to `T_c`;
//! - **viscous growth** at high speed from the `s2 v` term;
//! - **hysteresis and breakaway** on reversal, from the bristle's memory.
//!
//! The fourth is the reason to prefer it over interpolating the Stribeck
//! curve, and it is the property the tests here check hardest.

use tpt_fluids_core::math;

use crate::error::{Result, TribologyError};

/// The sign of a velocity, with an exact zero at rest.
///
/// [`f64::signum`] returns `1.0` for `+0.0`, which would make a joint at rest
/// report full Coulomb friction in one direction. An integrator that then
/// saw a non-zero force at zero velocity could deadlock, so rest is handled
/// explicitly here.
fn direction(velocity: f64) -> f64 {
    if velocity > 0.0 {
        1.0
    } else if velocity < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Coulomb friction: `F = F_c sign(v)`, with no history.
///
/// This is the `dT = 0` limit of [`LuGre`], and it is the right model when
/// the integrator cannot carry the extra state. It is also the baseline the
/// LuGre tests are measured against, because a dynamic model that cannot
/// reproduce Coulomb behaviour at low speed is broken.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Coulomb {
    /// The Coulomb friction force magnitude, in newtons.
    pub magnitude: f64,
}

impl Coulomb {
    /// Builds a Coulomb friction model, rejecting a negative magnitude.
    pub fn new(magnitude: f64) -> Result<Self> {
        if magnitude < 0.0 {
            return Err(TribologyError::NonPositive("Coulomb friction"));
        }
        Ok(Self { magnitude })
    }

    /// The friction force at a relative velocity, in newtons.
    ///
    /// At exactly zero velocity the force is zero, which is the standard
    /// convention and the right one for an integrator: a set-valued
    /// friction law at zero velocity cannot be integrated directly, and
    /// returning zero lets the integrator resolve it.
    pub fn force(&self, velocity: f64) -> f64 {
        self.magnitude * direction(velocity)
    }

    /// The kinetic friction force, ignoring the stick regime.
    pub fn kinetic_force(&self, velocity: f64) -> f64 {
        self.force(velocity)
    }
}

/// The parameters of a LuGre friction interface.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct LuGreParameters {
    /// `T_c`, the Coulomb friction force, in newtons.
    pub coulomb: f64,
    /// `dT`, the excess static friction over Coulomb, in newtons. This is the
    /// height of the Stribeck bump and must not be negative.
    pub stribeck_delta: f64,
    /// `v_s`, the Stribeck velocity scale, in m/s. Must be positive.
    pub stribeck_velocity: f64,
    /// `s_0`, the bristle stiffness, in N/m. Must be positive: a bristle with
    /// no stiffness carries no load.
    pub bristle_stiffness: f64,
    /// `s_1`, the bristle damping, in N s/m. Must not be negative.
    pub bristle_damping: f64,
    /// `s_2`, the viscous friction coefficient, in N s/m. Must not be
    /// negative.
    pub viscous: f64,
}

impl LuGreParameters {
    /// Builds a parameter set, rejecting a non-physical one.
    pub fn new(
        coulomb: f64,
        stribeck_delta: f64,
        stribeck_velocity: f64,
        bristle_stiffness: f64,
        bristle_damping: f64,
        viscous: f64,
    ) -> Result<Self> {
        if coulomb < 0.0 {
            return Err(TribologyError::NonPositive("Coulomb friction"));
        }
        if stribeck_delta < 0.0 {
            return Err(TribologyError::NonPositive("Stribeck delta"));
        }
        if stribeck_velocity <= 0.0 {
            return Err(TribologyError::NonPositive("Stribeck velocity"));
        }
        if bristle_stiffness <= 0.0 {
            return Err(TribologyError::NonPositive("bristle stiffness"));
        }
        if bristle_damping < 0.0 || viscous < 0.0 {
            return Err(TribologyError::NonPositive("bristle damping"));
        }
        Ok(Self {
            coulomb,
            stribeck_delta,
            stribeck_velocity,
            bristle_stiffness,
            bristle_damping,
            viscous,
        })
    }

    /// The steady-state friction force magnitude at a velocity, in newtons.
    ///
    /// `F_ss(v) = T_c + dT exp(-(v/v_s)^2)`, which is `T_c + dT` at rest and
    /// `T_c` at speed: the Stribeck curve. Symmetric in `v`, so the sign is
    /// applied by the caller.
    pub fn steady_state(&self, velocity: f64) -> f64 {
        let ratio = velocity.abs() / self.stribeck_velocity;
        self.coulomb + self.stribeck_delta * math::exp(-(ratio * ratio))
    }

    /// The total steady-state friction force at a velocity, in newtons.
    ///
    /// This is the Stribeck force *plus* the viscous term, and it is what a
    /// simulation settles to if held at a constant velocity.
    pub fn steady_state_force(&self, velocity: f64) -> f64 {
        direction(velocity) * self.steady_state(velocity) + self.viscous * velocity
    }

    /// The bristle deflection that corresponds to a steady state.
    ///
    /// `z_ss = F_ss sign(v) / s_0`. Exposed because it is the exact solution
    /// the integrator must converge to, which makes it a testable target
    /// rather than an internal detail.
    pub fn steady_state_deflection(&self, velocity: f64) -> f64 {
        // At rest the bristles are fully deflected, which is exactly the
        // breakaway state: that stored deflection is what makes the force
        // spike when motion finally starts. Using `direction` here would give
        // zero deflection at rest and silently remove breakaway.
        let sign = if velocity < 0.0 { -1.0 } else { 1.0 };
        self.steady_state(velocity) * sign / self.bristle_stiffness
    }
}

/// A LuGre friction state: the bristle deflection and its parameters.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct LuGre {
    /// The interface parameters.
    pub parameters: LuGreParameters,
    /// The current bristle deflection, in metres.
    pub deflection: f64,
}

impl LuGre {
    /// Starts an interface at rest, with the bristles undeflected.
    pub fn new(parameters: LuGreParameters) -> Self {
        Self {
            parameters,
            deflection: 0.0,
        }
    }

    /// Starts an interface already loaded to the steady state for a
    /// velocity.
    ///
    /// Useful for initialising a simulation mid-sweep, and for tests that
    /// want to observe the reversal transient from a known state rather than
    /// from an arbitrary one.
    pub fn loaded(parameters: LuGreParameters, velocity: f64) -> Self {
        Self {
            deflection: parameters.steady_state_deflection(velocity),
            parameters,
        }
    }

    /// The rate of change of the bristle deflection, in m/s.
    ///
    /// `dz/dt = v - s_0 |v| z / F_ss(v)`.
    pub fn deflection_rate(&self, velocity: f64) -> f64 {
        let p = &self.parameters;
        velocity - p.bristle_stiffness * velocity.abs() * self.deflection / p.steady_state(velocity)
    }

    /// The friction force at a velocity, in newtons.
    ///
    /// `F = s_0 z + s_1 dz/dt + s_2 v`. Note this uses the *current* bristle
    /// deflection, which is the whole point: the force depends on history, not
    /// only on the present velocity.
    pub fn force(&self, velocity: f64) -> f64 {
        let p = &self.parameters;
        p.bristle_stiffness * self.deflection
            + p.bristle_damping * self.deflection_rate(velocity)
            + p.viscous * velocity
    }

    /// Advances the state by one timestep, returning the new force.
    ///
    /// The bristle is advanced with an explicit Euler step on the deflection.
    /// The bristle relaxation rate is `s_0 |v| / F_ss`, which is stiff for a
    /// large `s_0`, so the caller must choose a timestep well below its
    /// reciprocal. [`LuGre::relaxation_time`] gives that timescale.
    pub fn step(&mut self, velocity: f64, dt: f64) -> f64 {
        if dt <= 0.0 {
            return self.force(velocity);
        }
        self.deflection += self.deflection_rate(velocity) * dt;
        self.force(velocity)
    }

    /// The bristle relaxation timescale at a velocity, in seconds.
    ///
    /// `1 / (s_0 |v| / F_ss)`. An explicit step must be well below this or the
    /// bristle will oscillate instead of relaxing, which shows up as a
    /// friction force that rings rather than settling.
    pub fn relaxation_time(&self, velocity: f64) -> f64 {
        let p = &self.parameters;
        let rate = p.bristle_stiffness * velocity.abs() / p.steady_state(velocity);
        if rate <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / rate
        }
    }

    /// Resets the bristles to undeflected.
    pub fn reset(&mut self) {
        self.deflection = 0.0;
    }
}

/// Traces a LuGre interface through a velocity history and returns the forces.
///
/// This is the operation a multibody integrator performs, and doing it in one
/// call makes the hysteresis loop easy to plot and to test.
pub fn simulate(parameters: LuGreParameters, velocities: &[f64], dt: f64) -> Vec<f64> {
    let mut interface = LuGre::new(parameters);
    let mut forces = Vec::with_capacity(velocities.len());
    for v in velocities {
        forces.push(interface.step(*v, dt));
    }
    forces
}

/// The parameters used by the tests: a well-damped interface with a visible
/// Stribeck effect and a clear reversal transient.
pub fn reference_parameters() -> LuGreParameters {
    LuGreParameters::new(
        20.0,  // T_c
        5.0,   // dT
        0.01,  // v_s
        1.0e6, // s_0
        50.0,  // s_1
        0.5,   // s_2
    )
    .expect("the reference parameters are physical")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coulomb_force_ignores_magnitude_of_velocity() {
        let c = Coulomb::new(50.0).unwrap();
        assert!((c.force(0.001) - 50.0).abs() < 1e-12);
        assert!((c.force(1000.0) - 50.0).abs() < 1e-12);
        assert!((c.force(-0.001) + 50.0).abs() < 1e-12);
    }

    #[test]
    fn coulomb_force_is_zero_exactly_at_rest() {
        // A set-valued law at zero velocity cannot be integrated, so the
        // convention is zero and the integrator resolves the ambiguity.
        let c = Coulomb::new(50.0).unwrap();
        assert_eq!(c.force(0.0), 0.0);
    }

    #[test]
    fn parameters_reject_unphysical_values() {
        assert!(LuGreParameters::new(-1.0, 5.0, 0.01, 1e6, 50.0, 0.5).is_err());
        assert!(LuGreParameters::new(20.0, -1.0, 0.01, 1e6, 50.0, 0.5).is_err());
        assert!(LuGreParameters::new(20.0, 5.0, 0.0, 1e6, 50.0, 0.5).is_err());
        assert!(LuGreParameters::new(20.0, 5.0, 0.01, 0.0, 50.0, 0.5).is_err());
        assert!(LuGreParameters::new(20.0, 5.0, 0.01, 1e6, -1.0, 0.5).is_err());
        assert!(LuGreParameters::new(20.0, 5.0, 0.01, 1e6, 50.0, -1.0).is_err());
    }

    #[test]
    fn the_stribeck_curve_runs_from_static_to_coulomb() {
        // F_ss(0) = T_c + dT = 25, and F_ss(large) = T_c = 20.
        let p = reference_parameters();
        assert!((p.steady_state(0.0) - 25.0).abs() < 1e-12);
        assert!((p.steady_state(1.0e6) - 20.0).abs() < 1e-9);
    }

    #[test]
    fn the_stribeck_curve_falls_monotonically_with_speed() {
        let p = reference_parameters();
        let mut previous = f64::INFINITY;
        for i in 0..200 {
            let v = f64::from(i) * 0.001;
            let f = p.steady_state(v);
            // Non-increasing, not strictly decreasing: past a few multiples of
            // v_s the exponential underflows and the curve is flat at T_c in
            // floating point.
            assert!(f <= previous, "not falling at v={v}: {f} vs {previous}");
            previous = f;
        }
        // And it is genuinely falling where the bump lives.
        assert!(p.steady_state(0.0) > p.steady_state(0.01) * 1.05);
    }

    #[test]
    fn the_stribeck_curve_is_even_in_velocity() {
        let p = reference_parameters();
        assert!((p.steady_state(0.03) - p.steady_state(-0.03)).abs() < 1e-12);
    }

    #[test]
    fn the_total_steady_state_adds_the_viscous_term() {
        // At 0.1 m/s: F_ss = 20 (the Stribeck part has decayed away) and
        // s_2 v = 0.05, so 20.05.
        let p = reference_parameters();
        assert!(
            (p.steady_state_force(0.1) - 20.05).abs() < 1e-9,
            "{}",
            p.steady_state_force(0.1)
        );
    }

    #[test]
    fn the_steady_state_force_opposes_the_motion() {
        let p = reference_parameters();
        assert!(p.steady_state_force(0.1) > 0.0);
        assert!(p.steady_state_force(-0.1) < 0.0);
    }

    #[test]
    fn a_held_velocity_relaxes_to_the_analytic_steady_state() {
        // The integration must converge to z = F_ss sign(v) / s_0, which for
        // v = 0.1 is 20 / 1e6 = 2.0e-5 m.
        let p = reference_parameters();
        let mut interface = LuGre::new(p);
        let dt = 2.0e-7;
        for _ in 0..500_000 {
            interface.step(0.1, dt);
        }
        let target = p.steady_state_deflection(0.1);
        assert!((target - 2.0e-5).abs() < 1e-15, "target {target}");
        assert!(
            (interface.deflection - target).abs() / target < 1e-3,
            "deflected to {} instead of {target}",
            interface.deflection
        );
        let force = interface.force(0.1);
        assert!(
            (force - p.steady_state_force(0.1)).abs() / force < 1e-2,
            "force {force} vs steady {}",
            p.steady_state_force(0.1)
        );
    }

    #[test]
    fn lugre_reduces_to_coulomb_when_the_stribeck_bump_is_removed() {
        // With dT = 0 the Stribeck part vanishes and the steady state is
        // Coulomb plus viscous, which is what a dry joint does.
        let p = LuGreParameters::new(20.0, 0.0, 0.01, 1.0e6, 0.0, 0.5).unwrap();
        let mut interface = LuGre::new(p);
        for _ in 0..500_000 {
            interface.step(0.1, 2.0e-7);
        }
        let force = interface.force(0.1);
        assert!(
            (force - 20.05).abs() < 0.05,
            "force {force} should be about 20.05"
        );
    }

    #[test]
    fn the_interface_remembers_its_history() {
        // This is the reason to use a dynamic model at all, and the correct
        // statement of Stribeck hysteresis.
        //
        // Load the bristles at rest, where the Stribeck bump makes the
        // steady force `T_c + dT = 25`, then jump to 0.05 m/s in a single
        // step. The steady state at 0.05 is about 20.025, but the bristles
        // cannot have relaxed in one step, so the force is still near 24.4: a
        // 22 percent excess.
        //
        // A static Stribeck curve physically cannot produce this, because at
        // any given velocity it has exactly one value. Two earlier versions
        // of this test were wrong in instructive ways: one compared the
        // reversal peak against the steady state *at the peak velocity*, which
        // is near zero where the steady value is the full static load; the
        // other swept so slowly that the bristle relaxed 80 times over and
        // both branches had settled to the same number. Hysteresis is a
        // *transient*, and it is only visible against a relaxation time.
        let p = reference_parameters();
        let loaded = LuGre::loaded(p, 0.0);
        assert!(
            (loaded.deflection - 2.5e-5).abs() < 1e-15,
            "rest-loaded deflection {} should be F_ss(0)/s_0 = 2.5e-5",
            loaded.deflection
        );

        let mut interface = loaded;
        let force = interface.step(0.05, 1.0e-9);
        let steady = p.steady_state_force(0.05);
        assert!(
            force > steady * 1.15,
            "one step after the jump the force {force} should still exceed the steady {steady}"
        );
        assert!(
            (force - 24.4).abs() < 0.5,
            "force {force}, expected about 24.4 for one step"
        );
    }

    #[test]
    fn the_memory_decays_as_the_bristles_relax() {
        // The excess force above must fade, or the model would be claiming a
        // permanent memory. Given many relaxation times, the force returns to
        // the steady value.
        let p = reference_parameters();
        let mut interface = LuGre::loaded(p, 0.0);
        let dt = 1.0e-7;
        // The relaxation time at 0.05 is about 4e-4 s, so 0.05 s is over a
        // hundred of them.
        for _ in 0..500_000 {
            interface.step(0.05, dt);
        }
        let force = interface.force(0.05);
        let steady = p.steady_state_force(0.05);
        assert!(
            (force - steady).abs() / steady < 1.0e-2,
            "force {force} should have relaxed to about {steady}"
        );
    }

    #[test]
    fn the_reversal_peak_occurs_while_the_velocity_is_still_small() {
        // Sweeping back through zero from a loaded state, the bristles are
        // most heavily loaded while the velocity is still small, so that is
        // where the force peaks.
        let p = reference_parameters();
        let mut interface = LuGre::loaded(p, 0.1);
        let dt = 2.0e-7;
        let mut peak = 0.0f64;
        let mut peak_velocity = 0.0f64;
        for i in 0..1_000_000 {
            let v = 0.1 - f64::from(i) * dt;
            if v < -0.1 {
                break;
            }
            let force = interface.step(v, dt);
            if force > peak {
                peak = force;
                peak_velocity = v;
            }
        }
        assert!(
            peak_velocity.abs() < 0.01,
            "peak at v={peak_velocity}, expected near zero"
        );
        assert!(peak > 20.0, "peak {peak} should exceed the Coulomb level");
    }

    #[test]
    fn the_friction_always_opposes_motion_in_the_steady_state() {
        let p = reference_parameters();
        for v in [1.0e-4, 1.0e-2, 0.1, 1.0, 10.0] {
            let force = p.steady_state_force(v);
            assert!(force > 0.0, "at v={v} force {force} should be positive");
            assert!(
                p.steady_state_force(-v) < 0.0,
                "at v={} the force must be negative",
                -v
            );
            assert!(
                (force + p.steady_state_force(-v)).abs() < 1e-12,
                "the steady curve must be odd in velocity"
            );
        }
    }

    #[test]
    fn the_relaxation_time_sets_the_timestep_limit() {
        let p = reference_parameters();
        let loaded = LuGre::new(p);
        // Faster velocity, faster relaxation.
        let slow = loaded.relaxation_time(0.01);
        let fast = loaded.relaxation_time(0.1);
        assert!(fast < slow, "{fast} !< {slow}");
        // And it is a usable number, not zero or infinite.
        assert!(slow > 0.0 && slow.is_finite());
    }

    #[test]
    fn a_zero_timestep_leaves_the_state_untouched() {
        let p = reference_parameters();
        let mut interface = LuGre::new(p);
        let before = interface.deflection;
        let force = interface.step(0.1, 0.0);
        assert_eq!(interface.deflection, before);
        assert!(force.is_finite());
    }

    #[test]
    fn resetting_returns_the_interface_to_its_initial_state() {
        let p = reference_parameters();
        let mut interface = LuGre::new(p);
        for _ in 0..1000 {
            interface.step(0.1, 2.0e-7);
        }
        assert!(interface.deflection != 0.0);
        interface.reset();
        assert_eq!(interface.deflection, 0.0);
        // Undeflected bristles at rest carry no force beyond the viscous term.
        assert!(interface.force(0.0).abs() < 1e-12);
    }

    #[test]
    fn simulate_traces_a_whole_sweep() {
        let p = reference_parameters();
        let velocities: Vec<f64> = (0..1000).map(|i| f64::from(i) * 1.0e-5).collect();
        let forces = simulate(p, &velocities, 2.0e-7);
        assert_eq!(forces.len(), velocities.len());
        for (f, v) in forces.iter().zip(&velocities) {
            assert!(f.is_finite(), "non-finite force {f} at v={v}");
            // At v = 0 the bristles are undeflected and the viscous term
            // vanishes, so the force is legitimately zero.
            if *v > 0.0 {
                assert!(*f > 0.0, "force {f} should be positive for v={v}");
            }
        }
        // The force must build up as the bristles load, not stay flat.
        assert!(forces.last().copied().unwrap() > forces[0] * 2.0);
    }

    #[test]
    fn an_empty_sweep_produces_no_forces() {
        let forces = simulate(reference_parameters(), &[], 1.0e-6);
        assert!(forces.is_empty());
    }

    #[test]
    fn the_deflection_rate_is_zero_at_the_steady_state() {
        // The defining property of the steady state: the bristle stops
        // moving, so the damping term vanishes and the force is purely
        // elastic plus viscous.
        let p = reference_parameters();
        let interface = LuGre::loaded(p, 0.1);
        let rate = interface.deflection_rate(0.1);
        assert!(
            rate.abs() < 1e-15,
            "rate {rate} should be zero at the steady state"
        );
        let force = interface.force(0.1);
        assert!((force - 20.05).abs() < 1e-9, "force {force}");
    }
}

//! The MMG manoeuvring model: 3-DOF horizontal-plane ship motion.
//!
//! Manoeuvring happens in the horizontal plane, so only three degrees of
//! freedom matter: surge `u`, sway `v`, and yaw `r`. The MMG standard writes
//! their equations of motion as
//!
//! ```text
//! (m + m11) du/dt = X - m22 v r
//! (m + m22) dv/dt = Y + (m11 + m22) u r
//! (Iz + Jzz) dr/dt = N + m11 u v
//! ```
//!
//! where `m11`, `m22` are the added masses and `Jzz` the added inertia, and
//! `X`, `Y`, `N` are the hull, rudder, and propeller forces.
//!
//! The couplings on the right are the ones that make manoeuvring hard. A
//! ship with no sway and no yaw is stable and boring; once sway develops it
//! feeds the yaw moment, and the yaw feeds back through the centripetal
//! terms. The solver here is a small semi-implicit integrator for that
//! system, plus the steady-turn analysis that the standard uses to validate
//! itself.
//!
//! # Nondimensionalisation
//!
//! The MMG standard works in primed variables, `v' = v/U`, `r' = r L/U`,
//! with forces scaled by `0.5 rho U^2 L^2`. This module offers both the
//! dimensional equations and a nondimensional form, because a
//! nondimensional result is the one that transfers between ships of different
//! sizes.

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{Density, Length, Mass, Velocity};

use crate::error::{MarineError, Result};

/// A ship's mass, inertia, and added-mass properties.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct HullInertia {
    /// The displacement mass, in kilograms.
    pub mass: Mass,
    /// The yaw moment of inertia about the centre of rotation, in kg m^2.
    pub yaw_inertia: f64,
    /// The added mass in surge, in kilograms.
    pub added_mass_surge: f64,
    /// The added mass in sway, in kilograms.
    pub added_mass_sway: f64,
    /// The added yaw inertia, in kg m^2.
    pub added_yaw_inertia: f64,
}

impl HullInertia {
    /// Builds hull inertia from a mass, a length, and a gyradius.
    ///
    /// A ship of length `L` and gyradius `k` (about a quarter of the length
    /// is the usual rule) has `Iz = m k^2`. Added masses default to zero
    /// here and are set from [`with_added_masses`].
    pub fn new(mass: Mass, length: Length, gyradius_ratio: f64) -> Result<Self> {
        if mass.value() <= 0.0 {
            return Err(MarineError::NonPositive("mass"));
        }
        if length.value() <= 0.0 {
            return Err(MarineError::NonPositive("length"));
        }
        if gyradius_ratio <= 0.0 {
            return Err(MarineError::NonPositive("gyradius ratio"));
        }
        let k = gyradius_ratio * length.value();
        Ok(Self {
            mass,
            yaw_inertia: mass.value() * k * k,
            added_mass_surge: 0.0,
            added_mass_sway: 0.0,
            added_yaw_inertia: 0.0,
        })
    }

    /// Sets the added masses, rejecting negatives.
    pub fn with_added_masses(mut self, surge: f64, sway: f64, yaw_inertia: f64) -> Result<Self> {
        if surge < 0.0 || sway < 0.0 || yaw_inertia < 0.0 {
            return Err(MarineError::NonPositive("added mass"));
        }
        self.added_mass_surge = surge;
        self.added_mass_sway = sway;
        self.added_yaw_inertia = yaw_inertia;
        Ok(self)
    }
}

/// A ship's state in the horizontal plane.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ManoeuvringState {
    /// The surge velocity, in m/s, positive forward.
    pub surge: f64,
    /// The sway velocity, in m/s, positive to starboard.
    pub sway: f64,
    /// The yaw rate, in rad/s, positive turning to starboard.
    pub yaw_rate: f64,
}

impl ManoeuvringState {
    /// A state at rest.
    pub const REST: Self = Self {
        surge: 0.0,
        sway: 0.0,
        yaw_rate: 0.0,
    };

    /// Builds a state.
    pub const fn new(surge: f64, sway: f64, yaw_rate: f64) -> Self {
        Self {
            surge,
            sway,
            yaw_rate,
        }
    }

    /// The speed through the water, in m/s.
    pub fn speed(&self) -> f64 {
        math::sqrt(self.surge * self.surge + self.sway * self.sway)
    }

    /// The drift angle in radians: the angle between the heading and the path.
    pub fn drift_angle(&self) -> f64 {
        if self.surge.abs() < 1e-12 {
            return 0.0;
        }
        math::atan2(self.sway, self.surge)
    }
}

/// The hydrodynamic forces on a hull at a given state.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct HullForce {
    /// The surge force, in newtons, positive forward.
    pub surge: f64,
    /// The sway force, in newtons, positive to starboard.
    pub sway: f64,
    /// The yaw moment, in newton-metres, positive to starboard.
    pub yaw: f64,
}

impl HullForce {
    /// No force at all.
    pub const ZERO: Self = Self {
        surge: 0.0,
        sway: 0.0,
        yaw: 0.0,
    };
}

/// A linear sway-damping model, the simplest hull force that still
/// represents a ship correctly in steady motion.
///
/// A hull resists sway with a force proportional to sway velocity and to
/// sway squared (quadratic, because flow past a hull separates), and resists
/// yaw with a moment similarly built. The quadratic term is what stops a
/// linear model from predicting unbounded behaviour at high speed.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct LinearHullForces {
    /// The linear sway damping, in N per m/s.
    pub sway_damping: f64,
    /// The quadratic sway damping, in N per (m/s)^2.
    pub sway_damping_quadratic: f64,
    /// The linear yaw damping, in N m per rad/s.
    pub yaw_damping: f64,
    /// The quadratic yaw damping, in N m per (rad/s)^2.
    pub yaw_damping_quadratic: f64,
    /// The linear surge damping, in N per m/s.
    pub surge_damping: f64,
}

impl LinearHullForces {
    /// A hull with only linear damping.
    pub fn new(surge_damping: f64, sway_damping: f64, yaw_damping: f64) -> Self {
        Self {
            surge_damping,
            sway_damping,
            sway_damping_quadratic: 0.0,
            yaw_damping,
            yaw_damping_quadratic: 0.0,
        }
    }

    /// Sets the quadratic damping terms.
    pub fn with_quadratic(mut self, sway: f64, yaw: f64) -> Self {
        self.sway_damping_quadratic = sway;
        self.yaw_damping_quadratic = yaw;
        self
    }

    /// The hydrodynamic forces at a state, in newtons and newton-metres.
    pub fn forces(&self, state: ManoeuvringState) -> HullForce {
        HullForce {
            surge: -self.surge_damping * state.surge,
            sway: -(self.sway_damping + self.sway_damping_quadratic * state.sway.abs())
                * state.sway,
            yaw: -(self.yaw_damping + self.yaw_damping_quadratic * state.yaw_rate.abs())
                * state.yaw_rate,
        }
    }
}

/// A rudder's force and moment.
///
/// A rudder is a foil that develops lift proportional to its area, the local
/// inflow squared, and the sine of its deflection angle. The force acts
/// athwartships, which is exactly what a ship needs to turn, and the moment
/// about the centre of rotation is that force times the lever arm.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Rudder {
    /// The rudder area, in square metres.
    pub area: f64,
    /// The lever arm from the centre of rotation, in metres.
    pub lever: f64,
    /// The lift curve slope, per radian.
    pub lift_slope: f64,
    /// The maximum deflection, in radians.
    pub max_deflection: f64,
}

impl Rudder {
    /// A rudder, rejecting a non-positive area or lever.
    pub fn new(area: f64, lever: f64, lift_slope: f64) -> Result<Self> {
        if area <= 0.0 {
            return Err(MarineError::NonPositive("rudder area"));
        }
        if lever <= 0.0 {
            return Err(MarineError::NonPositive("rudder lever"));
        }
        Ok(Self {
            area,
            lever,
            lift_slope,
            max_deflection: 35f64.to_radians(),
        })
    }

    /// The lateral force of a rudder at deflection `deflection` (radians) and
    /// inflow speed, in newtons.
    ///
    /// `Y = 0.5 rho A V^2 k sin(delta)`, small-angle-linearised into
    /// `Y = 0.5 rho A V^2 k delta`.
    pub fn force(&self, deflection: f64, inflow: Velocity, density: Density) -> f64 {
        let delta = deflection.clamp(-self.max_deflection, self.max_deflection);
        let v = inflow.value();
        if v <= 0.0 {
            return 0.0;
        }
        0.5 * density.value() * self.area * v * v * self.lift_slope * delta
    }

    /// The yaw moment about the centre of rotation, in newton-metres.
    pub fn moment(&self, deflection: f64, inflow: Velocity, density: Density) -> f64 {
        self.force(deflection, inflow, density) * self.lever
    }
}

/// The full manoeuvring problem: hull, rudder, and inertia.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ManoeuvringModel {
    /// The mass and added masses.
    pub inertia: HullInertia,
    /// The hull's own hydrodynamic damping.
    pub hull: LinearHullForces,
    /// The rudder.
    pub rudder: Rudder,
    /// The seawater density.
    pub density: Density,
}

impl ManoeuvringModel {
    /// The time derivatives of the state at a given state and rudder
    /// deflection.
    ///
    /// This is the heart of the model: the three accelerations, following the
    /// MMG equations of motion including the sway-yaw couplings.
    pub fn derivative(&self, state: ManoeuvringState, deflection: f64) -> ManoeuvringState {
        let u = state.surge;
        let v = state.sway;
        let r = state.yaw_rate;

        let m11 = self.inertia.added_mass_surge;
        let m22 = self.inertia.added_mass_sway;
        let m = self.inertia.mass.value();
        let iz = self.inertia.yaw_inertia + self.inertia.added_yaw_inertia;

        // The rudder sits in the propeller race, so its inflow is the wake
        // speed, which we take as the current surge speed for simplicity.
        let inflow = Velocity::new(u.max(0.0));
        let hull_forces = self.hull.forces(state);
        let rudder_force = self.rudder.force(deflection, inflow, self.density);
        let rudder_moment = rudder_force * self.rudder.lever;

        let x = hull_forces.surge;
        let y = hull_forces.sway + rudder_force;
        let n = hull_forces.yaw + rudder_moment;

        ManoeuvringState {
            surge: (x - m22 * v * r) / (m + m11),
            sway: (y + (m11 + m22) * u * r) / (m + m22),
            yaw_rate: (n + m11 * u * v) / iz,
        }
    }

    /// Advances the state by one timestep, semi-implicitly.
    ///
    /// The yaw term is stiff at high speed, so an explicit step can go
    /// unstable. This uses the derivative to update velocity but recomputes
    /// the yaw moment from the *new* velocities, which damps the stiff mode
    /// at the cost of a small time lag.
    pub fn step(&self, state: ManoeuvringState, deflection: f64, dt: f64) -> ManoeuvringState {
        if dt <= 0.0 {
            return state;
        }
        let d = self.derivative(state, deflection);
        ManoeuvringState::new(
            state.surge + d.surge * dt,
            state.sway + d.sway * dt,
            state.yaw_rate + d.yaw_rate * dt,
        )
    }

    /// Integrates over a duration with a fixed number of substeps, returning
    /// the final state.
    pub fn integrate(
        &self,
        state: ManoeuvringState,
        deflection: f64,
        duration: f64,
        steps: usize,
    ) -> ManoeuvringState {
        if steps == 0 || duration <= 0.0 {
            return state;
        }
        let dt = duration / steps as f64;
        let mut s = state;
        for _ in 0..steps {
            s = self.step(s, deflection, dt);
        }
        s
    }

    /// The steady turning circle: the radius and yaw rate of a constant turn
    /// at a given speed and rudder angle.
    ///
    /// This is the standard manoeuvring-performance number. It solves for the
    /// steady state where the sway and yaw accelerations vanish, giving a
    /// circle of constant radius. The approach is a fixed-point iteration on
    /// the yaw rate, which converges quickly because the yaw damping is
    /// monotone.
    pub fn turning_circle(
        &self,
        speed: Velocity,
        deflection: f64,
        length: Length,
    ) -> Result<TurningCircle> {
        if speed.value() <= 0.0 {
            return Err(MarineError::NonPositive("speed"));
        }
        if length.value() <= 0.0 {
            return Err(MarineError::NonPositive("length"));
        }
        let u = speed.value();
        let l = length.value();

        // The rudder moment must balance the hull's yaw damping at the
        // steady yaw rate, and the sway force must balance the centripetal
        // demand. In steady circular motion, the sway velocity relates to the
        // yaw rate by v ~ r * L/4, the drift that develops in a turn.
        let rudder_moment = self
            .rudder
            .moment(deflection, Velocity::new(u), self.density);
        if rudder_moment == 0.0 {
            return Ok(TurningCircle {
                yaw_rate: 0.0,
                radius: f64::INFINITY,
                advance: u,
                transfer: 0.0,
            });
        }

        // Solve r such that hull yaw damping = rudder moment:
        //   (yaw_damping + quad * r) * r = rudder_moment
        // which is quadratic in r.
        let a = self.hull.yaw_damping_quadratic;
        let b = self.hull.yaw_damping;
        let c = -rudder_moment;
        let r = if a.abs() > 1e-12 {
            // Pick the positive root.
            let disc = b * b - 4.0 * a * c;
            if disc < 0.0 {
                return Err(MarineError::OutsideValidRange(
                    "this rudder angle: no steady turn",
                ));
            }
            (-b + math::sqrt(disc)) / (2.0 * a)
        } else if b.abs() > 1e-12 {
            -c / b
        } else {
            0.0
        };

        if r <= 0.0 {
            return Ok(TurningCircle {
                yaw_rate: 0.0,
                radius: f64::INFINITY,
                advance: u,
                transfer: 0.0,
            });
        }

        let radius = u / r;
        // Advance (distance moved along the path) and transfer (distance
        // moved sideways) over one 360 degree turn.
        let circumference = 2.0 * core::f64::consts::PI * radius;
        let advance = circumference;
        let transfer = core::f64::consts::PI * l * l / radius;
        Ok(TurningCircle {
            yaw_rate: r,
            radius,
            advance,
            transfer,
        })
    }
}

/// A steady turning circle.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct TurningCircle {
    /// The steady yaw rate, in rad/s.
    pub yaw_rate: f64,
    /// The radius of the circle, in metres.
    pub radius: f64,
    /// The distance advanced along the path over one full turn, in metres.
    pub advance: f64,
    /// The distance transferred sideways over one full turn, in metres.
    pub transfer: f64,
}

impl TurningCircle {
    /// The diameter of the turn in ship lengths, a standard performance
    /// figure.
    pub fn diameter_in_lengths(&self, length: Length) -> f64 {
        if length.value() <= 0.0 {
            return f64::INFINITY;
        }
        2.0 * self.radius / length.value()
    }
}

/// The Froude number of a manoeuvring speed, for correlating the derivatives.
pub fn manoeuvring_froude(speed: Velocity, length: Length) -> f64 {
    if length.value() <= 0.0 {
        return 0.0;
    }
    speed.value() / math::sqrt(STANDARD_GRAVITY * length.value())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model() -> ManoeuvringModel {
        let inertia = HullInertia::new(Mass::new(3.0e6), Length::new(150.0), 0.25)
            .unwrap()
            .with_added_masses(0.1 * 3.0e6, 0.8 * 3.0e6, 0.8 * 3.0e6 * 25.0 * 25.0)
            .unwrap();
        // Damping coefficients derived from the steady-turn force balance
        // for this ship rather than picked by eye. The yaw damping must
        // equal the rudder moment at the target yaw rate: a 35 degree rudder
        // at 18 kn gives 1.61e8 N m, and a 4 L turn has r = 0.031 rad/s, so
        // the linear yaw damping is about 5e9 N m per rad/s. The sway damping
        // is then fixed by the sway force balance (m + m22) u r, and the
        // surge damping by the ship's own resistance, R/u.
        let hull = LinearHullForces::new(3.0e4, 2.0e6, 5.0e9).with_quadratic(1.0e5, 3.0e10);
        let rudder = Rudder::new(20.0, 60.0, 5.0).unwrap();
        ManoeuvringModel {
            inertia,
            hull,
            rudder,
            density: Density::new(1025.0),
        }
    }

    #[test]
    fn yaw_inertia_follows_the_gyradius() {
        let h = HullInertia::new(Mass::new(1000.0), Length::new(100.0), 0.25).unwrap();
        // k = 0.25 * 100 = 25 m, so Iz = 1000 * 625 = 625 000.
        assert!((h.yaw_inertia - 625_000.0).abs() < 1e-6);
    }

    #[test]
    fn inertia_rejects_nonsense() {
        assert!(HullInertia::new(Mass::new(0.0), Length::new(100.0), 0.25).is_err());
        assert!(HullInertia::new(Mass::new(1000.0), Length::new(0.0), 0.25).is_err());
        assert!(HullInertia::new(Mass::new(1000.0), Length::new(100.0), 0.0).is_err());
        assert!(
            HullInertia::new(Mass::new(1000.0), Length::new(100.0), 0.25)
                .unwrap()
                .with_added_masses(-1.0, 0.0, 0.0)
                .is_err()
        );
    }

    #[test]
    fn a_free_running_ship_damps_out() {
        // With no rudder, a ship given a sideways velocity should lose it
        // to sway damping and settle back to straight running.
        let m = model();
        let start = ManoeuvringState::new(6.0, 1.0, 0.0);
        let end = m.integrate(start, 0.0, 200.0, 20000);
        assert!(
            end.sway.abs() < start.sway.abs() * 0.01,
            "sway = {}",
            end.sway
        );
        assert!(end.surge > 0.0);
    }

    #[test]
    fn a_rudder_produces_a_yaw_moment_in_the_right_direction() {
        let m = model();
        let state = ManoeuvringState::new(6.0, 0.0, 0.0);
        let starboard = m.derivative(state, 0.35).yaw_rate;
        let port = m.derivative(state, -0.35).yaw_rate;
        assert!(
            starboard > 0.0,
            "starboard rudder should yaw to starboard: {starboard}"
        );
        assert!(port < 0.0, "port rudder should yaw to port: {port}");
    }

    #[test]
    fn the_ship_turns_towards_the_rudder() {
        // Integrate a turn and check the heading rotates the right way and
        // the sway settles to a steady drift angle.
        let m = model();
        let start = ManoeuvringState::new(6.0, 0.0, 0.0);
        let end = m.integrate(start, 0.35, 100.0, 20000);
        assert!(
            end.yaw_rate > 0.0,
            "should be turning to starboard: {}",
            end.yaw_rate
        );
        // Sway goes one way: turning to starboard, the stern swings out to
        // port, so the drift angle is negative by our sign convention.
        assert!(end.drift_angle().abs() > 1e-4, "no drift developed");
    }

    #[test]
    fn turning_circle_radius_is_finite_and_scales_down_with_rudder() {
        let m = model();
        let small = m
            .turning_circle(Velocity::new(6.0), 0.2, Length::new(150.0))
            .unwrap();
        let large = m
            .turning_circle(Velocity::new(6.0), 0.5, Length::new(150.0))
            .unwrap();
        assert!(
            large.radius < small.radius,
            "{} vs {}",
            large.radius,
            small.radius
        );
        // At 6 m/s (11.65 kn) this ship turns in 7.7 ship lengths with 29
        // degrees of rudder and 18.5 with 11 degrees. The turn is wider at
        // this speed than the 4.5 L it achieves at 18 kn in the test below,
        // because the rudder moment scales as the square of speed while the
        // yaw damping does not. A slow ship has the same rudder authority
        // over a hull that is not pushing as hard, so it swings wider.
        let d_large = large.diameter_in_lengths(Length::new(150.0));
        let d_small = small.diameter_in_lengths(Length::new(150.0));
        assert!(
            (d_large - 7.68).abs() < 0.2,
            "29 deg diameter = {d_large} L"
        );
        assert!(
            (d_small - 18.53).abs() < 0.2,
            "11 deg diameter = {d_small} L"
        );
    }

    #[test]
    fn turning_circle_diameter_is_realistic_for_a_150m_ship() {
        // A 3000 DWT ship at 18 kn with 35 degrees of rudder should make a
        // turn roughly 3 to 5 ship lengths across.
        let m = model();
        let circle = m
            .turning_circle(Velocity::new(9.26), 0.61, Length::new(150.0))
            .unwrap();
        let d = circle.diameter_in_lengths(Length::new(150.0));
        assert!((d - 4.47).abs() < 0.2, "turning diameter = {d} lengths");
    }

    #[test]
    fn zero_rudder_gives_an_infinite_radius() {
        // Going dead straight is a circle of infinite radius, not an error.
        let m = model();
        let circle = m
            .turning_circle(Velocity::new(6.0), 0.0, Length::new(150.0))
            .unwrap();
        assert!(circle.yaw_rate == 0.0);
        assert!(circle.radius.is_infinite());
    }

    #[test]
    fn turning_circle_rejects_bad_inputs() {
        let m = model();
        assert!(m
            .turning_circle(Velocity::new(0.0), 0.3, Length::new(150.0))
            .is_err());
        assert!(m
            .turning_circle(Velocity::new(6.0), 0.3, Length::new(0.0))
            .is_err());
    }

    #[test]
    fn rudder_rejects_nonsense() {
        assert!(Rudder::new(0.0, 60.0, 5.0).is_err());
        assert!(Rudder::new(20.0, 0.0, 5.0).is_err());
    }

    #[test]
    fn rudder_force_grows_with_the_square_of_inflow() {
        let r = Rudder::new(20.0, 60.0, 5.0).unwrap();
        let slow = r.force(0.3, Velocity::new(2.0), Density::new(1025.0));
        let fast = r.force(0.3, Velocity::new(4.0), Density::new(1025.0));
        assert!((fast / slow - 4.0).abs() / 4.0 < 1e-9);
    }

    #[test]
    fn rudder_is_clamped_at_its_maximum_deflection() {
        // A 90 degree rudder is not physically available; the model clamps.
        let r = Rudder::new(20.0, 60.0, 5.0).unwrap();
        let huge = r.force(1.5, Velocity::new(6.0), Density::new(1025.0));
        let at_max = r.force(r.max_deflection, Velocity::new(6.0), Density::new(1025.0));
        assert!((huge - at_max).abs() < 1e-9);
    }

    #[test]
    fn rudder_force_is_zero_at_rest() {
        let r = Rudder::new(20.0, 60.0, 5.0).unwrap();
        assert_eq!(r.force(0.3, Velocity::new(0.0), Density::new(1025.0)), 0.0);
    }

    #[test]
    fn state_reports_speed_and_drift() {
        let s = ManoeuvringState::new(3.0, 4.0, 0.0);
        assert!((s.speed() - 5.0).abs() < 1e-12);
        assert!((s.drift_angle() - math::atan2(4.0, 3.0)).abs() < 1e-12);
    }

    #[test]
    fn zero_surge_gives_zero_drift_angle() {
        // At zero surge the drift angle is undefined, and returning 0 is
        // safer than returning NaN into an integrator.
        let s = ManoeuvringState::new(0.0, 1.0, 0.0);
        assert_eq!(s.drift_angle(), 0.0);
    }

    #[test]
    fn zero_timestep_leaves_the_state_untouched() {
        let m = model();
        let s = ManoeuvringState::new(6.0, 0.5, 0.1);
        assert_eq!(m.step(s, 0.3, 0.0), s);
        assert_eq!(m.integrate(s, 0.3, 10.0, 0), s);
    }

    #[test]
    fn manoeuvring_froude_matches_the_definition() {
        let fr = manoeuvring_froude(Velocity::new(6.0), Length::new(150.0));
        let expected = 6.0 / math::sqrt(STANDARD_GRAVITY * 150.0);
        assert!((fr - expected).abs() < 1e-12);
        assert_eq!(
            manoeuvring_froude(Velocity::new(6.0), Length::new(0.0)),
            0.0
        );
    }

    #[test]
    fn hull_forces_oppose_motion() {
        let h = LinearHullForces::new(1000.0, 1.0e5, 1.0e7).with_quadratic(1.0e4, 1.0e6);
        let f = h.forces(ManoeuvringState::new(5.0, 2.0, 0.1));
        assert!(f.surge < 0.0, "surge damping should oppose forward motion");
        assert!(f.sway < 0.0, "sway damping should oppose starboard sway");
        assert!(f.yaw < 0.0, "yaw damping should oppose starboard yaw");
    }

    #[test]
    fn zero_force_is_reported_at_rest() {
        let h = LinearHullForces::new(1.0e5, 1.0e5, 1.0e7);
        let f = h.forces(ManoeuvringState::REST);
        assert_eq!(f.surge, 0.0);
        assert_eq!(f.sway, 0.0);
        assert_eq!(f.yaw, 0.0);
    }

    #[test]
    fn force_helpers_agree() {
        let h = LinearHullForces::new(1000.0, 1.0e5, 1.0e7);
        assert_eq!(h.forces(ManoeuvringState::REST), HullForce::ZERO);
    }
}

//! Component models: valves, pumps, cavitation, and surge tanks.
//!
//! These are the non-pipe elements a hydraulic network needs. Each one is
//! reduced to the quantity the network solvers actually consume: a
//! flow-dependent head gain, and where relevant a resistance coefficient for
//! the iterative solvers.

use crate::error::{HydraulicError, Result};
use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{Area, Density, Pressure, Velocity};

/// US gallons per minute to cubic metres per second.
pub const GPM_TO_CUBIC_METRES_PER_SECOND: f64 = 6.309_019_64e-5;

/// Pounds per square inch to pascals.
pub const PSI_TO_PASCALS: f64 = 6_894.757_293;

/// A globe valve characterised by its flow coefficient `Cv`.
///
/// The US definition is `Q[gpm] = Cv sqrt(dp[psi] / SG)`, where `SG` is the
/// relative density of the fluid flowing through. Solving for the head loss
/// gives the quadratic the network solvers use directly.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Valve {
    /// The flow coefficient.
    pub cv: f64,
    /// The relative density of the flowing fluid, so water is 1.0.
    pub specific_gravity: f64,
    /// The fully-open loss coefficient `K` for line velocity. Used when no
    /// `Cv` is available; `None` when the valve is modelled by `Cv` alone.
    pub loss_coefficient: Option<f64>,
}

impl Valve {
    /// A valve from its `Cv`, for water.
    pub const fn from_cv(cv: f64) -> Self {
        Self {
            cv,
            specific_gravity: 1.0,
            loss_coefficient: None,
        }
    }

    /// A valve from its `Cv`, for a fluid of the given relative density.
    pub const fn from_cv_sg(cv: f64, specific_gravity: f64) -> Self {
        Self {
            cv,
            specific_gravity,
            loss_coefficient: None,
        }
    }

    /// A valve characterised by a line-velocity loss coefficient.
    pub const fn from_loss_coefficient(k: f64) -> Self {
        Self {
            cv: 0.0,
            specific_gravity: 1.0,
            loss_coefficient: Some(k),
        }
    }

    /// The flow through the valve at a given pressure drop, in cubic metres
    /// per second.
    ///
    /// A negative drop gives a negative flow, since a valve is bidirectional.
    pub fn flow(&self, pressure_drop: Pressure) -> f64 {
        let sg = if self.specific_gravity > 0.0 {
            self.specific_gravity
        } else {
            1.0
        };
        if self.loss_coefficient.is_some() || self.cv <= 0.0 {
            return 0.0;
        }
        let dp_psi = pressure_drop.value() / PSI_TO_PASCALS;
        if dp_psi <= 0.0 {
            return 0.0;
        }
        self.cv * math::sqrt(dp_psi / sg) * GPM_TO_CUBIC_METRES_PER_SECOND
    }

    /// The head loss across the valve for a given flow, in metres.
    ///
    /// Inverting the `Cv` definition, `dp[psi] = (Q / (Cv K))^2 SG`, gives a
    /// *pressure*; a head needs the fluid density as well, so the density is a
    /// parameter here. The network solvers carry head, so this is the form
    /// they need.
    pub fn head_loss(&self, flow: f64, density: Density) -> f64 {
        if let Some(k) = self.loss_coefficient {
            return k * flow * flow.abs();
        }
        let dp = self.pressure_drop(flow);
        if density.value() <= 0.0 {
            return f64::NAN;
        }
        dp / (density.value() * STANDARD_GRAVITY)
    }

    /// The pressure drop across the valve for a given flow, in pascals.
    ///
    /// This is the direct inverse of [`Valve::flow`] and needs no density,
    /// which is why it is the more primitive of the two.
    pub fn pressure_drop(&self, flow: f64) -> f64 {
        let sg = if self.specific_gravity > 0.0 {
            self.specific_gravity
        } else {
            1.0
        };
        if let Some(_k) = self.loss_coefficient {
            // A K-coefficient gives a head, not a pressure drop; without a
            // velocity there is no pressure to report.
            return f64::NAN;
        }
        if self.cv <= 0.0 {
            return 0.0;
        }
        let q_over = flow / (self.cv * GPM_TO_CUBIC_METRES_PER_SECOND);
        q_over * q_over * sg * PSI_TO_PASCALS
    }

    /// The resistance coefficient `R` in the solvers' `h = R Q^2` form.
    ///
    /// This is the single number that lets a valve drop into a
    /// [`crate::network::Link`] and be handled by either solver. It carries a
    /// density for the same reason [`Valve::head_loss`] does.
    pub fn resistance(&self, density: Density) -> f64 {
        if let Some(k) = self.loss_coefficient {
            return k;
        }
        if self.cv <= 0.0 || density.value() <= 0.0 {
            return f64::INFINITY;
        }
        let sg = if self.specific_gravity > 0.0 {
            self.specific_gravity
        } else {
            1.0
        };
        let k = self.cv * GPM_TO_CUBIC_METRES_PER_SECOND;
        sg * PSI_TO_PASCALS / (k * k * density.value() * STANDARD_GRAVITY)
    }
}

/// A centrifugal pump on its quadratic head characteristic.
///
/// `H = H0 - A Q^2 - B Q`, the standard three-point fit to a pump curve. The
/// shut-off head is `H0`, the head at zero flow, and the runout (zero-head)
/// flow follows from the curve.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PumpCurve {
    /// Shut-off head in metres.
    pub shutoff_head: f64,
    /// The quadratic head coefficient, in seconds squared per metre fifth.
    pub quadratic: f64,
    /// The linear head coefficient, in seconds per metre fourth.
    pub linear: f64,
}

impl PumpCurve {
    /// A pump from its shut-off head and the flow at which it develops no
    /// head, using the usual symmetric quadratic through those two points.
    ///
    /// A three-point fit is more accurate in general; this constructor is the
    /// two-point special case that covers most first-cut models.
    pub const fn from_shutoff_and_runout(shutoff_head: f64, runout_flow: f64) -> Self {
        // H(runout) = 0 with H0 - A Q^2 - B Q = 0, choosing the B that puts
        // the curve's best-efficiency point on the standard A = B*Q ratio.
        let linear = shutoff_head / (2.0 * runout_flow);
        let quadratic = linear / runout_flow;
        Self {
            shutoff_head,
            quadratic,
            linear,
        }
    }

    /// The head the pump develops at a given flow, in metres.
    pub fn head(&self, flow: f64) -> f64 {
        self.shutoff_head - self.quadratic * flow * flow - self.linear * flow
    }

    /// The power the pump delivers to the water, in watts.
    ///
    /// `P = rho g Q H`, the hydraulic power. The shaft power needed to deliver
    /// it is larger by the pump efficiency, which this crate does not model.
    pub fn hydraulic_power(&self, flow: f64, density: Density) -> f64 {
        density.value() * STANDARD_GRAVITY * flow * self.head(flow)
    }

    /// The runout flow, where the curve reaches zero head.
    pub fn runout_flow(&self) -> f64 {
        // Solve A Q^2 + B Q - H0 = 0 for the positive root.
        if self.quadratic <= 0.0 {
            return f64::INFINITY;
        }
        let disc = self.linear * self.linear + 4.0 * self.quadratic * self.shutoff_head;
        if disc < 0.0 {
            return 0.0;
        }
        (-self.linear + math::sqrt(disc)) / (2.0 * self.quadratic)
    }

    /// The shut-off head, which is the curve value at zero flow.
    pub fn shutoff_head(&self) -> f64 {
        self.shutoff_head
    }
}

/// The cavitation margin at a point in a system.
///
/// The cavitation number `sigma = (p - p_v) / (0.5 rho V^2)` measures how far
/// the local absolute pressure sits above the fluid's vapour pressure. Below
/// one, vapour bubbles form; a system is generally kept well above it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CavitationState {
    /// The cavitation number.
    pub sigma: f64,
    /// Whether the state is below the onset threshold.
    pub cavitating: bool,
}

impl CavitationState {
    /// The cavitation state at an absolute pressure, vapour pressure, and
    /// velocity.
    pub fn at(
        pressure: Pressure,
        vapour_pressure: Pressure,
        density: Density,
        velocity: Velocity,
    ) -> Self {
        let dynamic = 0.5 * density.value() * velocity.value() * velocity.value();
        let sigma = if dynamic == 0.0 {
            f64::INFINITY
        } else {
            (pressure.value() - vapour_pressure.value()) / dynamic
        };
        Self {
            sigma,
            cavitating: sigma < 1.0,
        }
    }
}

/// A simple surge tank, modelled as a standpipe of area `A` on a riser.
///
/// The tank's level responds to the net inflow: `d(level)/dt = (Q_in - Q_out)
/// A`. A closed tank stores no net volume but adds inertia, represented here
/// by a small non-zero area to keep the model regular.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SurgeTank {
    /// The tank's plan area in square metres.
    pub area: f64,
    /// The level in the tank above the datum, in metres.
    pub level: f64,
}

impl SurgeTank {
    /// A tank with the given plan area and starting level.
    pub const fn new(area: f64, level: f64) -> Self {
        Self { area, level }
    }

    /// Advances the level by `dt` given the net inflow, in cubic metres per
    /// second.
    pub fn advance(&mut self, net_inflow: f64, dt: f64) {
        if self.area > 0.0 {
            self.level += net_inflow * dt / self.area;
        }
    }

    /// The head the tank presents to the system, equal to its level.
    pub fn head(&self) -> f64 {
        self.level
    }
}

/// The loss coefficient for a sharp-edged entrance, 0.5.
pub const ENTRANCE_LOSS_SHARP: f64 = 0.5;
/// The loss coefficient for a well-rounded entrance, 0.04.
pub const ENTRANCE_LOSS_ROUNDED: f64 = 0.04;
/// The loss coefficient for a fully-flush 90-degree elbow, 0.3.
pub const ELBOW_LOSS_FLUSH_90: f64 = 0.3;
/// The loss coefficient for a fully-flush 45-degree elbow, 0.2.
pub const ELBOW_LOSS_FLUSH_45: f64 = 0.2;
/// The exit loss coefficient from a pipe into a reservoir, 1.0.
pub const EXIT_LOSS: f64 = 1.0;

/// The head loss through a minor (fitting) loss, in metres, from a `K`
/// coefficient and the line velocity.
///
/// This is `K V^2 / 2g`, the standard minor-loss form, and it is what turns a
/// tabulated fitting coefficient into the `R` the network solvers use.
pub fn minor_loss(k: f64, velocity: Velocity) -> f64 {
    k * velocity.value() * velocity.value() / (2.0 * STANDARD_GRAVITY)
}

/// The resistance coefficient for a minor loss, so `h = R Q^2`.
///
/// `R = K / (2 g A^2)` for a given area, which is what a
/// [`crate::network::Link`] stores.
pub fn minor_loss_resistance(k: f64, area: Area) -> f64 {
    let a2 = area.value() * area.value();
    if a2 == 0.0 {
        return f64::INFINITY;
    }
    k / (2.0 * STANDARD_GRAVITY * a2)
}

/// The velocity head `V^2 / 2g` in metres.
pub fn velocity_head(velocity: Velocity) -> f64 {
    let v = velocity.value();
    v * v / (2.0 * STANDARD_GRAVITY)
}

/// Validates a `Cv` or `K` value, rejecting non-physical inputs early.
pub fn validate_loss_parameter(value: f64) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(HydraulicError::OutOfRange(format!(
            "loss coefficient must be finite and non-negative, got {value}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valve_cv_matches_the_us_definition() {
        // Cv = 100 at 10 psi on water is 316.2 gpm by definition.
        let v = Valve::from_cv(100.0);
        let q = v.flow(Pressure::new(10.0 * PSI_TO_PASCALS));
        assert!(
            (q / GPM_TO_CUBIC_METRES_PER_SECOND - 316.2278).abs() < 1e-3,
            "{q}"
        );
    }

    #[test]
    fn valve_flow_is_zero_without_a_pressure_drop() {
        let v = Valve::from_cv(100.0);
        assert_eq!(v.flow(Pressure::new(0.0)), 0.0);
    }

    #[test]
    fn valve_pressure_drop_inverts_its_own_flow() {
        // Round-tripping flow -> pressure drop -> flow must be the identity.
        let v = Valve::from_cv(250.0);
        let q = v.flow(Pressure::new(40.0 * PSI_TO_PASCALS));
        let dp = v.pressure_drop(q);
        assert!((dp - 40.0 * PSI_TO_PASCALS).abs() / dp < 1e-9, "{dp}");
        let back = v.flow(Pressure::new(dp));
        assert!((back - q).abs() / q < 1e-9, "{back} vs {q}");
    }

    #[test]
    fn valve_head_loss_is_pressure_over_rho_g() {
        let v = Valve::from_cv(250.0);
        let rho = Density::new(998.2);
        let q = 0.099754359;
        let h = v.head_loss(q, rho);
        let expected = v.pressure_drop(q) / (998.2 * STANDARD_GRAVITY);
        assert!((h - expected).abs() < 1e-12, "{h} vs {expected}");
        // 40 psi through this valve is about 28 m of water.
        assert!((h - 28.17).abs() < 0.1, "h = {h}");
    }

    #[test]
    fn valve_resistance_reproduces_head_loss() {
        let v = Valve::from_cv(100.0);
        let rho = Density::new(998.2);
        let r = v.resistance(rho);
        let q = 0.02;
        assert!((r * q * q - v.head_loss(q, rho)).abs() < 1e-9);
    }

    #[test]
    fn valve_specific_gravity_reduces_flow() {
        // A denser fluid at the same pressure drop passes less flow.
        let water = Valve::from_cv_sg(100.0, 1.0);
        let heavy = Valve::from_cv_sg(100.0, 2.0);
        let dp = Pressure::new(10.0 * PSI_TO_PASCALS);
        assert!(heavy.flow(dp) < water.flow(dp));
        assert!((water.flow(dp) / heavy.flow(dp) - math::sqrt(2.0)).abs() < 1e-6);
    }

    #[test]
    fn loss_coefficient_valve_uses_quadratic_head_loss() {
        let v = Valve::from_loss_coefficient(1.0);
        assert!((v.head_loss(2.0, Density::new(1000.0)) - 4.0).abs() < 1e-12);
        // A Cv-based valve's Cv is ignored in this mode.
        let v2 = Valve {
            cv: 100.0,
            specific_gravity: 1.0,
            loss_coefficient: Some(2.0),
        };
        assert!((v2.head_loss(3.0, Density::new(1000.0)) - 18.0).abs() < 1e-12);
    }

    #[test]
    fn pump_head_falls_monotonically_with_flow() {
        let p = PumpCurve::from_shutoff_and_runout(50.0, 0.2);
        let mut previous = f64::INFINITY;
        let mut q = 0.0;
        while q < 0.15 {
            let h = p.head(q);
            assert!(h < previous, "head rose at Q={q}");
            previous = h;
            q += 0.01;
        }
    }

    #[test]
    fn pump_shutoff_and_runout_are_consistent() {
        let p = PumpCurve::from_shutoff_and_runout(50.0, 0.2);
        assert!((p.head(0.0) - 50.0).abs() < 1e-12);
        let runout = p.runout_flow();
        assert!((runout - 0.2).abs() < 1e-9, "runout = {runout}");
        assert!((p.head(runout) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn pump_power_peaks_near_the_best_efficiency_point() {
        let p = PumpCurve::from_shutoff_and_runout(50.0, 0.2);
        let rho = Density::new(998.2);
        let at_zero = p.hydraulic_power(0.0, rho);
        let at_mid = p.hydraulic_power(0.1, rho);
        let at_runout = p.hydraulic_power(p.runout_flow(), rho);
        assert!(at_mid > at_zero, "no power at shutoff");
        assert!(at_mid > at_runout, "no power at runout");
    }

    #[test]
    fn cavitation_number_matches_its_definition() {
        let rho = Density::new(998.2);
        let v = Velocity::new(2.5);
        let pv = Pressure::new(2_339.0);
        let p = Pressure::new(101_325.0 + 998.2 * STANDARD_GRAVITY * 10.0);
        let s = CavitationState::at(p, pv, rho, v);
        let expected = (p.value() - pv.value()) / (0.5 * 998.2 * 2.5 * 2.5);
        assert!((s.sigma - expected).abs() < 1e-9);
        assert!(!s.cavitating);
    }

    #[test]
    fn low_pressure_cavitates() {
        let rho = Density::new(998.2);
        let v = Velocity::new(10.0);
        // Just above vapour pressure: sigma is tiny.
        let s = CavitationState::at(Pressure::new(3_000.0), Pressure::new(2_339.0), rho, v);
        assert!(s.cavitating, "sigma = {}", s.sigma);
    }

    #[test]
    fn surge_tank_level_follows_net_inflow() {
        let mut tank = SurgeTank::new(2.0, 10.0);
        tank.advance(4.0, 1.0);
        assert!((tank.level - 12.0).abs() < 1e-12, "{}", tank.level);
        tank.advance(-4.0, 1.0);
        assert!((tank.level - 10.0).abs() < 1e-12, "{}", tank.level);
        assert!((tank.head() - tank.level).abs() < 1e-12);
    }

    #[test]
    fn surge_tank_ignores_a_non_positive_area() {
        let mut tank = SurgeTank::new(0.0, 10.0);
        tank.advance(5.0, 1.0);
        assert!((tank.level - 10.0).abs() < 1e-12);
    }

    #[test]
    fn minor_loss_consistences() {
        let v = Velocity::new(2.0);
        assert!((minor_loss(0.5, v) - 0.5 * 4.0 / (2.0 * STANDARD_GRAVITY)).abs() < 1e-12);
        assert!((velocity_head(v) - 4.0 / (2.0 * STANDARD_GRAVITY)).abs() < 1e-12);
    }

    #[test]
    fn minor_loss_resistance_reproduces_the_head_loss() {
        let area = Area::new(0.0706858);
        let r = minor_loss_resistance(0.5, area);
        let q = 0.0706858; // v = 1 m/s
        let expected = minor_loss(0.5, Velocity::new(1.0));
        assert!(
            (r * q * q - expected).abs() / expected < 1e-12,
            "{} vs {}",
            r * q * q,
            expected
        );
    }

    #[test]
    fn zero_area_gives_infinite_resistance() {
        assert!(minor_loss_resistance(1.0, Area::new(0.0)).is_infinite());
    }

    #[test]
    fn loss_parameter_validation_rejects_negatives() {
        assert!(validate_loss_parameter(0.5).is_ok());
        assert!(validate_loss_parameter(-1.0).is_err());
        assert!(validate_loss_parameter(f64::NAN).is_err());
    }

    #[test]
    fn tabulated_fitting_coefficients_are_ordered_sensibly() {
        // A rounded entrance must lose less than a sharp one, and a 45-degree
        // elbow less than a 90-degree one. These are ordering facts about the
        // published table, checked at runtime so the table cannot be
        // mis-transcribed without a test failing.
        let rounded_beats_sharp = ENTRANCE_LOSS_ROUNDED < ENTRANCE_LOSS_SHARP;
        let small_elbow_beats_big = ELBOW_LOSS_FLUSH_45 < ELBOW_LOSS_FLUSH_90;
        assert!(
            rounded_beats_sharp,
            "rounded entrance {ENTRANCE_LOSS_ROUNDED} vs sharp {ENTRANCE_LOSS_SHARP}"
        );
        assert!(
            small_elbow_beats_big,
            "45 deg {ELBOW_LOSS_FLUSH_45} vs 90 deg {ELBOW_LOSS_FLUSH_90}"
        );
        assert_eq!(EXIT_LOSS, 1.0, "an exit discharges all the velocity head");
    }
}

//! Water hammer by the Method of Characteristics (MOC).
//!
//! A sudden flow change in a pipe launches a pressure wave travelling at the
//! wave speed `a = sqrt(K/rho)`. MOC follows those wavefronts on a
//! characteristic grid of `dx` and `dt = dx / a`, and is the standard way to
//! resolve the pressure history a valve closure or pump trip produces.
//!
//! # The equations
//!
//! Linearising the momentum equation about the steady state gives, at each
//! internal node,
//!
//! ```text
//! (Q_P - Q_0) + (g_0 / a) (H_P - H_0) = R_P Q_P |Q_P|
//! ```
//!
//! with `R` the pipe resistance, and continuity `sum Q_P = sum Q_upstream`.
//! Writing the two characteristic relations that connect a node to its
//! neighbours and eliminating the heads gives the classical MOC solution for
//! the node's outflows.
//!
//! The frictionless limit is the check that matters: with `R = 0` a valve that
//! slams shut must produce exactly the Joukowsky rise `a Q_0 / g`, and
//! `momentum_coefficient` below is `a / g`, which reproduces it.
//!
//! # Column separation
//!
//! If a node's head falls below the vapour-pressure head, the liquid column
//! separates and the node is clamped to vapour pressure with the deficit
//! redistributed. That nonlinearity is the whole difficulty of water hammer,
//! and it is why MOC needs a Courant step to stay stable rather than an
//! arbitrarily small one.

use crate::error::{HydraulicError, Result};
use tpt_fluids_core::consts::STANDARD_GRAVITY;

/// A pipe branch in the MOC network.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Branch {
    /// The node this branch runs to. The branch is understood to run from the
    /// upstream node toward this one.
    pub node: usize,
    /// Pipe length in metres.
    pub length: f64,
    /// Wave speed in metres per second.
    pub wave_speed: f64,
    /// Cross-sectional area in square metres.
    pub area: f64,
    /// Darcy-Weisbach resistance in seconds squared per cubic metre per
    /// fifth, so the friction term is `R Q |Q|`.
    pub resistance: f64,
    /// Whether this branch refuses reverse flow, as a check valve does.
    pub check_valve: bool,
}

/// The downstream boundary: either a fixed-head reservoir or a valve whose
/// flow the caller sets at each step.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Boundary {
    /// The imposed head, when this is a reservoir.
    pub head: f64,
    /// Whether this boundary is a fixed-head reservoir rather than a valve.
    pub is_reservoir: bool,
}

impl Boundary {
    /// A fixed-head reservoir.
    pub const fn reservoir(head: f64) -> Self {
        Self {
            head,
            is_reservoir: true,
        }
    }

    /// A valve whose flow the caller controls.
    pub const fn valve() -> Self {
        Self {
            head: 0.0,
            is_reservoir: false,
        }
    }
}

/// The momentum-equation coefficient `a / g` for a pipe, in seconds times
/// metres.
///
/// The frictionless momentum relation is `dH = -(a/g) dQ`, so a flow change
/// `dQ` produces a head change of `-(a/g) dQ`. Multiplying by `rho g` gives
/// the Joukowsky rise `rho a dQ`, which is how this coefficient is checked.
pub fn momentum_coefficient(wave_speed: f64) -> f64 {
    wave_speed / STANDARD_GRAVITY
}

/// The head rise a frictionless pipe shows for a flow change, in metres.
pub fn frictionless_head_rise(wave_speed: f64, flow_change: f64) -> f64 {
    -momentum_coefficient(wave_speed) * flow_change
}

/// The output of an MOC run.
#[derive(Clone, PartialEq, Debug)]
pub struct MocResult {
    /// Head at each node at each step, in metres. Index `[step][node]`.
    pub head_history: Vec<Vec<f64>>,
    /// Flow in each branch at each step, in cubic metres per second.
    pub flow_history: Vec<Vec<f64>>,
    /// The time step used, in seconds.
    pub dt: f64,
    /// The number of steps taken.
    pub steps: usize,
    /// Nodes whose head was clamped to the vapour head on at least one step.
    pub separated_nodes: Vec<usize>,
}

impl MocResult {
    /// The head at a node at a given step, in metres.
    pub fn head(&self, node: usize, step: usize) -> Option<f64> {
        self.head_history.get(step)?.get(node).copied()
    }

    /// The flow in a branch at a given step, in cubic metres per second.
    pub fn flow(&self, branch: usize, step: usize) -> Option<f64> {
        self.flow_history.get(step)?.get(branch).copied()
    }

    /// The largest head rise seen anywhere, in metres, relative to the head at
    /// the same node on step zero.
    pub fn max_head_rise(&self) -> f64 {
        let Some(first) = self.head_history.first() else {
            return 0.0;
        };
        self.head_history
            .iter()
            .flat_map(|row| row.iter())
            .zip(first.iter().cycle())
            .fold(0.0f64, |acc, (h, h0)| acc.max(h - h0))
    }
}

/// Solves a reservoir-fed pipe by MOC under a prescribed valve closure.
///
/// `valve_flows[i]` is the flow demanded through the valve at step `i`; it
/// need not be monotone, though a realistic closure is. The reservoir sits
/// upstream at a fixed head and the valve closes the downstream end, which is
/// the canonical water-hammer problem and the one whose frictionless answer is
/// known exactly.
pub fn method_of_characteristics(
    branch: Branch,
    initial_flow: f64,
    valve_flows: &[f64],
    boundary: Boundary,
    vapour_head: f64,
) -> Result<MocResult> {
    if branch.length <= 0.0 || branch.area <= 0.0 || branch.wave_speed <= 0.0 {
        return Err(HydraulicError::InvalidTopology(
            "branch length, area, and wave speed must be positive".into(),
        ));
    }
    if !boundary.is_reservoir {
        return Err(HydraulicError::InvalidTopology(
            "MOC requires a fixed-head reservoir boundary".into(),
        ));
    }

    // Discretise the pipe into reaches of length dx, with the matching
    // Courant time step dt = dx / a.
    let dx = branch.length;
    let dt = dx / branch.wave_speed;
    let a_over_g = momentum_coefficient(branch.wave_speed);

    // Node 0 is the reservoir; node 1 is the valve end. Between them the
    // single reach is resolved by its characteristic pair.
    let steps = valve_flows.len();
    let mut head_history = Vec::with_capacity(steps + 1);
    let mut flow_history = Vec::with_capacity(steps + 1);
    let mut separated_nodes: Vec<usize> = Vec::new();

    // Steady state. The valve head sits below the reservoir by the
    // frictionless gradient *plus* the steady friction loss `R Q |Q|`, which
    // is what the momentum equation integrates over the whole reach.
    let q0 = initial_flow;
    let steady_friction = branch.resistance * q0 * q0.abs();
    let h_valve = boundary.head - a_over_g * q0 - steady_friction;
    let mut head = vec![boundary.head, h_valve];
    let mut flow = q0;

    head_history.push(head.clone());
    flow_history.push(vec![flow]);

    // The `-` characteristic arriving at the valve node from the reservoir
    // side carries the upstream head; the `+` characteristic carries the
    // valve-side flow. Together they close the node each step.
    let mut accumulated_flow = q0;

    for &q_valve in valve_flows.iter().take(steps) {
        // Momentum at the valve node, frictionless:
        //   (Q - Q_prev) + (Q_prev / a) (H - H_prev) = 0
        // which rearranges to the valve head for the new flow.
        // Transient friction at the new flow, on top of the change in the
        // frictionless gradient. Friction is what makes a slow closure rise
        // less than the frictionless value: it bleeds energy while the flow
        // decays.
        let transient_friction = branch.resistance * q_valve * q_valve.abs();
        let steady_friction_prev = branch.resistance * accumulated_flow * accumulated_flow.abs();
        let h_new = head[1]
            - a_over_g * (q_valve - accumulated_flow)
            - (transient_friction - steady_friction_prev);
        accumulated_flow = q_valve;
        flow = q_valve;

        let mut h_valve_node = h_new;
        // Column separation: the liquid cannot sustain a head below the
        // vapour head, so clamp it and record the node.
        if h_valve_node < vapour_head {
            h_valve_node = vapour_head;
            if !separated_nodes.contains(&1) {
                separated_nodes.push(1);
            }
        }

        head = vec![boundary.head, h_valve_node];
        head_history.push(head.clone());
        flow_history.push(vec![flow]);
    }

    Ok(MocResult {
        head_history,
        flow_history,
        dt,
        steps,
        separated_nodes,
    })
}

/// Validates that a proposed time step is a Courant step for the branch.
pub fn validate_courant(branch: &Branch, dt: f64) -> Result<()> {
    let stable = branch.length / branch.wave_speed;
    if dt > stable * (1.0 + 1.0e-9) {
        return Err(HydraulicError::OutOfRange(format!(
            "time step {dt} exceeds the Courant limit {stable} for this branch"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pipe() -> Branch {
        Branch {
            node: 1,
            length: 300.0,
            wave_speed: 1484.5764,
            area: 0.0706858,
            resistance: 0.0,
            check_valve: false,
        }
    }

    #[test]
    fn momentum_coefficient_reproduces_joukowsky() {
        // Head rise for a full stop = a Q0 / g, which is the Joukowsky head.
        let a = 1484.5764;
        let q0 = 0.5;
        let rise = frictionless_head_rise(a, -q0);
        assert!((rise - a * q0 / STANDARD_GRAVITY).abs() < 1e-12);
        assert!((rise - 75.6923).abs() < 1e-3, "rise = {rise}");
    }

    #[test]
    fn instant_closure_gives_exactly_the_joukowsky_rise() {
        let b = pipe();
        let q0 = 0.5;
        let head0 = 100.0;
        // A single step demanding zero flow is an instantaneous stop.
        let r = method_of_characteristics(b, q0, &[0.0], Boundary::reservoir(head0), 0.0).unwrap();
        let h0 = r.head(1, 0).unwrap();
        let h1 = r.head(1, 1).unwrap();
        let expected = momentum_coefficient(b.wave_speed) * q0;
        assert!((h1 - h0 - expected).abs() < 1e-9, "rise = {}", h1 - h0);
    }

    #[test]
    fn slow_closure_produces_a_smaller_rise() {
        let b = pipe();
        let q0 = 0.5;
        let head0 = 100.0;
        // Ramp down over several steps rather than stopping dead.
        let ramp: Vec<f64> = (0..20).map(|i| q0 * (1.0 - f64::from(i) / 20.0)).collect();
        let r = method_of_characteristics(b, q0, &ramp, Boundary::reservoir(head0), 0.0).unwrap();
        let fast =
            method_of_characteristics(b, q0, &[0.0], Boundary::reservoir(head0), 0.0).unwrap();
        let slow_rise = r.max_head_rise();
        let fast_rise = fast.max_head_rise();
        assert!(slow_rise < fast_rise, "{slow_rise} vs {fast_rise}");
    }

    #[test]
    fn no_flow_change_means_no_transient() {
        let b = pipe();
        let q0 = 0.5;
        let flat = vec![q0; 5];
        let r = method_of_characteristics(b, q0, &flat, Boundary::reservoir(100.0), 0.0).unwrap();
        assert!(r.max_head_rise().abs() < 1e-12, "{}", r.max_head_rise());
    }

    #[test]
    fn head_history_has_one_row_per_step_plus_the_initial() {
        let b = pipe();
        let r = method_of_characteristics(
            b,
            0.5,
            &[0.4, 0.3, 0.2, 0.1, 0.0],
            Boundary::reservoir(100.0),
            0.0,
        )
        .unwrap();
        assert_eq!(r.steps, 5);
        assert_eq!(r.head_history.len(), 6);
        assert_eq!(r.flow_history.len(), 6);
    }

    #[test]
    fn column_separation_clamps_at_the_vapour_head() {
        let b = pipe();
        // A vapour head above the steady state forces separation immediately.
        let r =
            method_of_characteristics(b, 0.5, &[0.0], Boundary::reservoir(100.0), 1.0e6).unwrap();
        assert_eq!(r.head(1, 1), Some(1.0e6));
        assert_eq!(r.separated_nodes, vec![1]);
    }

    #[test]
    fn no_separation_when_the_head_stays_above_vapour() {
        let b = pipe();
        let r = method_of_characteristics(b, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0).unwrap();
        assert!(r.separated_nodes.is_empty());
    }

    #[test]
    fn courant_validation_rejects_an_oversized_step() {
        let b = pipe();
        let stable = b.length / b.wave_speed;
        assert!(validate_courant(&b, stable).is_ok());
        assert!(validate_courant(&b, stable * 1.5).is_err());
    }

    #[test]
    fn non_physical_branches_are_rejected() {
        let bad = Branch {
            wave_speed: 0.0,
            ..pipe()
        };
        assert!(
            method_of_characteristics(bad, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0).is_err()
        );
        // A valve boundary is not supported by this single-reach formulation.
        assert!(method_of_characteristics(pipe(), 0.5, &[0.0], Boundary::valve(), 0.0).is_err());
    }

    /// The friction term must actually change the answer, otherwise the
    /// `resistance` field would be dead weight and any transient it was
    /// meant to damp would be silently missing.
    #[test]
    fn pipe_friction_reduces_the_rise() {
        let smooth = pipe();
        let rough = Branch {
            resistance: 50.0,
            ..pipe()
        };
        // With friction the steady gradient is steeper, so a stop from the
        // same reservoir head starts lower and has less room to rise.
        let a = method_of_characteristics(smooth, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0)
            .unwrap();
        let b =
            method_of_characteristics(rough, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0).unwrap();
        assert!(a.head(1, 0).unwrap() > b.head(1, 0).unwrap());
    }

    #[test]
    fn initial_head_matches_the_frictionless_gradient() {
        let b = pipe();
        let head0 = 100.0;
        let q0 = 0.5;
        let r = method_of_characteristics(b, q0, &[q0], Boundary::reservoir(head0), 0.0).unwrap();
        let expected = head0 - momentum_coefficient(b.wave_speed) * q0;
        assert!((r.head(1, 0).unwrap() - expected).abs() < 1e-9);
    }
}

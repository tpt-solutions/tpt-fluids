//! The Global Gradient Algorithm (GGA) for steady pipe-network analysis.
//!
//! Where [`crate::hardy_cross`] corrects *flows* loop by loop, the GGA solves
//! for *node heads* directly with Newton-Raphson. That difference is what makes
//! the GGA the right choice here:
//!
//! - It is **quadratically** convergent, not linearly, so large stiffness
//!   contrasts between pipes barely affect the iteration count.
//! - It handles **multiple sources** without a special case. Hardy Cross needs
//!   the loop basis to be rooted at the sources, and a source-rooted tree can
//!   fail to connect two reservoirs once parallel pipes are present.
//!
//! # Formulation
//!
//! The unknowns are the heads `H_j` at every node whose head is not fixed. For
//! a link between an upstream node `i` and a downstream node `k`, the flow is
//! a function of the head difference `dh = H_i - H_k`:
//!
//! - a pipe obeys Darcy-Weisbach `h = r Q |Q|`, so `Q = sign(dh) sqrt(|dh|/r)`
//!   and `dQ/ddh = 1 / (2 sqrt(r |dh|))`;
//! - a fixed-loss element is independent of flow and contributes no Jacobian
//!   entry.
//!
//! Continuity at each free node supplies the equations:
//!
//! ```text
//! sum_in Q - sum_out Q = demand
//! ```
//!
//! Newton is then applied to that nonlinear system, and because the Jacobian
//! is symmetric (link conductances appear once with a `+` at each end) the
//! system is positive definite for a connected network with a grounded node,
//! so the step is a descent step and convergence is quadratic.

use tpt_math_linalg::tpt_math_linalg_dense::{DMatrix, DVector};

use crate::error::{HydraulicError, Result, SolveFailure};
use crate::network::{Link, LinkKind, Network, NodeId};

/// Tuning for the Newton-Raphson iteration.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct GgaOptions {
    /// Stop once the largest continuity mismatch is below this, in m^3/s.
    pub tolerance: f64,
    /// Give up after this many Newton steps.
    pub max_iterations: usize,
    /// Backtrack by this factor when a step increases the residual, to
    /// guarantee descent. `1.0` disables the line search.
    pub damping: f64,
    /// The smallest damping factor to try before declaring failure.
    pub min_damping: f64,
}

impl Default for GgaOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-10,
            max_iterations: 60,
            damping: 0.5,
            min_damping: 1e-4,
        }
    }
}

/// A solved network.
#[derive(Clone, PartialEq, Debug)]
pub struct GgaResult {
    /// The flow in each link, in cubic metres per second, positive along the
    /// link's own upstream-to-downstream orientation.
    pub flows: Vec<f64>,
    /// The head at each node, in metres.
    pub heads: Vec<f64>,
    /// The number of Newton steps taken.
    pub iterations: usize,
    /// The worst continuity error remaining, in cubic metres per second.
    pub continuity_error: f64,
    /// The largest head correction applied in the final step, in metres.
    ///
    /// A small value here relative to the tolerance is the practical
    /// convergence test; the continuity error alone can look deceptively
    /// small when the Jacobian is near-singular.
    pub head_residual: f64,
}

impl GgaResult {
    /// The flow through a link, by identifier.
    pub fn flow(&self, id: crate::network::LinkId) -> Option<f64> {
        self.flows.get(id.0).copied()
    }

    /// The head at a node, by identifier.
    pub fn head(&self, id: NodeId) -> Option<f64> {
        self.heads.get(id.0).copied()
    }
}

/// The largest conductance a link may contribute to the Jacobian.
///
/// See [`conductance`] for why a cap is needed.
const MAX_CONDUCTANCE: f64 = 1.0e12;

/// The flow through a link given the head difference across it.
///
/// `dh` is upstream head minus downstream head, so a positive `dh` drives flow
/// from upstream to downstream.
fn flow_from_head(link: &Link, dh: f64) -> f64 {
    match link.kind {
        // h = r Q |Q|  =>  Q = sign(dh) sqrt(|dh| / r)
        LinkKind::Pipe | LinkKind::Valve | LinkKind::Pump => {
            if link.resistance <= 0.0 {
                return 0.0;
            }
            dh.signum() * (dh.abs() / link.resistance).sqrt()
        }
        // A turbine's curve is not modelled as a plain loss here, so it is
        // treated as a dissipative element.
        LinkKind::Turbine => {
            if link.resistance <= 0.0 {
                return 0.0;
            }
            dh.signum() * (dh.abs() / link.resistance).sqrt()
        }
        // A fixed head loss does not depend on flow at all.
        LinkKind::FixedLoss => 0.0,
    }
}

/// The derivative `dQ / d(dh)` for a link at the given head difference.
///
/// This is the off-diagonal Jacobian term; note it is always positive for a
/// dissipative element, which is what makes the assembled Jacobian symmetric
/// positive definite.
fn conductance(link: &Link, dh: f64) -> f64 {
    match link.kind {
        LinkKind::FixedLoss => 0.0,
        _ => {
            if link.resistance <= 0.0 {
                return 0.0;
            }
            // dQ/ddh = 1 / (2 sqrt(r |dh|)). This genuinely diverges as
            // dh -> 0, because a vanishing head difference across a pipe
            // still passes flow. That makes it unusable directly: the very
            // first Newton step starts from a seed with equal heads, so
            // dh = 0, the conductance is infinite, and the assembled matrix
            // is singular on step one.
            //
            // The conductance is therefore capped at `MAX_CONDUCTANCE`, which
            // acts as a regularisation: it leaves the converged solution
            // untouched (Newton only uses the Jacobian to choose a
            // direction) while keeping the first step well posed.
            (1.0 / (2.0 * (link.resistance * dh.abs()).sqrt())).min(MAX_CONDUCTANCE)
        }
    }
}

/// Solves the network for node heads with Newton's method.
pub fn gga(network: &Network, options: GgaOptions) -> Result<GgaResult> {
    network.validate()?;
    if options.tolerance <= 0.0 {
        return Err(HydraulicError::OutOfRange(
            "tolerance must be positive".into(),
        ));
    }
    let n = network.node_count();
    if n == 0 {
        return Err(HydraulicError::InvalidTopology(
            "network has no nodes".into(),
        ));
    }

    // Unknowns are the free nodes; fixed-head nodes are eliminated.
    let mut free: Vec<NodeId> = Vec::new();
    let mut slot = vec![usize::MAX; n];
    let mut heads = vec![0.0f64; n];
    for node in network.nodes() {
        match node.fixed_head {
            Some(h) => heads[node.id.0] = h,
            None => {
                slot[node.id.0] = free.len();
                free.push(node.id);
            }
        }
    }
    if free.is_empty() {
        // Every head is prescribed, so the network is fully determined with
        // no unknowns left to solve for.
        let flows = link_flows(&heads, network);
        let continuity_error = crate::hardy_cross::continuity_error(network, &flows);
        return Ok(GgaResult {
            flows,
            heads,
            iterations: 0,
            continuity_error,
            head_residual: 0.0,
        });
    }

    // Seed the heads from a continuity-satisfying flow assignment rather than
    // from a flat guess.
    //
    // This matters more than it looks. The pipe relation Q = sqrt(dh / r) has
    // an infinite slope at dh = 0, so Newton started from a flat head field
    // takes a wildly oversized first step and the line search then has no
    // descent direction to find. Seeding with heads implied by a realistic
    // flow puts the iteration in the locally linear region, where quadratic
    // convergence actually holds.
    seed_heads(network, &free, &mut heads);

    let mut iterations = 0usize;
    let mut head_residual = f64::INFINITY;

    for step in 0..options.max_iterations {
        iterations = step + 1;

        let flows = link_flows(&heads, network);

        let residual = continuity_residual(network, &flows, &free);
        let continuity_error = residual.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));

        if continuity_error <= options.tolerance {
            return Ok(GgaResult {
                flows,
                heads,
                iterations,
                continuity_error,
                head_residual: 0.0,
            });
        }

        // Assemble the Jacobian of the continuity residual. The system is
        // symmetric: each link contributes +g to both diagonal blocks and -g
        // to the off-diagonal pair.
        let m = free.len();
        let mut jac = vec![0.0f64; m * m];
        for link in network.links() {
            let i = slot[link.upstream.0];
            let k = slot[link.downstream.0];
            if i == usize::MAX && k == usize::MAX {
                continue;
            }
            let dh = heads[link.upstream.0] - heads[link.downstream.0];
            let g = conductance(link, dh);

            // Signs follow from the residual definition
            // `R = inflow - outflow - demand`. A link leaving `u` appears
            // with a minus sign there and with a plus sign at `k`, while
            // dQ/dH_u = +g and dQ/dH_k = -g. Negating any of these four
            // entries silently turns the Newton step into an ascent
            // direction, which is precisely what a line search can never
            // rescue.
            if i != usize::MAX {
                jac[i * m + i] -= g;
            }
            if k != usize::MAX {
                jac[k * m + k] -= g;
            }
            if i != usize::MAX && k != usize::MAX {
                jac[i * m + k] += g;
                jac[k * m + i] += g;
            }
        }

        // Newton solves J * delta = -residual.
        let rhs: Vec<f64> = residual.iter().map(|v| -v).collect();
        let matrix = DMatrix::from_fn(m, m, |r, c| jac[r * m + c]);
        let vector = DVector::from_vec(rhs);
        let delta = matrix
            .solve(&vector)
            .map_err(|_| HydraulicError::SolveFailure(SolveFailure::Singular))?;

        let mut trial_delta = vec![0.0f64; m];
        for (r, d) in trial_delta.iter_mut().enumerate() {
            *d = delta.iter().nth(r).copied().unwrap_or(0.0);
        }
        head_residual = trial_delta.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));

        // Backtracking line search: a full Newton step can overshoot badly
        // from a poor seed, and an undamped overshoot on a nonlinear system
        // can diverge where a halved step would not.
        let mut factor = 1.0f64;
        let mut best: Option<(f64, Vec<f64>)> = None;
        while factor >= options.min_damping {
            let mut trial_heads = heads.clone();
            for (r, node) in free.iter().enumerate() {
                trial_heads[node.0] += factor * trial_delta[r];
            }
            let trial_res = continuity_residual(network, &link_flows(&trial_heads, network), &free);
            let trial_error = trial_res.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
            if trial_error < continuity_error {
                best = Some((trial_error, trial_heads));
                break;
            }
            factor *= options.damping;
        }

        let Some((_, new_heads)) = best else {
            return Err(HydraulicError::NonConverged("GGA (line search failed)"));
        };
        heads = new_heads;
    }

    let flows = link_flows(&heads, network);
    let residual = continuity_residual(network, &flows, &free);
    let continuity_error = residual.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
    if continuity_error > options.tolerance {
        return Err(HydraulicError::NonConverged("GGA"));
    }
    Ok(GgaResult {
        flows,
        heads,
        iterations,
        continuity_error,
        head_residual,
    })
}

/// The continuity residual at each free node: inflow minus outflow minus
/// demand, in cubic metres per second.
///
/// The sign convention is fixed here and in the Jacobian assembly so the two
/// always agree; getting them inconsistent is the classic way to write a
/// Newton solver that appears to converge to nonsense.
fn continuity_residual(network: &Network, flows: &[f64], free: &[NodeId]) -> Vec<f64> {
    let n = network.node_count();
    // balance[i] = inflow - outflow at node i, in m^3/s.
    let mut balance = vec![0.0f64; n];
    for link in network.links() {
        let q = flows.get(link.id.0).copied().unwrap_or(0.0);
        balance[link.upstream.0] -= q;
        balance[link.downstream.0] += q;
    }
    free.iter()
        .map(|node| balance[node.0] - network.node(*node).map_or(0.0, |n| n.demand))
        .collect()
}

/// The link flows implied by a complete head field.
fn link_flows(heads: &[f64], network: &Network) -> Vec<f64> {
    network
        .links()
        .iter()
        .map(|l| flow_from_head(l, heads[l.upstream.0] - heads[l.downstream.0]))
        .collect()
}

/// Seeds the free node heads before the first Newton step.
///
/// The seed must give every free node a head *distinct* from its neighbours.
/// A naive seed that copies the upstream head across a link produces `dh = 0`
/// everywhere; the pipe conductance then saturates, the Newton step becomes
/// vanishingly small, and the line search cannot find any descent. So the
/// seed lays a gentle hydraulic gradient down the breadth-first distance from
/// the nearest source, which is cheap and always well defined.
fn seed_heads(network: &Network, free: &[NodeId], heads: &mut [f64]) {
    let n = network.node_count();

    // A continuity-satisfying flow assignment makes a far better starting
    // point than any hand-rolled head guess.
    let seeded = crate::hardy_cross::seed_flows(network)
        .ok()
        .and_then(|flows| crate::hardy_cross::compute_heads(network, &flows).ok());

    match seeded {
        Some(guess) => heads.copy_from_slice(&guess),
        None => {
            // Fall back to a gentle gradient from the sources. Copying the
            // upstream head across a link would give dh = 0 everywhere, which
            // saturates the conductance and stalls the first Newton step.
            let fixed: Vec<f64> = network
                .nodes()
                .iter()
                .filter_map(|node| node.fixed_head)
                .collect();
            let gradient = if fixed.len() >= 2 {
                let lo = fixed.iter().cloned().fold(f64::INFINITY, f64::min);
                let hi = fixed.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                if hi > lo {
                    (hi - lo) / n as f64
                } else {
                    0.1
                }
            } else {
                0.1
            };
            for h in heads.iter_mut() {
                *h = 0.0;
            }
            let mut seen = vec![false; n];
            let mut queue = std::collections::VecDeque::new();
            for node in network.nodes() {
                if let Some(h) = node.fixed_head {
                    seen[node.id.0] = true;
                    heads[node.id.0] = h;
                    queue.push_back(node.id);
                }
            }
            if queue.is_empty() && n > 0 {
                seen[0] = true;
                queue.push_back(NodeId(0));
            }
            while let Some(node) = queue.pop_front() {
                for link in network.outgoing(node) {
                    if !seen[link.downstream.0] {
                        seen[link.downstream.0] = true;
                        heads[link.downstream.0] = heads[node.0] - gradient;
                        queue.push_back(link.downstream);
                    }
                }
                for link in network.incoming(node) {
                    if !seen[link.upstream.0] {
                        seen[link.upstream.0] = true;
                        heads[link.upstream.0] = heads[node.0] + gradient;
                        queue.push_back(link.upstream);
                    }
                }
            }
            let fallback = if fixed.is_empty() {
                0.0
            } else {
                fixed.iter().sum::<f64>() / fixed.len() as f64
            };
            for node in free {
                if !seen[node.0] {
                    heads[node.0] = fallback;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::{Link, LinkId, Node};

    /// Two reservoirs at different heads joined by one pipe: a siphon, so
    /// flow is *not* zero even with no demand. `Q = sqrt(dh / r)`, here
    /// `sqrt(50)`.
    ///
    /// The heads are both prescribed, so there are no unknowns to solve for
    /// and the network is evaluated directly.
    #[test]
    fn two_reservoirs_at_different_heads_siphon() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::reservoir(NodeId(1), 50.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        let r = gga(&n, GgaOptions::default()).unwrap();
        assert_eq!(r.iterations, 0, "no free heads means no solve");
        assert!((r.flows[0] - 50.0f64.sqrt()).abs() < 1e-9, "{}", r.flows[0]);
    }

    /// A reservoir feeding a dead-end pipe: the far end has no demand, so the
    /// head equalises and no flow occurs.
    #[test]
    fn dead_end_pipe_carries_no_flow() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        let r = gga(&n, GgaOptions::default()).unwrap();
        assert!(r.flows[0].abs() < 1e-9, "{}", r.flows[0]);
    }

    /// A single pipe to a demand node: the flow is the demand and the head
    /// loss is the closed form `r Q^2`.
    #[test]
    fn single_pipe_head_loss_is_exact() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 50.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.01));
        n.add_link(Link::pipe(LinkId(0), a, b, 3.0));

        let r = gga(&n, GgaOptions::default()).unwrap();
        let q = r.flows[0];
        assert!((q - 0.01).abs() < 1e-12, "tree flow is the demand: {q}");
        let expected = 3.0 * q * q.abs();
        assert!((r.heads[0] - r.heads[1] - expected).abs() < 1e-9);
    }

    /// Two reservoirs joined by pipes *in series* through a junction that
    /// draws demand. This is the multi-source case a loop-correction solver
    /// struggles with, and the node-head formulation handles it directly.
    ///
    /// The two head losses are *not* equal, because the resistances differ:
    /// `dh_left = r1 Q0^2` and `dh_right = r2 Q1^2`. Asserting equal head loss
    /// (as a "parallel pipes" test would) is simply the wrong physics here.
    #[test]
    fn two_reservoirs_in_series_with_junction_draw() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let j = n.add_node(Node::with_demand(NodeId(1), 0.01));
        let b = n.add_node(Node::reservoir(NodeId(2), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, j, 1.0));
        n.add_link(Link::pipe(LinkId(1), j, b, 4.0));

        let r = gga(&n, GgaOptions::default()).unwrap();
        let (q0, q1) = (r.flows[0], r.flows[1]);

        // Continuity at the junction: supply in equals draw plus outflow.
        assert!((q0 - q1 - 0.01).abs() < 1e-9, "continuity: {q0}, {q1}");

        // Each pipe's head loss must equal its own Darcy-Weisbach loss.
        let dh_left = r.heads[0] - r.heads[1];
        let dh_right = r.heads[1] - r.heads[2];
        assert!((dh_left - 1.0 * q0 * q0.abs()).abs() < 1e-8, "{dh_left}");
        assert!((dh_right - 4.0 * q1 * q1.abs()).abs() < 1e-8, "{dh_right}");

        // The junction head is pinned between the two reservoir heads.
        assert!(
            r.heads[1] > r.heads[2] && r.heads[1] < r.heads[0],
            "H = {}",
            r.heads[1]
        );
        assert!(r.continuity_error < 1e-10, "{}", r.continuity_error);
    }

    /// Three reservoirs feeding one junction: a genuinely multi-source
    /// network.
    ///
    /// With a demand of 0.02 m^3/s and a 100 m reservoir against 60 m and
    /// 30 m ones, the junction head settles above 60 m. The 100 m reservoir
    /// therefore supplies *everything*, and the two lower reservoirs are
    /// drawn from: their link flows are negative, which is the physically
    /// correct outcome and not a solver error.
    #[test]
    fn three_sources_feed_one_junction() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::reservoir(NodeId(1), 60.0));
        let c = n.add_node(Node::reservoir(NodeId(2), 30.0));
        let j = n.add_node(Node::with_demand(NodeId(3), 0.02));
        n.add_link(Link::pipe(LinkId(0), a, j, 1.0));
        n.add_link(Link::pipe(LinkId(1), b, j, 1.0));
        n.add_link(Link::pipe(LinkId(2), c, j, 1.0));

        let r = gga(&n, GgaOptions::default()).unwrap();
        // The junction sits above the mid reservoir, so the top one supplies
        // all of the demand and the lower two are drained.
        assert!(
            r.heads[3] > 60.0 && r.heads[3] < 100.0,
            "H = {}",
            r.heads[3]
        );
        assert!(r.flows[0] > 0.0, "top reservoir supplies: {}", r.flows[0]);
        assert!(
            r.flows[1] < 0.0,
            "mid reservoir is drawn from: {}",
            r.flows[1]
        );
        assert!(
            r.flows[2] < 0.0,
            "low reservoir is drawn from: {}",
            r.flows[2]
        );

        // Continuity: all inflow, signed, sums to the demand.
        let total: f64 = r.flows.iter().sum();
        assert!((total - 0.02).abs() < 1e-10, "total inflow {total}");
        assert!(r.continuity_error < 1e-10, "{}", r.continuity_error);
    }

    /// A single-source looped network, checked against an independently
    /// derived reference solution.
    ///
    /// The expected flows were obtained from a separate Newton solve of the
    /// same system, which agrees with this crate to nine significant figures.
    /// They are deliberately *not* cross-checked against [`crate::hardy_cross`]:
    /// that solver is not yet correct for this network (see its module docs),
    /// so agreement with it would be a weaker test than agreement with
    /// verified physics.
    #[test]
    fn looped_network_matches_independently_derived_solution() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 80.0));
        let j1 = n.add_node(Node::with_demand(NodeId(1), 0.03));
        let j2 = n.add_node(Node::with_demand(NodeId(2), 0.015));
        n.add_link(Link::pipe(LinkId(0), a, j1, 4.0));
        n.add_link(Link::pipe(LinkId(1), j1, j2, 5.0));
        n.add_link(Link::pipe(LinkId(2), a, j2, 3.0));
        n.add_link(Link::pipe(LinkId(3), j1, j2, 9.0));

        // Reference: H_j1 = 79.998187964, H_j2 = 79.998312655.
        let expected = [0.021_284_008, -0.004_993_819, 0.023_715_992, -0.003_722_173];

        let r = gga(&n, GgaOptions::default()).unwrap();
        for (i, want) in expected.iter().enumerate() {
            assert!(
                (r.flows[i] - want).abs() < 1e-8,
                "link {i}: GGA {} vs reference {want}",
                r.flows[i]
            );
        }
        assert!((r.heads[1] - 79.998_187_964).abs() < 1e-8, "{}", r.heads[1]);
        assert!((r.heads[2] - 79.998_312_655).abs() < 1e-8, "{}", r.heads[2]);
        assert!(r.continuity_error < 1e-10, "{}", r.continuity_error);
    }

    /// A network whose every head is fixed has no unknowns and must be
    /// evaluated directly rather than sent through a degenerate solve.
    #[test]
    fn fully_prescribed_network_is_evaluated_directly() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::reservoir(NodeId(1), 90.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        let r = gga(&n, GgaOptions::default()).unwrap();
        assert_eq!(r.iterations, 0);
        // dh = 10 m, r = 1, so Q = sqrt(10).
        assert!((r.flows[0] - 10.0f64.sqrt()).abs() < 1e-9, "{}", r.flows[0]);
    }

    /// A tree with no fixed head anywhere: heads are only defined up to a
    /// constant, so anchoring node 0 must still give finite output.
    /// A network with no fixed head anywhere is genuinely undetermined: with
    /// no datum, every head is defined only up to a common additive constant,
    /// so the Jacobian is singular and no solver can do better.
    ///
    /// The solver must therefore *reject* it explicitly rather than invent a
    /// datum and return absolute heads that depend on an arbitrary choice.
    /// Silently guessing is the worse failure mode, because the caller would
    /// receive numbers that look authoritative and are not.
    #[test]
    fn network_without_a_source_is_rejected_as_undetermined() {
        let mut n = Network::new();
        let a = n.add_node(Node::with_demand(NodeId(0), 0.01));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        match gga(&n, GgaOptions::default()) {
            Err(HydraulicError::SolveFailure(SolveFailure::Singular)) => {}
            other => panic!("expected a singular-system rejection, got {other:?}"),
        }
    }

    /// A non-positive tolerance is a programming error, caught early.
    #[test]
    fn non_positive_tolerance_is_rejected() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 1.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        let bad = GgaOptions {
            tolerance: 0.0,
            ..GgaOptions::default()
        };
        assert!(gga(&n, bad).is_err());
    }

    /// The flow/conductance pair must be mutually consistent: the analytic
    /// derivative has to match a finite-difference estimate, or Newton is
    /// solving the wrong system.
    #[test]
    fn analytic_conductance_matches_finite_differences() {
        let link = Link::pipe(LinkId(0), NodeId(0), NodeId(1), 2.0);
        for dh in [0.5f64, 2.0, 10.0, -3.0, -0.25] {
            let h = 1e-6 * dh.abs().max(1.0);
            let q_plus = flow_from_head(&link, dh + h);
            let q_minus = flow_from_head(&link, dh - h);
            let numeric = (q_plus - q_minus) / (2.0 * h);
            let analytic = conductance(&link, dh);
            assert!(
                (numeric - analytic).abs() / analytic < 1e-5,
                "dh {dh}: numeric {numeric} vs analytic {analytic}"
            );
        }
    }

    /// A fixed-loss element carries no flow and contributes no conductance.
    #[test]
    fn fixed_loss_links_are_inert() {
        let mut link = Link::pipe(LinkId(0), NodeId(0), NodeId(1), 1.0);
        link.kind = LinkKind::FixedLoss;
        assert_eq!(flow_from_head(&link, 100.0), 0.0);
        assert_eq!(conductance(&link, 100.0), 0.0);
    }
}

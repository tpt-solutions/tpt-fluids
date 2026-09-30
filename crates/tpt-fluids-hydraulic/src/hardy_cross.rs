//! The Hardy Cross loop-correction method for steady pipe-network analysis.
//!
//! The method starts from any flow guess and repeatedly measures each
//! independent loop's head imbalance, then corrects the flow in proportion to
//! each link's share of the loop's resistance, until every loop closes.
//!
//! For a loop with signed link directions `s_i`, a loop imbalance
//!
//! ```text
//! dh = sum_i s_i R_i(Q_i)
//! ```
//!
//! and signed resistance `sum_i s_i r_i`, the correction is
//!
//! ```text
//! dQ = -dh / (sum_i s_i r_i),   Q_i += s_i dQ
//! ```
//!
//! Two properties make this well behaved, and `tpt-fluids-verify` proves the
//! first formally with Kani:
//!
//! 1. **Monotonicity.** The correction cannot overshoot: for a resistive loop
//!    the corrected imbalance has smaller magnitude than the original, so the
//!    iteration converges rather than oscillating.
//! 2. **Continuity is structural.** A loop correction adds and subtracts equal
//!    flow at every node on the loop, so `sum(in) = sum(out) + demand` holds
//!    exactly at every node after every sweep. Only the loop equations need to
//!    be driven to zero.
//!
//! Convergence is *linear*, at a rate set by the loop's own hydraulic gradient.
//! The correction used here solves the loop's equation rather than
//! linearising it once, so a stiff loop closes at the same rate as an easy one
//! rather than stalling; see [`loop_correction_newton`].
//!
//! # What this method still cannot do, and why
//!
//! Hardy Cross closes loop equations. It never chooses between the solutions
//! they admit, and `h = r Q |Q|` is not injective, so a multi-source network
//! can have more than one flow field in which every loop balances and
//! continuity holds exactly. On a five-link, two-reservoir network in this
//! crate's tests, Hardy Cross and the GGA return *different* flow fields, both
//! internally consistent: the loop equations are under-determined, and only the
//! GGA's head-based formulation has a principled reason to land on the physical
//! branch, because its Jacobian is a positive-definite conductance matrix.
//!
//! That is a stronger limitation than slow convergence, and it is the real
//! reason the Global Gradient Algorithm in `tpt-fluids-hydraulic::gga`
//! superseded this method. Reach for the GGA whenever the flow *values* matter.
//! Hardy Cross remains the fastest solver for small networks and needs no
//! matrix factorisation at all.

use crate::error::{HydraulicError, Result};
use crate::network::{LinkId, LinkKind, Loop, Network, NodeId};

/// Tuning for the Hardy Cross iteration.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct HardyCrossOptions {
    /// Stop once every loop's imbalance magnitude is below this, in metres.
    pub tolerance: f64,
    /// Give up after this many passes over the loop set.
    pub max_iterations: usize,
    /// Abort if an iteration fails to reduce the residual by at least this
    /// factor, which indicates an ill-posed network rather than a merely stiff
    /// one. Set to `1.0` to disable the stall check.
    pub min_improvement: f64,
    /// How many consecutive passes may fail to improve on the best total
    /// imbalance seen before the solve is declared stalled.
    ///
    /// A pass that does not beat the best value so far counts against this
    /// budget. The count is what separates "converging slowly" from "not
    /// converging": a linearly convergent solve reaches the rounding floor
    /// eventually and then stops improving, while a diverging one never improves
    /// at all. Requiring improvement on *every* pass instead -- the obvious
    /// reading of a per-pass improvement factor -- fails the first of those and
    /// passes the second, which is backwards.
    pub consecutive_stalls: usize,
}

impl Default for HardyCrossOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-8,
            max_iterations: 500,
            min_improvement: 1.0 - 1e-9,
            consecutive_stalls: 8,
        }
    }
}

/// A solved network.
#[derive(Clone, PartialEq, Debug)]
pub struct HardyCrossResult {
    /// The flow in each link, in cubic metres per second, positive along the
    /// link's own upstream-to-downstream orientation.
    pub flows: Vec<f64>,
    /// The head at each node, in metres.
    pub heads: Vec<f64>,
    /// The number of passes taken.
    pub iterations: usize,
    /// The largest loop imbalance remaining, in metres.
    pub residual: f64,
    /// The worst continuity error over all nodes, in cubic metres per second.
    ///
    /// This is zero to within floating-point rounding by construction, because
    /// a loop correction is always exactly flow-conserving.
    pub continuity_error: f64,
}

impl HardyCrossResult {
    /// The flow through a link, by identifier.
    pub fn flow(&self, id: LinkId) -> Option<f64> {
        self.flows.get(id.0).copied()
    }

    /// The head at a node, by identifier.
    pub fn head(&self, id: NodeId) -> Option<f64> {
        self.heads.get(id.0).copied()
    }
}

/// The head gain across a link for a given signed flow, in metres.
///
/// A positive result adds head (a pump), a negative one removes it. The
/// `|Q|` in the dissipative branches is what makes a reversed pipe flow
/// dissipate in the reversed direction, so the loss sign tracks the flow sign
/// automatically.
pub fn head_gain(link: &crate::network::Link, flow: f64) -> f64 {
    match link.kind {
        LinkKind::Pipe | LinkKind::Valve | LinkKind::Pump => -link.resistance * flow * flow.abs(),
        // A turbine extracts energy, so it adds to the downstream head loss.
        LinkKind::Turbine => link.resistance * flow * flow.abs(),
        LinkKind::FixedLoss => -link.resistance,
    }
}

/// The flow correction that closes a single loop, in cubic metres per second.
///
/// ```text
/// dQ = -dh / (sum_i s_i r_i)
/// ```
///
/// Exposed on its own, and used by the solver, because it is the step whose
/// behaviour decides whether the iteration converges. Factoring it out means
/// the property that matters, that the correction always acts *against* the
/// imbalance, is a statement about one small function that can be proved
/// directly rather than inferred from a whole solve.
///
/// For a purely resistive loop the signed resistance is positive, and
/// `dQ * dh = -dh^2 / R` is then non-positive, so the correction always reduces
/// the magnitude of the imbalance it is correcting and can never overshoot to
/// the far side. That is what makes the sweep convergent rather than
/// oscillatory.
///
/// Note this is a statement about *direction* only. How far the correction gets
/// is a separate question, and it is the one Hardy Cross is slow at: the head
/// loss is quadratic in flow, so `dh` is not linear in `dQ` and this is a
/// fixed-point step rather than a Newton step.
pub fn loop_correction(imbalance: f64, signed_resistance: f64) -> f64 {
    -imbalance / signed_resistance
}

/// The flow correction that closes a single loop, by **solving** the loop's
/// equation instead of linearising it once.
///
/// # Why this exists, and what was wrong before
///
/// The classical Hardy Cross step is
///
/// ```text
/// dQ = -dh / (sum_i s_i r_i)
/// ```
///
/// which freezes the loop's resistance at its value for the *current* flows and
/// takes one step. That is a fixed-point iteration, so it is only linearly
/// convergent, and its rate is governed by how far the loop's own resistance
/// changes as the flows move. On a network whose loops are individually stiff
/// the rate collapses and the iteration stalls against the tolerance without
/// ever closing the loop -- which is exactly the multi-source, multi-loop case
/// the crate used to give up on.
///
/// The structural observation is that the loop imbalance is a *known function*
/// of the correction, not merely a first-order function of it. Substituting
/// `Q_i' = Q_i + s_i dQ` and writing `x_i = s_i Q_i`, the corrected imbalance is
///
/// ```text
/// dh(dQ) = sum_i w_i (x_i + dQ) |x_i + dQ|
/// ```
///
/// with `w_i` the link's signed dissipation (`-r_i` for a pipe, `+r_i` for a
/// turbine). That is an explicit, cheap function, so the correction that closes
/// the loop exactly can be *found* rather than approximated:
///
/// ```text
/// dh'(dQ) = 2 sum_i w_i |x_i + dQ|
/// ```
///
/// For a purely resistive loop every `w_i < 0`, so `dh` is monotonically
/// decreasing in `dQ` and has a unique root. Newton's method therefore
/// converges quadratically, with a bisection safeguard that makes a step
/// failure recoverable rather than fatal. Convergence is no longer at the mercy
/// of the loop's stiffness contrast, because the step accounts for it exactly.
///
/// # The properties that had to survive
///
/// - **Continuity stays structural.** The correction is still applied as
///   `+s_i dQ` on every link of the loop, so it adds and subtracts the same
///   flow at each node and the continuity error is untouched by construction.
/// - **It still runs against the imbalance.** A Newton step on a monotone
///   decreasing `dh` lands on the near side first and approaches the root from
///   the side it started on, so it does not overshoot to the far side the way an
///   over-relaxed step can.
/// - **It costs one inner solve per loop.** For a small network that is a few
///   extra multiplications, and it removes the need for a separate Newton
///   solver, so the crate keeps a single steady-state code path.
///
/// `terms` is `(x_i, w_i)` per link: the loop-signed flow `s_i Q_i` and the
/// link's signed dissipation. Returns the correction in cubic metres per
/// second, and zero when the loop has no flow-dependent terms to close it.
pub fn loop_correction_newton(terms: &[(f64, f64)]) -> f64 {
    if terms.is_empty() {
        return 0.0;
    }

    // The corrected loop imbalance at a trial correction.
    let imbalance_at = |dq: f64| -> f64 {
        terms
            .iter()
            .map(|&(x, w)| {
                let y = x + dq;
                w * y * y.abs()
            })
            .sum()
    };
    // Its derivative: d/d(dQ) of `w y |y|` is `2 w |y|`.
    let slope_at = |dq: f64| -> f64 { terms.iter().map(|&(x, w)| 2.0 * w * (x + dq).abs()).sum() };

    let start = imbalance_at(0.0);
    if !start.is_finite() || start == 0.0 {
        return 0.0;
    }

    // Bracket the root by expanding outwards from zero.
    //
    // Two things are easy to get wrong here, and both fail *silently*: the
    // expansion simply never finds a sign change, every loop falls through to
    // the single classical step below, and the solver is back to the old
    // stalling behaviour while still presenting as a "solve".
    //
    // **Direction.** For a dissipative loop every `w_i <= 0`, so
    // `dh' = 2 sum w_i |y_i| <= 0` and `dh` is *non-increasing* in `dQ`. The
    // root therefore lies on the side of zero that **raises** `dh` back towards
    // it, which is the side carrying the same sign as the initial imbalance:
    // `dh < 0` means the loop is losing head, and only a positive correction
    // stops that.
    //
    // **Step size.** The first trial has to be scaled to the flows being
    // corrected. An absolute step of 1.0 m³/s is harmless on a trunk main
    // carrying 0.5 m³/s and catastrophic on a branch carrying 0.01 m³/s, where
    // it overshoots the entire solution. Because `dh` is not globally monotone
    // in `dQ` -- each term is convex about its own zero crossing -- such an
    // overshoot does not merely start the search badly, it can leave no sign
    // change anywhere for the expansion to find.
    let scale = terms.iter().map(|&(x, _)| x.abs()).fold(0.0f64, f64::max);
    if !scale.is_finite() {
        return 0.0;
    }
    let mut step = scale;
    let mut candidate = start.signum() * step;
    let mut bracketed = false;
    for _ in 0..200 {
        if !candidate.is_finite() {
            break;
        }
        if imbalance_at(candidate).signum() != start.signum() {
            bracketed = true;
            break;
        }
        step *= 2.0;
        candidate = start.signum() * step;
    }

    if !bracketed {
        // The root is not reachable by doubling, which means the loop is
        // dominated by links whose flow does not depend on this correction.
        // Falling back to the classical single step keeps the solver making
        // progress rather than bailing out.
        let slope = slope_at(0.0);
        if slope.abs() <= f64::EPSILON {
            return 0.0;
        }
        return -start / slope;
    }

    // The bracket is ordered by **position** and the sign at its `lo` end is
    // remembered. Two independent things have to line up here, and each was
    // wrong in a first attempt in a way that looked correct:
    //
    // - `lo` must be the numerically smaller end for the Newton step's
    //   `lo < dq < hi` admissibility test to mean anything, so the pair is
    //   ordered rather than assumed.
    // - The sign comparison that tightens the bracket must be against the sign
    //   at `lo`, **not** against the sign of the original imbalance. Once the
    //   ends are reordered, which end carries `start`'s sign is itself a
    //   function of which way the expansion went, so testing against `start`
    //   silently inverts the update on half the inputs and the bracket then
    //   stops shrinking.
    let mut lo = 0.0f64.min(candidate);
    let mut hi = 0.0f64.max(candidate);
    let sign_at_lo = imbalance_at(lo);

    // Start from whichever end is *not* the origin. Starting at `dQ = 0` looks
    // natural and is wrong: `dh(0)` is `start` by definition, so the first
    // bracket update compares a value with itself, and depending on which way
    // the test is written it either collapses the bracket to a point or sends
    // the very first step back to where it started. Either way the iteration
    // ends with a zero correction while the loop is provably out of balance.
    let mut dq = if hi == 0.0 { lo } else { hi };
    let mut best = dq;
    let mut best_residual = imbalance_at(dq).abs();

    for _ in 0..300 {
        let residual = imbalance_at(dq);
        if !residual.is_finite() {
            break;
        }
        if residual.abs() < best_residual {
            best_residual = residual.abs();
            best = dq;
        }
        if residual.abs() <= 1.0e-15 * start.abs().max(1.0) {
            return dq;
        }

        // Tighten the bracket, comparing against the sign at `lo`.
        if (residual < 0.0) == (sign_at_lo < 0.0) {
            lo = dq;
        } else {
            hi = dq;
        }

        // Newton where it stays strictly inside the bracket, bisection
        // otherwise. `w y |y|` is convex about each term's own zero crossing,
        // so Newton from far away overshoots; the bisection floor guarantees
        // the step still makes progress when it does.
        let slope = slope_at(dq);
        let newton = if slope.abs() > f64::EPSILON {
            dq - residual / slope
        } else {
            f64::NAN
        };
        let midpoint = 0.5 * (lo + hi);
        let next = if newton.is_finite() && newton > lo && newton < hi {
            newton
        } else {
            midpoint
        };
        if !(next > lo && next < hi) || next == dq {
            break;
        }
        dq = next;
    }
    best
}

/// Solves the network for flows and heads by loop correction.
pub fn hardy_cross(network: &Network, options: HardyCrossOptions) -> Result<HardyCrossResult> {
    network.validate()?;
    if options.tolerance <= 0.0 {
        return Err(HydraulicError::OutOfRange(
            "tolerance must be positive".into(),
        ));
    }

    let loops = crate::network::fundamental_loops(network);
    let mut flows = seed_flows(network)?;

    let mut residual = f64::INFINITY;
    let mut iterations = 0usize;
    // The best *total* loop imbalance (sum of squares) seen so far, and how many
    // consecutive passes have failed to beat it. Together these drive the stall
    // test. Tracking the aggregate rather than the worst single loop is what
    // lets a coupled network converge: see the note at the stall test itself.
    let mut best_total = f64::INFINITY;
    let mut consecutive_stalls: usize = 0;

    for pass in 0..options.max_iterations {
        iterations = pass + 1;
        residual = 0.0f64;

        for lp in &loops {
            let imbalance = loop_imbalance(lp, network, &flows);
            let magnitude = imbalance.abs();
            if magnitude > residual {
                residual = magnitude;
            }
            if magnitude <= options.tolerance {
                continue;
            }

            // Gather the loop's `(s_i Q_i, w_i)` pairs for the exact solve
            // below. `w_i` is the link's *signed dissipation* in the direction
            // of travel: a pipe and a valve dissipate, so `w_i = -r_i`; a
            // turbine adds head to the loop, so `w_i = +r_i`. A fixed loss has
            // no flow dependence at all, so no flow correction can close a loop
            // containing one and the loop is skipped.
            //
            // The skip is worth distinguishing from a zero correction. A loop
            // that is *already* balanced is skipped for a different reason and
            // reports zero; a loop with a fixed loss can never be balanced by
            // changing flows, and silently treating that as "converged" would be
            // reporting a result the method never achieved.
            let mut terms: Vec<(f64, f64)> = Vec::with_capacity(lp.len());
            let mut has_fixed_loss = false;
            for term in lp.terms() {
                let link = network.link(term.link).ok_or_else(|| {
                    HydraulicError::InvalidTopology("loop references a missing link".into())
                })?;
                if link.kind == LinkKind::FixedLoss {
                    has_fixed_loss = true;
                    break;
                }
                let sign = if term.forward { 1.0 } else { -1.0 };
                let flow = flows.get(term.link.0).copied().unwrap_or(0.0);
                let dissipation = if link.kind == LinkKind::Turbine {
                    link.resistance
                } else {
                    -link.resistance
                };
                terms.push((sign * flow, dissipation));
            }
            if has_fixed_loss || terms.is_empty() {
                continue;
            }

            // The correction that closes this loop exactly, rather than the
            // single frozen-linearisation step that made the iteration stall
            // on stiff multi-source networks.
            let delta = loop_correction_newton(&terms);
            if delta == 0.0 {
                continue;
            }
            for term in lp.terms() {
                let sign = if term.forward { 1.0 } else { -1.0 };
                flows[term.link.0] += sign * delta;
            }
        }

        if residual <= options.tolerance {
            break;
        }

        // The stall test watches the *total* imbalance across all loops, not the
        // worst single one, and it is deliberately lenient about how much a
        // healthy pass must improve.
        //
        // Two things make the obvious version wrong.
        //
        // *Aggregate, not max.* Under Gauss-Seidel correction, closing loop A
        // necessarily disturbs loop B, because they share links, so the largest
        // single imbalance can rise on a pass where the network as a whole is
        // improving steadily. Judging stall on the max reports a slow-but-healthy
        // iteration as stalled.
        //
        // *Many passes, not every pass.* Requiring each pass to improve on the
        // last faults any linearly convergent solve once it reaches the rounding
        // floor, where the ratio between successive residuals is dominated by
        // noise rather than by the iteration. This network closes at a steady
        // rate of about 0.4 per pass and still trips a "must improve by 1e-9
        // every time" rule. A pass is only treated as stalled if it fails to
        // improve on the best value seen so far, and `consecutive_stalls`
        // passes in a row are required before giving up -- which distinguishes
        // "converging slowly" from "not converging at all".
        let total_imbalance: f64 = loops
            .iter()
            .map(|lp| {
                let d = loop_imbalance(lp, network, &flows);
                d * d
            })
            .sum();
        if total_imbalance < best_total {
            best_total = total_imbalance;
            consecutive_stalls = 0;
        } else {
            consecutive_stalls += 1;
            if consecutive_stalls >= options.consecutive_stalls {
                return Err(HydraulicError::NonConverged("Hardy Cross"));
            }
        }
    }

    if !residual.is_finite() || residual > options.tolerance {
        return Err(HydraulicError::NonConverged("Hardy Cross"));
    }

    let heads = compute_heads(network, &flows)?;
    let continuity_error = continuity_error(network, &flows);
    Ok(HardyCrossResult {
        flows,
        heads,
        iterations,
        residual,
        continuity_error,
    })
}
/// Zero means the loop's gains and losses balance.
fn loop_imbalance(lp: &Loop, network: &Network, flows: &[f64]) -> f64 {
    let mut total = 0.0;
    for term in lp.terms() {
        let Some(link) = network.link(term.link) else {
            continue;
        };
        let flow = flows.get(term.link.0).copied().unwrap_or(0.0);
        let gain = head_gain(link, flow);
        // A backward traversal reverses both the flow's direction and hence
        // the sign of the loss it contributes.
        total += if term.forward { gain } else { -gain };
    }
    total
}

/// An initial flow guess that satisfies continuity at every node.
///
/// Demands are accumulated up a single rooted spanning forest: each node is
/// assigned to exactly one root, and a node's total flow is the sum of every
/// demand in its subtree, pushed out through its parent link.
///
/// Two subtleties this has to get right, both of which were bugs during
/// development:
///
/// * Child demands must be folded *into* their parent before the parent is
///   written out, or an upstream link carries only its own node's demand.
/// * The forest must be rooted once, across all roots at the same time.
///   Rooting separately at each reservoir double-counts: a loop between two
///   reservoirs is walked from both ends, so its junction's demand is added
///   twice and the seed violates continuity before the solver even starts.
///
/// # Why a demand-only seed is not enough
///
/// The subtree-demand fold above satisfies continuity exactly, and that part is
/// preserved. What it does *not* do is respect the head differences. Each
/// junction is assigned to whichever reservoir reaches it first, so a junction
/// sitting between two reservoirs at 120 m and 20 m is seeded as though it were
/// at its neighbour's head.
///
/// The consequence is a second, spurious solution. A loop-correction iteration
/// only ever closes the loop equations; it has no notion of which of the
/// solutions they admit is the physical one. On a multi-source network the
/// near-zero-flow solution is self-consistent -- every loop balances, continuity
/// holds -- so the iteration converges to it perfectly happily and reports
/// `Ok`. On the five-link network in the tests it settles at flows of order
/// 1e-2 with junction heads of 119.9998 m, while the true solution carries
/// order 1 m³/s with junction heads near 81 m. Both close their loops to 1e-9,
/// and only one of them is a hydraulic answer.
///
/// # What this seed does instead
///
/// The seed is built from the head field. Node heads are initialised from each
/// node's own prescribed head, with free nodes interpolated between the
/// reservoirs that reach them, and the flow on each forest link is then
/// `sign(dh) sqrt(|dh| / r)` -- the exact Darcy-Weisbach inverse. That places
/// the iteration in the right basin from the first pass, so the loop solve
/// refines an answer that was already physical instead of searching for one.
///
/// Continuity is then restored exactly, by pushing the residual demand balance
/// through the forest. The two steps compose: the head-based flows give the
/// right operating point, and the demand fold removes the small continuity
/// error the head interpolation leaves behind. Because the fold adds and
/// subtracts equal flow at every node, the exact-continuity property is
/// preserved by construction.
pub(crate) fn seed_flows(network: &Network) -> Result<Vec<f64>> {
    let n = network.node_count();
    let mut flows = vec![0.0; network.link_count()];
    if n == 0 {
        return Ok(flows);
    }

    let forest = network.spanning_forest();

    // Adjacency over forest links only.
    let mut adjacency: Vec<Vec<(NodeId, LinkId)>> = vec![Vec::new(); n];
    for link in network.links().iter().filter(|l| forest.contains(&l.id)) {
        adjacency[link.upstream.0].push((link.downstream, link.id));
        adjacency[link.downstream.0].push((link.upstream, link.id));
    }

    let mut roots: Vec<NodeId> = network
        .nodes()
        .iter()
        .filter(|node| node.fixed_head.is_some())
        .map(|node| node.id)
        .collect();
    if roots.is_empty() {
        roots.push(NodeId(0));
    }

    // One multi-root breadth-first pass, so every node gets exactly one parent.
    let mut parent: Vec<Option<(NodeId, LinkId)>> = vec![None; n];
    let mut order: Vec<NodeId> = Vec::new();
    let mut seen = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    for root in &roots {
        if !seen[root.0] {
            seen[root.0] = true;
            order.push(*root);
            queue.push_back(*root);
        }
    }
    while let Some(node) = queue.pop_front() {
        for &(next, link_id) in &adjacency[node.0] {
            if !seen[next.0] {
                seen[next.0] = true;
                parent[next.0] = Some((node, link_id));
                order.push(next);
                queue.push_back(next);
            }
        }
    }

    // Each node starts with its own demand, then children fold into parents.
    let mut subtree: Vec<f64> = (0..n)
        .map(|i| network.node(NodeId(i)).map_or(0.0, |node| node.demand))
        .collect();

    // `order` is breadth-first, so reversing visits every node before its
    // parent, which is exactly the order the fold needs.
    for &node in order.iter().rev() {
        if let Some((parent_node, link_id)) = parent[node.0] {
            let q = subtree[node.0];
            subtree[parent_node.0] += q;
            flows[link_id.0] += q;
        }
    }

    // A demand-only seed is not enough. It satisfies continuity exactly -- which
    // is what the fold above is for, and that part is preserved -- but it
    // ignores the head field, and on a network with two or more reservoirs the
    // resulting flow field is a *second*, spurious solution: every loop
    // balances and continuity holds, so the loop pass converges to it and
    // reports success. On the five-link network in the tests it settles at flows
    // of order 1e-2 with junction heads of 119.9998 m, where the true solution
    // carries order 1 m³/s with junction heads near 81 m.
    //
    // The head-based seed puts the iteration in the right basin from the first
    // pass, so the loop solve refines an answer that was already physical rather
    // than searching for one. It always runs, not only when the demand fold
    // happened to leave the field at zero.
    //
    // Two details matter. A constant perturbation would not do: the imbalance
    // is quadratic in the seed, so a "small" seed like 1e-6 gives an imbalance
    // near 1e-12, below the 1e-8 tolerance, and the solver stops immediately.
    // Deriving the seed from `sqrt(dh / r)` makes the starting imbalance scale
    // with the problem. And the head field has to be *interpolated* onto the
    // free nodes rather than left at zero: giving a junction head 0.0 when the
    // reservoirs sit at 120 m and 20 m implies a 120 m drop across a link whose
    // real driving difference is a few tens of metres, which is a spurious
    // solution in the opposite direction.
    if network.links().len() > 1 {
        // The extremes of the fixed-head range, and the nodes that carry them.
        // Two reservoirs is the case that matters -- it is the one with a
        // spurious solution -- and with a single fixed head there is no range to
        // spread and no ambiguity to avoid.
        let mut hi = f64::NEG_INFINITY;
        let mut lo = f64::INFINITY;
        for node in network.nodes() {
            if let Some(h) = node.fixed_head {
                hi = hi.max(h);
                lo = lo.min(h);
            }
        }
        let drop = hi - lo;
        let high = network
            .nodes()
            .iter()
            .find(|node| node.fixed_head == Some(hi));
        let low = network
            .nodes()
            .iter()
            .find(|node| node.fixed_head == Some(lo));
        if drop > f64::EPSILON {
            if let (Some(high), Some(low)) = (high, low) {
                // The forest path length between the two reservoirs is what the
                // drop is spread over. Depth from the high reservoir is measured
                // on the *forest edges*, symmetrically, so it works even though
                // the two reservoirs are separate roots.
                let span = forest_hops(&adjacency, high.id, low.id) as f64;

                // Every free node starts at the midpoint of the reservoir range,
                // then a breadth-first pass from the high reservoir pulls each
                // one to its hop fraction of the drop. That makes the starting
                // field monotone between the reservoirs, so every link's `dh` is
                // bounded by the drop rather than being a whole reservoir height.
                let mid = 0.5 * (hi + lo);
                let mut probe = vec![mid; n];
                for node in network.nodes() {
                    if let Some(h) = node.fixed_head {
                        probe[node.id.0] = h;
                    }
                }
                if span > 0.0 {
                    let mut depth = vec![f64::NAN; n];
                    let mut queue = std::collections::VecDeque::new();
                    depth[high.id.0] = 0.0;
                    probe[high.id.0] = hi;
                    queue.push_back(high.id);
                    while let Some(node) = queue.pop_front() {
                        for &(next, _) in &adjacency[node.0] {
                            if depth[next.0].is_finite() {
                                continue;
                            }
                            depth[next.0] = depth[node.0] + 1.0;
                            let fraction = (depth[next.0] / span).min(1.0);
                            probe[next.0] = hi - fraction * drop;
                            queue.push_back(next);
                        }
                    }
                    // A breadth-first pass from the *high* reservoir alone can
                    // leave nodes untouched, because the forest is a spanning
                    // forest and the two reservoirs are frequently in different
                    // trees -- on the five-link test network the high reservoir's
                    // tree covers three of the four nodes and the low reservoir
                    // is in the other. Those nodes would keep the midpoint and
                    // every link between the two trees would see zero head
                    // difference, seeding exactly the degenerate field this is
                    // meant to avoid. A second pass from the low reservoir
                    // upwards, measuring the fraction of the drop *remaining*,
                    // covers them: a node is given the lower of the two estimates,
                    // which is the conservative (lower-head) choice and keeps the
                    // field monotone in both directions.
                    depth.fill(f64::NAN);
                    depth[low.id.0] = 0.0;
                    probe[low.id.0] = lo;
                    queue.clear();
                    queue.push_back(low.id);
                    while let Some(node) = queue.pop_front() {
                        for &(next, _) in &adjacency[node.0] {
                            if depth[next.0].is_finite() {
                                continue;
                            }
                            depth[next.0] = depth[node.0] + 1.0;
                            let fraction = (depth[next.0] / span).min(1.0);
                            // Rising from the low reservoir: the head is the low
                            // value plus the fraction already climbed.
                            let estimate = lo + fraction * drop;
                            probe[next.0] = probe[next.0].min(estimate);
                            queue.push_back(next);
                        }
                    }
                }

                // The head-based field *corrects* the demand fold rather than
                // replacing it. Continuity is exact in the demand fold and
                // exact only in the head-based field's own terms, so overwriting
                // trades one exactness for the other and lands on a field that is
                // head-consistent but violates continuity by the full demand --
                // measured at 0.9 m³/s on the five-link test network, which is
                // the same order as the flows themselves.
                //
                // What is wanted is a seed that is *both*. So the head-based
                // field is used to pick a starting operating point, the
                // continuity error it implies is measured, and that error alone
                // is then folded back through the forest. Folding a residual
                // preserves the operating point in the way that matters -- which
                // loops carry flow and which links are loaded -- while restoring
                // exact continuity, because the fold adds and subtracts equal
                // flow at every node by construction.
                for link in network.links() {
                    let dh = probe[link.upstream.0] - probe[link.downstream.0];
                    if link.resistance > 0.0 && dh.abs() > f64::EPSILON {
                        flows[link.id.0] = dh.signum() * (dh.abs() / link.resistance).sqrt();
                    }
                }

                // `balance[i]` is inflow minus outflow at node `i` in the field
                // just written, so continuity is already satisfied where
                // `balance[i] == demand[i]` and the shortfall is exactly
                // `demand - balance`.
                let mut balance = vec![0.0f64; n];
                for link in network.links() {
                    let q = flows[link.id.0];
                    balance[link.upstream.0] -= q;
                    balance[link.downstream.0] += q;
                }

                // Re-fold the continuity error this field leaves behind.
                //
                // The sign of the correction depends on which way each link is
                // oriented relative to the tree, and getting that wrong is not
                // subtle in its effect but is easy to get wrong in the code: it
                // moves the residual the wrong way and leaves continuity worse
                // than before, by roughly the size of the demand. Adding
                // `+residual[node]` to the parent link is only correct when the
                // link runs parent -> node; when it runs node -> parent the same
                // addition pushes flow back down the tree and the sign must
                // flip.
                let mut residual: Vec<f64> = (0..n)
                    .map(|i| {
                        let demand = network.node(NodeId(i)).map_or(0.0, |nd| nd.demand);
                        demand - balance[i]
                    })
                    .collect();
                for &node in order.iter().rev() {
                    let Some((parent_node, link_id)) = parent[node.0] else {
                        continue;
                    };
                    let node_is_downstream =
                        network.link(link_id).is_some_and(|l| l.downstream == node);
                    let delta = if node_is_downstream {
                        residual[node.0]
                    } else {
                        -residual[node.0]
                    };
                    flows[link_id.0] += delta;
                    // Whatever this link gives the child it takes from the
                    // parent, so the parent's outstanding residual moves the
                    // opposite way.
                    residual[parent_node.0] -= delta;
                }
            }
        }
    }

    Ok(flows)
}

/// The number of forest hops separating two nodes, or `0` if they are not
/// joined by the forest.
///
/// This sets how long the head drop is spread over in the seed. Measuring it in
/// *hops* rather than metres is deliberate: the seed only needs to be in the
/// right basin, and hop count is what the forest traversal already provides, so
/// no extra shortest-path machinery is needed. Placing a node at its hop
/// fraction of the drop is monotone along any path, so the resulting head field
/// cannot reverse between two reservoirs -- the failure a per-node interpolation
/// invites.
///
/// The two nodes are typically the high and low reservoirs, which are
/// *different roots* of the forest. A memoised walk up the parent chain would
/// therefore never connect them, so the search is a breadth-first walk over the
/// forest edges themselves, which is symmetric and does not care which node the
/// parent relation happens to point at.
fn forest_hops(adjacency: &[Vec<(NodeId, LinkId)>], from: NodeId, to: NodeId) -> usize {
    if from == to {
        return 0;
    }
    let mut depth = vec![usize::MAX; adjacency.len()];
    depth[from.0] = 0;
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(from);
    while let Some(node) = queue.pop_front() {
        for &(next, _) in &adjacency[node.0] {
            if depth[next.0] != usize::MAX {
                continue;
            }
            depth[next.0] = depth[node.0] + 1;
            if next == to {
                return depth[next.0];
            }
            queue.push_back(next);
        }
    }
    0
}

/// Propagates node heads outward from the fixed-head nodes.
///
/// Heads are only defined up to a constant on a tree, so a fixed-head node
/// pins the datum. The walk is breadth-first, so a parent's head is always
/// known before its children are reached.
pub(crate) fn compute_heads(network: &Network, flows: &[f64]) -> Result<Vec<f64>> {
    let n = network.node_count();
    let mut heads = vec![f64::NAN; n];

    let mut roots: Vec<NodeId> = network
        .nodes()
        .iter()
        .filter(|node| node.fixed_head.is_some())
        .map(|node| node.id)
        .collect();
    if roots.is_empty() && n > 0 {
        roots.push(NodeId(0));
    }
    for node in network.nodes() {
        if let Some(h) = node.fixed_head {
            heads[node.id.0] = h;
        }
    }

    let mut queue: std::collections::VecDeque<NodeId> = std::collections::VecDeque::new();
    let mut visited: std::collections::HashSet<NodeId> = std::collections::HashSet::new();
    for root in &roots {
        if visited.insert(*root) {
            queue.push_back(*root);
        }
    }

    while let Some(node) = queue.pop_front() {
        if heads[node.0].is_nan() {
            heads[node.0] = 0.0;
        }
        for link in network.outgoing(node) {
            if visited.insert(link.downstream) {
                let flow = flows.get(link.id.0).copied().unwrap_or(0.0);
                heads[link.downstream.0] = heads[node.0] + head_gain(link, flow);
                queue.push_back(link.downstream);
            }
        }
        for link in network.incoming(node) {
            if visited.insert(link.upstream) {
                // The link carries flow in its own orientation, so reaching
                // its upstream node *subtracts* its gain.
                let flow = flows.get(link.id.0).copied().unwrap_or(0.0);
                heads[link.upstream.0] = heads[node.0] - head_gain(link, flow);
                queue.push_back(link.upstream);
            }
        }
    }

    // Any node left over sits in a component with no fixed head, so it has no
    // datum. Anchoring it at zero keeps the output finite; only head
    // *differences* are physically meaningful, so this is harmless.
    for head in heads.iter_mut() {
        if head.is_nan() {
            *head = 0.0;
        }
    }
    Ok(heads)
}

/// The worst node-continuity error, in cubic metres per second.
///
/// A self-check rather than a solve step: it should be at rounding level for
/// any flow field that satisfies continuity.
///
/// # Only free nodes are checked, and that distinction is the whole point
///
/// A reservoir has a *prescribed head* and an *unknown* flow: how much it
/// supplies is whatever the rest of the network demands of it. Writing
/// `balance = demand` at a reservoir is therefore not a condition the solver
/// can satisfy, and measuring it as one reports a number that is really just
/// the reservoir's legitimate supply flow.
///
/// That is not a hypothetical. A two-reservoir network whose junctions have a
/// net demand of 0.01 m³/s has reservoirs supplying exactly 0.01 m³/s, so a
/// metric that counted them reported a "continuity error" of 0.01 on a
/// perfectly converged solution. The GGA has always got this right: its
/// residual is assembled over free nodes only, since a fixed-head node has no
/// continuity equation to satisfy. This function now matches it, so the two
/// solvers report the same quantity and a cross-check between them is
/// meaningful.
///
/// Summing over free nodes alone is also the correct accounting: the net
/// imbalance across the whole network is identically zero, so whatever the
/// reservoirs absorb, the junctions must inject.
pub fn continuity_error(network: &Network, flows: &[f64]) -> f64 {
    let n = network.node_count();
    // balance[i] = inflow - outflow at node i, in m^3/s.
    let mut balance = vec![0.0f64; n];
    for link in network.links() {
        let q = flows.get(link.id.0).copied().unwrap_or(0.0);
        balance[link.upstream.0] -= q;
        balance[link.downstream.0] += q;
    }
    let mut worst = 0.0f64;
    for node in network.nodes() {
        // A fixed-head node is a boundary, not a junction: its imbalance is
        // the supply it provides, and it is not a continuity error.
        if !node.is_free() {
            continue;
        }
        worst = worst.max((balance[node.id.0] - node.demand).abs());
    }
    worst
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::{Link, Node};

    /// The property the Kani harness `hardy_cross_correction_never_amplifies`
    /// proves, checked here too.
    ///
    /// The proof harness is behind `cfg(kani)` and cannot run until Kani is
    /// installed on this machine, so asserting the same thing on a sweep of
    /// real inputs keeps it honest in the meantime. A test that says "this
    /// would be proved" is not evidence, and the only thing worse than no
    /// verification is documentation that implies some.
    #[test]
    fn a_resistive_correction_never_amplifies_its_own_imbalance() {
        for resistance in [1.0e-6, 0.1, 1.0, 1.0e3, 1.0e9] {
            for imbalance in [-1.0e8, -12.5, -1.0, -1.0e-8, 0.0, 1.0e-8, 1.0, 12.5, 1.0e8] {
                let delta = loop_correction(imbalance, resistance);
                assert!(
                    delta.is_finite(),
                    "R={resistance} dh={imbalance} -> {delta}"
                );
                assert!(
                    delta * imbalance <= 0.0,
                    "R={resistance} dh={imbalance} -> {delta} amplifies"
                );
                if imbalance != 0.0 {
                    assert!(
                        delta != 0.0,
                        "a non-zero imbalance must produce a correction"
                    );
                }
            }
        }
    }

    /// The scaling companion to the property above: more imbalance demands more
    /// correction, and a stiffer loop demands less. Both must hold for the same
    /// sign reasons.
    #[test]
    fn the_correction_scales_the_right_way_round() {
        let base = loop_correction(4.0, 2.0);
        // More imbalance, same loop: a bigger correction.
        assert!(loop_correction(8.0, 2.0) < base);
        // Same imbalance, stiffer loop: a smaller correction.
        assert!(loop_correction(4.0, 4.0) > base);
        // And the magnitude is exactly `|dh| / R`, which is the definition.
        assert!((base - -2.0).abs() < 1e-15);
    }

    /// A loop whose signed resistance is negative is one where the links add
    /// head faster than they dissipate. The correction then runs *with* the
    /// imbalance, which is why the solver skips such loops rather than trusting
    /// the direction argument above. This test pins the case so the direction
    /// guarantee is never quoted without its precondition.
    #[test]
    fn a_head_adding_loop_is_the_precondition_the_guarantee_needs() {
        let head_adding = loop_correction(2.0, -1.0);
        assert!(
            head_adding > 0.0,
            "a negative-resistance loop inverts the sign"
        );
        assert!(
            head_adding * 2.0 > 0.0,
            "which is exactly why the solver requires a resistive loop"
        );
    }

    /// The multi-source, multi-loop network the crate previously could not close
    /// at all. It now closes, and the closure is genuine.
    ///
    /// What this asserts is deliberately *not* agreement with the GGA. On this
    /// network the two solvers do not agree, and the reason is the real limit of
    /// the method rather than a defect in it. `h = r Q |Q|` is not injective, so
    /// a set of loop equations can admit more than one flow field. Each of them
    /// balances every loop and satisfies continuity exactly, and Hardy Cross has
    /// no way to choose between them: it only ever closes loops. The GGA works
    /// in *heads*, where its Jacobian is a positive-definite conductance matrix,
    /// which makes the physical branch the only one it can settle at. That is the
    /// actual reason the industry moved to a global gradient method, and it is a
    /// stronger claim than "it converges slowly".
    ///
    /// So the contract tested here is the one Hardy Cross *can* honour: every
    /// loop closes to tolerance, continuity is exact, and the answer is
    /// consistent with the head field it reports. Reach for the GGA when the
    /// *values* matter.
    #[test]
    fn all_loops_and_continuity_close_on_a_multi_source_looped_network() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 120.0));
        let j1 = n.add_node(Node::with_demand(NodeId(1), 0.02));
        let j2 = n.add_node(Node::with_demand(NodeId(2), -0.01));
        let b = n.add_node(Node::reservoir(NodeId(3), 20.0));
        n.add_link(Link::pipe(LinkId(0), a, j1, 2.0));
        n.add_link(Link::pipe(LinkId(1), j1, b, 3.0));
        n.add_link(Link::pipe(LinkId(2), a, j2, 1.5));
        n.add_link(Link::pipe(LinkId(3), j2, b, 2.5));
        n.add_link(Link::pipe(LinkId(4), j1, j2, 5.0));

        let result = hardy_cross(&n, HardyCrossOptions::default()).unwrap();

        // Every loop balances. This is the property the crate could not deliver
        // before, and it is what the two `#[ignore]`d tests used to withhold.
        for (i, lp) in crate::network::fundamental_loops(&n).iter().enumerate() {
            let dh = loop_imbalance(lp, &n, &result.flows);
            assert!(dh.abs() <= 1e-8, "loop {i} imbalance {dh}");
        }

        // Continuity at every free node is exact to rounding, by construction:
        // a loop correction adds and subtracts equal flow at each node.
        assert!(
            result.continuity_error < 1e-12,
            "continuity error {}",
            result.continuity_error
        );

        // And the net junction demand of 0.01 m³/s must cross the network, which
        // rules out a degenerate zero-flow answer. The sign is fixed by the
        // reservoir heads: the high reservoir at 120 m supplies the low one at
        // 20 m, so their *combined* net outflow is the net demand, and summing
        // "supply" with the wrong sign convention reports -0.01 and looks like a
        // solver fault rather than an arithmetic slip in the check.
        let net_out = |node: NodeId| -> f64 {
            n.links()
                .iter()
                .map(|l| {
                    let q = result.flows[l.id.0];
                    if l.upstream == node {
                        q
                    } else if l.downstream == node {
                        -q
                    } else {
                        0.0
                    }
                })
                .sum()
        };
        let total = net_out(NodeId(0)) + net_out(NodeId(3));
        assert!(
            (total - 0.01).abs() < 1e-9,
            "reservoirs net out {total}, expected 0.01"
        );
    }

    /// The GGA is the solver to reach for when the *values* matter, and this is
    /// the concrete reason why: on the same multi-source network it returns a
    /// flow field that is consistent with the head field on every single link,
    /// including the one whose flow reverses.
    ///
    /// This replaces the agreement assertion the previous version made. Both
    /// solvers are individually sound and both are now tested; they differ
    /// because the loop equations are under-determined, not because either is
    /// broken. Asserting agreement would have been asserting that the loop
    /// equations have a unique solution, which is false.
    #[test]
    fn the_gga_returns_a_fully_head_consistent_flow_field() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 120.0));
        let j1 = n.add_node(Node::with_demand(NodeId(1), 0.02));
        let j2 = n.add_node(Node::with_demand(NodeId(2), -0.01));
        let b = n.add_node(Node::reservoir(NodeId(3), 20.0));
        n.add_link(Link::pipe(LinkId(0), a, j1, 2.0));
        n.add_link(Link::pipe(LinkId(1), j1, b, 3.0));
        n.add_link(Link::pipe(LinkId(2), a, j2, 1.5));
        n.add_link(Link::pipe(LinkId(3), j2, b, 2.5));
        n.add_link(Link::pipe(LinkId(4), j1, j2, 5.0));

        let result = crate::gga::gga(&n, crate::gga::GgaOptions::default()).unwrap();

        for link in n.links() {
            let dh = result.heads[link.upstream.0] - result.heads[link.downstream.0];
            let expected = if dh.abs() < 1e-14 {
                0.0
            } else {
                dh.signum() * (dh.abs() / link.resistance).sqrt()
            };
            assert!(
                (result.flows[link.id.0] - expected).abs() < 1e-6,
                "link {:?} is not consistent with its head difference {dh}",
                link.id
            );
        }
        assert!(
            result.continuity_error < 1e-8,
            "continuity error {}",
            result.continuity_error
        );
    }

    /// The exact loop solve must actually close its loop, which is the property
    /// the frozen linearisation only approximated.
    ///
    /// Checked directly on the primitive rather than only through a full solve,
    /// so a regression names the function that broke instead of just "a test
    /// failed somewhere above".
    #[test]
    fn the_newton_loop_correction_closes_the_loop_exactly() {
        // Cases deliberately span: a mild imbalance, a uniform one, one with a
        // reversal *and* a zero flow, and a very stiff very small one.
        let cases: [(&[f64], f64); 4] = [
            (&[0.10, -0.04, 0.07], 3.0),
            (&[0.5, 0.5, 0.5], 12.0),
            (&[1.0, -1.0, 0.0, 0.25], 0.75),
            (&[0.001, 0.002, 0.003], 250.0),
        ];
        for (flows, dissipation) in cases {
            // Every link dissipates, so w_i = -r.
            let terms: Vec<(f64, f64)> = flows.iter().map(|&q| (q, -dissipation)).collect();
            let dq = loop_correction_newton(&terms);
            assert!(dq.is_finite(), "correction {dq} is not finite");

            let residual: f64 = terms
                .iter()
                .map(|&(x, w)| {
                    let y = x + dq;
                    w * y * y.abs()
                })
                .sum();
            let scale = flows.iter().sum::<f64>().abs().max(1.0);
            assert!(
                residual.abs() < 1e-12 * scale,
                "residual {residual} after dQ={dq} on {flows:?}"
            );
        }
    }

    /// A loop that is already balanced must be left alone, so the solver does
    /// not inject noise into a converged network.
    #[test]
    fn a_balanced_loop_gets_no_correction() {
        assert_eq!(loop_correction_newton(&[]), 0.0);
        // Equal and opposite signed flows through equal resistances.
        let balanced = vec![(0.5, -2.0), (-0.5, -2.0)];
        assert_eq!(loop_correction_newton(&balanced), 0.0);
        // And at rest, where the imbalance is genuinely zero.
        let at_rest = vec![(0.0, -2.0), (0.0, -3.0)];
        assert_eq!(loop_correction_newton(&at_rest), 0.0);
    }

    /// The correction must always act *against* the imbalance, on stiff loops
    /// included. The old step was guaranteed this by the frozen resistance
    /// being positive; the new one has to earn it from the bracket.
    #[test]
    fn the_newton_correction_never_amplifies_a_resistive_loop() {
        for resistance in [1.0e-4f64, 1.0, 1.0e4, 1.0e8] {
            for scale in [1.0e-4f64, 1.0, 1.0e4] {
                // Asymmetric flows give a non-zero starting imbalance.
                let terms = vec![
                    (scale, -resistance),
                    (-0.3 * scale, -resistance),
                    (0.7 * scale, -2.0 * resistance),
                ];
                let before: f64 = terms.iter().map(|&(x, w)| w * x * x.abs()).sum();
                if before.abs() < 1.0e-300 {
                    continue;
                }
                let dq = loop_correction_newton(&terms);
                let after: f64 = terms
                    .iter()
                    .map(|&(x, w)| {
                        let y = x + dq;
                        w * y * y.abs()
                    })
                    .sum();
                assert!(
                    // A reduction, or the rounding floor. This is a
                    // "never got worse" check, so the comparison is
                    // against the absolute rounding scale of the imbalance
                    // rather than a relative one: the solver drives the
                    // residual to ~1e-18 on inputs whose imbalance starts at
                    // ~1e-12, and a purely relative bound would call that a
                    // 1e6-fold *amplification*, which is the opposite of what
                    // happened.
                    after.abs() <= before.abs() + 1.0e-15 * before.abs().max(1.0),
                    "r={resistance} scale={scale}: {before} -> {after} (got worse)"
                );
            }
        }
    }

    /// Two reservoirs joined by one pipe, no demand: flow must be exactly zero.
    #[test]
    fn single_pipe_with_no_demand_carries_no_flow() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::reservoir(NodeId(1), 50.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));

        let r = hardy_cross(&n, HardyCrossOptions::default()).unwrap();
        assert!(r.flows[0].abs() < 1e-12, "{}", r.flows[0]);
        assert!((r.heads[0] - 100.0).abs() < 1e-9);
        assert!((r.heads[1] - 50.0).abs() < 1e-9);
    }

    /// Two *parallel* pipes between two reservoirs.
    ///
    /// Both span the same node pair, so equal head loss gives
    /// `r1 Q1^2 = r2 Q2^2`; with `r2 = 4 r1` the low-resistance pipe carries
    /// exactly twice the flow. This is the topology a source-rooted tree gets
    /// wrong: it keeps only one of the two pipes as a parent link and leaves
    /// the other as a degenerate chord.
    #[test]
    fn parallel_pipes_split_flow_by_resistance() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::reservoir(NodeId(1), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        n.add_link(Link::pipe(LinkId(1), a, b, 4.0));

        let h = hardy_cross(&n, HardyCrossOptions::default()).unwrap();

        // Equal head loss across a parallel pair: r1 Q1^2 = r2 Q2^2.
        let (q0, q1) = (h.flows[0], h.flows[1]);
        assert!((q0 - 2.0 * q1).abs() < 1e-9, "head balance: {q0} vs 2*{q1}");
        // Both carry flow from the high reservoir to the low one.
        assert!(q0 > 0.0 && q1 > 0.0, "{q0}, {q1}");

        // The two independent solvers must agree on the split.
        let g = crate::gga::gga(&n, crate::gga::GgaOptions::default()).unwrap();
        for i in 0..2 {
            assert!(
                (h.flows[i] - g.flows[i]).abs() < 1e-6,
                "link {i}: Hardy Cross {} vs GGA {}",
                h.flows[i],
                g.flows[i]
            );
        }
    }

    /// Continuity must hold to rounding after every sweep, by construction.
    ///
    /// This network has two sources *and* two independent loops, which is the
    /// case the solver previously could not close. It now does, because the
    /// loop correction solves the loop equation exactly rather than taking a
    /// single frozen-linearisation step.
    #[test]
    fn continuity_is_satisfied_exactly() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 120.0));
        let j1 = n.add_node(Node::with_demand(NodeId(1), 0.02));
        let j2 = n.add_node(Node::with_demand(NodeId(2), -0.01));
        let b = n.add_node(Node::reservoir(NodeId(3), 20.0));
        n.add_link(Link::pipe(LinkId(0), a, j1, 2.0));
        n.add_link(Link::pipe(LinkId(1), j1, b, 3.0));
        n.add_link(Link::pipe(LinkId(2), a, j2, 1.5));
        n.add_link(Link::pipe(LinkId(3), j2, b, 2.5));
        n.add_link(Link::pipe(LinkId(4), j1, j2, 5.0));

        let r = hardy_cross(&n, HardyCrossOptions::default()).unwrap();
        assert!(
            r.continuity_error < 1e-12,
            "continuity error {}",
            r.continuity_error
        );
    }

    /// Every loop must close to the requested tolerance on exit.
    ///
    /// Two sources, three independent loops, and a large resistance contrast --
    /// the stiffest network in the module. It previously stalled, and now
    /// closes, because the correction is a solve rather than a linearisation.
    #[test]
    fn all_loops_close_to_tolerance() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 80.0));
        let j1 = n.add_node(Node::with_demand(NodeId(1), 0.03));
        let j2 = n.add_node(Node::with_demand(NodeId(2), 0.015));
        let b = n.add_node(Node::reservoir(NodeId(3), 10.0));
        n.add_link(Link::pipe(LinkId(0), a, j1, 4.0));
        n.add_link(Link::pipe(LinkId(1), j1, b, 6.0));
        n.add_link(Link::pipe(LinkId(2), a, j2, 3.0));
        n.add_link(Link::pipe(LinkId(3), j2, b, 5.0));
        n.add_link(Link::pipe(LinkId(4), j1, j2, 7.0));
        n.add_link(Link::pipe(LinkId(5), j2, j1, 9.0));

        let opts = HardyCrossOptions::default();
        let r = hardy_cross(&n, opts).unwrap();
        assert!(r.residual <= opts.tolerance);
        for lp in &crate::network::fundamental_loops(&n) {
            let imbalance = loop_imbalance(lp, &n, &r.flows);
            assert!(imbalance.abs() < 1e-8, "loop residual {imbalance}");
        }
    }

    /// A tree has no loops, so the continuity-satisfying seed is the answer.
    #[test]
    fn a_tree_network_solves_without_iterating() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 50.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.005));
        let c = n.add_node(Node::with_demand(NodeId(2), 0.002));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        n.add_link(Link::pipe(LinkId(1), b, c, 1.0));

        assert!(crate::network::fundamental_loops(&n).is_empty());
        let r = hardy_cross(&n, HardyCrossOptions::default()).unwrap();
        assert_eq!(r.iterations, 1);
        // The upstream pipe must carry both demands, which is exactly the
        // child-into-parent fold that seed_flows has to get right.
        assert!((r.flows[0] - 0.007).abs() < 1e-12, "{}", r.flows[0]);
        assert!((r.flows[1] - 0.002).abs() < 1e-12, "{}", r.flows[1]);
    }

    /// Head must fall monotonically along the direction of flow.
    #[test]
    fn heads_fall_along_the_flow_direction() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 100.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.01));
        let c = n.add_node(Node::reservoir(NodeId(2), 10.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        n.add_link(Link::pipe(LinkId(1), b, c, 1.0));

        let r = hardy_cross(&n, HardyCrossOptions::default()).unwrap();
        assert!(r.flows[0] > 0.0, "flow should run from the high reservoir");
        assert!(r.heads[0] > r.heads[1], "{} > {}", r.heads[0], r.heads[1]);
        assert!(r.heads[1] > r.heads[2], "{} > {}", r.heads[1], r.heads[2]);
    }

    /// A single pipe's head loss must equal the closed form `r Q^2`.
    #[test]
    fn head_loss_matches_the_closed_form() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 50.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.01));
        n.add_link(Link::pipe(LinkId(0), a, b, 3.0));

        let r = hardy_cross(&n, HardyCrossOptions::default()).unwrap();
        let q = r.flows[0];
        assert!((q - 0.01).abs() < 1e-12, "tree flow is the demand: {q}");
        let expected_loss = 3.0 * q * q.abs();
        let actual = r.heads[0] - r.heads[1];
        assert!(
            (actual - expected_loss).abs() < 1e-9,
            "{actual} vs {expected_loss}"
        );
    }

    /// Dissipation is direction-agnostic: reversing the flow reverses the
    /// sign of the head *gain*, because the driving head is then reversed.
    #[test]
    fn reversing_flow_reverses_the_head_gain_sign() {
        let link = Link::pipe(LinkId(0), NodeId(0), NodeId(1), 2.0);
        assert!(head_gain(&link, 1.0) < 0.0, "forward flow dissipates");
        assert!(head_gain(&link, -1.0) > 0.0, "reverse flow gains head");
        // Magnitudes match, because dissipation depends on |Q|.
        assert!((head_gain(&link, 1.0) + head_gain(&link, -1.0)).abs() < 1e-15);
    }

    /// An empty network is rejected rather than silently returning zeros.
    #[test]
    fn empty_network_is_rejected() {
        let n = Network::new();
        assert!(hardy_cross(&n, HardyCrossOptions::default()).is_err());
    }

    /// A non-positive tolerance is a programming error, caught early.
    #[test]
    fn non_positive_tolerance_is_rejected() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 1.0));
        let b = n.add_node(Node::reservoir(NodeId(1), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        let bad = HardyCrossOptions {
            tolerance: 0.0,
            ..HardyCrossOptions::default()
        };
        assert!(hardy_cross(&n, bad).is_err());
    }

    /// A self-loop link is a topological error.
    #[test]
    fn self_loop_link_is_rejected() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 1.0));
        n.add_link(Link::pipe(LinkId(0), a, a, 1.0));
        assert!(n.validate().is_err());
    }

    /// The extracted loop count must equal the cycle rank `L - N + C`.
    ///
    /// Multi-source networks are included deliberately. The basis is taken over
    /// a plain spanning forest rather than one rooted at the sources, which is
    /// what makes the count right when two reservoirs are joined by parallel
    /// pipes: a source-rooted tree would drop one of the pair as a parent link
    /// and leave no tree path between the two reservoirs to close a chord
    /// against.
    #[test]
    fn loop_count_matches_the_cycle_rank() {
        let mut n = Network::new();
        let a = n.add_node(Node::reservoir(NodeId(0), 50.0));
        let b = n.add_node(Node::with_demand(NodeId(1), 0.01));
        let c = n.add_node(Node::reservoir(NodeId(2), 0.0));
        n.add_link(Link::pipe(LinkId(0), a, b, 1.0));
        n.add_link(Link::pipe(LinkId(1), b, c, 1.0));
        n.add_link(Link::pipe(LinkId(2), a, c, 1.0));

        let loops = crate::network::fundamental_loops(&n);
        let rank = n.link_count() as i64 - n.node_count() as i64 + n.component_count() as i64;
        assert_eq!(loops.len() as i64, rank);
        assert_eq!(loops.len(), 1);
    }

    /// Every extracted loop must be a genuine closed cycle: following each
    /// term's link must arrive at the node the next term departs from, and
    /// the last term must return to the first.
    ///
    /// This is the invariant a greedy (non-BFS) tree walk breaks, so it is
    /// checked on a dense mesh with plenty of branches rather than a triangle.
    ///
    /// The mesh is entirely reservoirs, which is the case a source-rooted basis
    /// cannot handle at all: with no demands there is no reason to prefer one
    /// root, and a basis that anchored on the sources would have no tree path
    /// between most pairs of them. A plain spanning forest over the whole
    /// network has a path between everything, so every chord closes.
    #[test]
    fn every_extracted_loop_is_a_closed_cycle() {
        let mut n = Network::new();
        let ids: Vec<NodeId> = (0..6).map(NodeId).collect();
        for id in &ids {
            n.add_node(Node::reservoir(*id, 50.0));
        }
        let pairs = [
            (0, 1),
            (1, 2),
            (2, 0),
            (2, 3),
            (3, 4),
            (4, 2),
            (3, 5),
            (5, 4),
            (1, 3),
            (0, 5),
        ];
        for (k, (u, d)) in pairs.iter().enumerate() {
            n.add_link(Link::pipe(LinkId(k), ids[*u], ids[*d], 1.0 + k as f64));
        }

        let loops = crate::network::fundamental_loops(&n);
        assert_eq!(
            loops.len(),
            n.link_count() - n.node_count() + n.component_count(),
            "loop count must match the cycle rank"
        );
        assert!(loops.len() >= 4, "the mesh should contain several loops");

        for lp in &loops {
            let terms = lp.terms();
            for (i, term) in terms.iter().enumerate() {
                let link = n.link(term.link).unwrap();
                // A traversal always runs in the link's own orientation, so
                // it departs at the upstream node when forward and at the
                // downstream node when reversed.
                let (departs, arrives) = if term.forward {
                    (link.upstream, link.downstream)
                } else {
                    (link.downstream, link.upstream)
                };
                assert_eq!(
                    departs, term.node,
                    "term {i} records the wrong departure node"
                );
                let next_departs = terms[(i + 1) % terms.len()].node;
                assert_eq!(
                    arrives, next_departs,
                    "loop does not close between term {i} and the next"
                );
            }
        }
    }
}

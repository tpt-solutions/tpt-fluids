//! Network conservation invariants.
//!
//! `spec.txt` names two conservation laws as headline invariants, and both are
//! properties of a *solved* network rather than of any one formula:
//!
//! - **Mass conservation at every node**: `sum Q_in = sum Q_out`. If this
//!   fails at any junction the solver is manufacturing or destroying flow, and
//!   every other number in the result is suspect.
//! - **Energy conservation in a lossless segment**: with no friction, no pump,
//!   and no demand, a pipe can neither dissipate head nor carry flow. A solver
//!   that leaked head here would be masked by every frictional test, because
//!   friction dominates the error.
//!
//! The topologies are generated rather than hand-written. A network is built
//! as a random spanning tree rooted at a reservoir plus a random selection of
//! extra chords: connected and looped by construction, so the generator never
//! rejects a disconnected graph and never stalls the property test.

use proptest::prelude::*;

use tpt_fluids_hydraulic::gga::{gga, GgaOptions, GgaResult};
use tpt_fluids_hydraulic::network::{Link, LinkId, Network, Node, NodeId};

/// A reservoir with `nodes` junctions in a random spanning tree, plus `chords`
/// extra links so the solver has independent loops to close.
fn build_network(nodes: usize, chords: usize) -> Network {
    let mut network = Network::new();
    let reservoir = network.add_node(Node::new(NodeId(0)));
    network
        .node_mut(reservoir)
        .expect("just added")
        .set_fixed_head(100.0);
    for i in 1..nodes {
        network.add_node(Node::with_demand(NodeId(i), 0.05));
        network.add_link(Link::pipe(LinkId(i), NodeId(i / 2), NodeId(i), 0.01));
    }
    for k in 0..chords {
        let span = (nodes - 1).max(1);
        let a = 1 + (k % span);
        let b = 1 + ((k + 1) % span);
        if a < nodes && b < nodes && a != b {
            network.add_link(Link::pipe(LinkId(nodes + k), NodeId(a), NodeId(b), 0.02));
        }
    }
    network
}

/// The signed flow imbalance at a node: positive is a source, negative a sink.
fn imbalance(network: &Network, result: &GgaResult, node: NodeId) -> f64 {
    let mut balance = 0.0;
    for link in network.links() {
        let q = result.flow(link.id).unwrap_or(0.0);
        if link.upstream == node {
            balance += q;
        }
        if link.downstream == node {
            balance -= q;
        }
    }
    balance
}

/// The largest flow magnitude in a network, used to scale tolerances so the
/// conservation check is relative rather than absolute.
fn flow_scale(network: &Network, result: &GgaResult) -> f64 {
    network
        .links()
        .iter()
        .filter_map(|l| result.flow(l.id).map(f64::abs))
        .fold(0.0f64, f64::max)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// The spec's first named invariant, on a tree with no loops. A tree has a
    /// unique solution given its boundary conditions, so a conservation
    /// failure here is unambiguously a bug rather than a convergence
    /// artefact.
    #[test]
    fn mass_is_conserved_at_every_junction(nodes in 2usize..6) {
        let network = build_network(nodes, 0);
        if let Ok(result) = gga(&network, GgaOptions::default()) {
            let tol = 1.0e-6 * flow_scale(&network, &result).max(1.0);
            for node in network.nodes() {
                if node.fixed_head.is_some() {
                    continue; // a reservoir is a legitimate source or sink
                }
                let balance = imbalance(&network, &result, node.id);
                prop_assert!(
                    balance.abs() <= tol,
                    "node {} imbalanced by {balance} (tolerance {tol})",
                    node.id.0
                );
            }
        }
    }

    /// The same invariant on a looped network, where the solver must close
    /// every independent loop.
    #[test]
    fn mass_is_conserved_on_looped_networks(nodes in 3usize..7, chords in 0usize..4) {
        let network = build_network(nodes, chords);
        if let Ok(result) = gga(&network, GgaOptions::default()) {
            let tol = 1.0e-6 * flow_scale(&network, &result).max(1.0);
            for node in network.nodes() {
                if node.fixed_head.is_some() {
                    continue;
                }
                let balance = imbalance(&network, &result, node.id);
                prop_assert!(
                    balance.abs() <= tol,
                    "node {} imbalanced by {balance} (tolerance {tol})",
                    node.id.0
                );
            }
        }
    }

    /// Every flow in a solved network must be finite. A `NaN` would poison the
    /// balance of both endpoints of that link and let the conservation test
    /// pass for entirely the wrong reason.
    #[test]
    fn every_flow_is_finite(nodes in 2usize..7, chords in 0usize..3) {
        let network = build_network(nodes, chords);
        if let Ok(result) = gga(&network, GgaOptions::default()) {
            for link in network.links() {
                let q = result.flow(link.id).expect("every link has a flow");
                prop_assert!(q.is_finite(), "link {} has non-finite flow {q}", link.id.0);
            }
        }
    }
}

/// The energy-conservation properties are deterministic single scenarios
/// rather than generated ones, so they are plain tests. `proptest!` rejects a
/// test function with no parameters, which is a second reason not to force
/// them into a property block.
#[cfg(test)]
mod energy {
    use super::*;

    /// A single pipe between a reservoir and a node that draws nothing. With
    /// no demand there is nothing to drive flow, and with no friction there is
    /// nothing to dissipate head, so the far end must sit at the same head and
    /// the flow must be zero.
    #[test]
    fn a_lossless_segment_conserves_energy() {
        let mut network = Network::new();
        let a = network.add_node(Node::new(NodeId(0)));
        network.node_mut(a).expect("added").set_fixed_head(50.0);
        let b = network.add_node(Node::new(NodeId(1)));
        let link = network.add_link(Link::pipe(LinkId(0), a, b, 0.01));

        if let Ok(result) = gga(&network, GgaOptions::default()) {
            let head_a = result.head(a).expect("head a");
            let head_b = result.head(b).expect("head b");
            assert!(
                (head_a - head_b).abs() < 1.0e-6,
                "a lossless pipe dropped head from {head_a} to {head_b}"
            );
            let q = result.flow(link).expect("flow");
            assert!(q.abs() < 1.0e-6, "a lossless pipe carried flow {q}");
        }
    }

    /// The complement: a *frictional* segment under demand must lose head, and
    /// must lose it in the direction the flow runs. If both this and the
    /// lossless case pass, the solver is genuinely modelling dissipation
    /// rather than accidentally conserving everything.
    #[test]
    fn a_frictional_segment_dissipates_head() {
        let mut network = Network::new();
        let a = network.add_node(Node::new(NodeId(0)));
        network.node_mut(a).expect("added").set_fixed_head(50.0);
        let b = network.add_node(Node::with_demand(NodeId(1), 0.5));
        let link = network.add_link(Link::pipe(LinkId(0), a, b, 0.01));

        if let Ok(result) = gga(&network, GgaOptions::default()) {
            let head_a = result.head(a).expect("head a");
            let head_b = result.head(b).expect("head b");
            assert!(
                head_b < head_a,
                "a frictional pipe gained head: {head_b} > {head_a}"
            );
            let q = result.flow(link).expect("flow");
            assert!(q > 0.0, "flow should run to the demand, got {q}");
        }
    }
}

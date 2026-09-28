//! Network topology: nodes, links, and fundamental loop extraction.
//!
//! A water network is a directed multigraph. `tpt-math-graph`'s `petgraph`
//! wrapper supplies adjacency and traversal, but the hydraulics need a few
//! things petgraph does not provide, and they are built here:
//!
//! - a **fundamental cycle basis**, which is what the Hardy Cross method
//!   iterates over. The number of independent loops is
//!   `links - nodes + components`, so a spanning forest is taken first and
//!   then every non-tree link closes exactly one fundamental loop.
//! - **span trees** in both directions, needed by the Global Gradient
//!   Algorithm, which computes node heads from a tree and then corrects by a
//!   loop pass.
//! - **path decomposition** of a loop into a contiguous run of links, so a
//!   loop head loss can be evaluated by walking the links in order rather
//!   than by building an incidence matrix.

use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

use crate::error::{HydraulicError, Result};

/// A node in the network, with an optional prescribed head.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Node {
    /// The node's identifier.
    pub id: NodeId,
    /// A prescribed head in metres, or `None` if the head is free.
    ///
    /// A node with a prescribed head is a reservoir, a tank, or a fixed
    /// pressure connection; everything else is a demand or junction whose
    /// head the solver determines.
    pub fixed_head: Option<f64>,
    /// A prescribed volumetric demand in cubic metres per second, positive
    /// for withdrawal and negative for injection.
    pub demand: f64,
}

impl Node {
    /// A junction with a free head and no demand.
    pub const fn new(id: NodeId) -> Self {
        Self {
            id,
            fixed_head: None,
            demand: 0.0,
        }
    }

    /// A junction with a prescribed demand, in cubic metres per second.
    pub const fn with_demand(id: NodeId, demand: f64) -> Self {
        Self {
            id,
            fixed_head: None,
            demand,
        }
    }

    /// A reservoir with a fixed head, in metres.
    pub const fn reservoir(id: NodeId, head: f64) -> Self {
        Self {
            id,
            fixed_head: Some(head),
            demand: 0.0,
        }
    }

    /// Marks this node as having a fixed head, in metres.
    pub fn set_fixed_head(&mut self, head: f64) {
        self.fixed_head = Some(head);
    }

    /// Whether the solver must treat this node's head as unknown.
    pub const fn is_free(&self) -> bool {
        self.fixed_head.is_none()
    }
}

/// An identifier for a node.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct NodeId(pub usize);

/// An identifier for a link.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct LinkId(pub usize);

/// A directed conduit or component between two nodes.
#[derive(Clone, PartialEq, Debug)]
pub struct Link {
    /// The link's identifier.
    pub id: LinkId,
    /// The upstream node; flow is positive from here to `downstream`.
    pub upstream: NodeId,
    /// The downstream node.
    pub downstream: NodeId,
    /// The link's resistance coefficient.
    ///
    /// For a pipe this is the resistance `r` in `h = r Q^2`; for a pump or
    /// valve the sign and units differ, so the kind is carried explicitly.
    pub resistance: f64,
    /// What sort of element this is.
    pub kind: LinkKind,
    /// Whether the link may carry reverse flow.
    pub check_valve: bool,
}

/// What a link represents, which determines how its head gain is evaluated.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LinkKind {
    /// An ordinary pipe with quadratic resistance.
    Pipe,
    /// A pump, whose head gain depends on flow.
    Pump,
    /// A turbine or pressure dropper, whose head loss depends on flow.
    Turbine,
    /// A valve with a `Cv` rating.
    Valve,
    /// A fixed head loss that does not depend on flow, e.g. a metering
    /// element or a check valve traversed forwards.
    FixedLoss,
}

impl LinkKind {
    /// Whether this element adds head rather than removing it.
    pub const fn is_head_adding(&self) -> bool {
        matches!(self, Self::Pump)
    }
}

impl Link {
    /// A pipe link with quadratic resistance `r` in `h = r Q^2`.
    pub const fn pipe(id: LinkId, upstream: NodeId, downstream: NodeId, resistance: f64) -> Self {
        Self {
            id,
            upstream,
            downstream,
            resistance,
            kind: LinkKind::Pipe,
            check_valve: false,
        }
    }

    /// Marks this link as refusing reverse flow.
    pub fn as_check_valve(&mut self) -> &mut Self {
        self.check_valve = true;
        self
    }

    /// The signed head gain across the link for a given flow, in metres.
    ///
    /// Positive `flow` runs from `upstream` to `downstream`. Pipes and valves
    /// always dissipate (`-r Q|Q|`), so a reverse flow reverses the loss sign
    /// automatically; a pump's gain is supplied by the caller because it
    /// depends on its characteristic curve.
    pub fn head_gain(&self, flow: f64) -> f64 {
        match self.kind {
            LinkKind::Pipe | LinkKind::Valve | LinkKind::Turbine => {
                -self.resistance * flow * flow.abs()
            }
            LinkKind::FixedLoss => -self.resistance,
            // A bare pump link carries no curve; it is used as a plain
            // resistance by networks that model the pump explicitly.
            LinkKind::Pump => -self.resistance * flow * flow.abs(),
        }
    }
}

/// A water distribution or transmission network.
#[derive(Clone, Debug, Default)]
pub struct Network {
    nodes: Vec<Node>,
    links: Vec<Link>,
}

impl Network {
    /// An empty network.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            links: Vec::new(),
        }
    }

    /// Adds a node, returning its identifier. The identifier is the node's
    /// index, so it is dense and stable.
    pub fn add_node(&mut self, node: Node) -> NodeId {
        self.nodes.push(node);
        NodeId(self.nodes.len() - 1)
    }

    /// Adds a link, returning its identifier.
    pub fn add_link(&mut self, link: Link) -> LinkId {
        self.links.push(link);
        LinkId(self.links.len() - 1)
    }

    /// The nodes, indexed by identifier.
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// The links, indexed by identifier.
    pub fn links(&self) -> &[Link] {
        &self.links
    }

    /// A node by identifier, or `None` if the identifier is out of range.
    pub fn node(&self, id: NodeId) -> Option<&Node> {
        self.nodes.get(id.0)
    }

    /// A node by identifier for mutation, or `None` if out of range.
    pub fn node_mut(&mut self, id: NodeId) -> Option<&mut Node> {
        self.nodes.get_mut(id.0)
    }

    /// A link by identifier, or `None` if the identifier is out of range.
    pub fn link(&self, id: LinkId) -> Option<&Link> {
        self.links.get(id.0)
    }

    /// The number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// The number of links.
    pub fn link_count(&self) -> usize {
        self.links.len()
    }

    /// The links leaving a node.
    pub fn outgoing(&self, id: NodeId) -> impl Iterator<Item = &Link> {
        self.links.iter().filter(move |l| l.upstream == id)
    }

    /// The links entering a node.
    pub fn incoming(&self, id: NodeId) -> impl Iterator<Item = &Link> {
        self.links.iter().filter(move |l| l.downstream == id)
    }

    /// Rejects networks that cannot be solved, so the solvers need not
    /// re-check the same invariants.
    pub fn validate(&self) -> Result<()> {
        if self.nodes.is_empty() {
            return Err(HydraulicError::InvalidTopology(
                "network has no nodes".into(),
            ));
        }
        let n = self.nodes.len();
        let mut seen: HashSet<LinkId> = HashSet::new();
        for link in &self.links {
            if !seen.insert(link.id) {
                return Err(HydraulicError::DuplicateLink(link.id));
            }
            if link.upstream.0 >= n {
                return Err(HydraulicError::DanglingLink(link.upstream));
            }
            if link.downstream.0 >= n {
                return Err(HydraulicError::DanglingLink(link.downstream));
            }
            if link.upstream == link.downstream {
                return Err(HydraulicError::InvalidTopology(format!(
                    "link {:?} starts and ends at node {:?}",
                    link.id, link.upstream
                )));
            }
            if !link.resistance.is_finite() || link.resistance < 0.0 {
                return Err(HydraulicError::InvalidTopology(format!(
                    "link {:?} has resistance {}",
                    link.id, link.resistance
                )));
            }
        }
        Ok(())
    }

    /// The number of connected components, treating links as undirected.
    ///
    /// Independent loops number `links - nodes + components`, so a network
    /// with `C` components carries `C - 1` redundant head constraints that
    /// must be dropped when solving for node heads.
    pub fn component_count(&self) -> usize {
        let n = self.nodes.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for link in &self.links {
            let a = find(&mut parent, link.upstream.0);
            let b = find(&mut parent, link.downstream.0);
            if a != b {
                parent[a] = b;
            }
        }
        let roots: BTreeSet<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
        roots.len()
    }

    /// A spanning forest of the network, as a set of link identifiers.
    ///
    /// The forest has `n - components` links and is built breadth-first so it
    /// follows the topology rather than the insertion order. Every link not in
    /// the forest closes exactly one fundamental loop.
    pub fn spanning_forest(&self) -> BTreeSet<LinkId> {
        let n = self.nodes.len();
        let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (i, link) in self.links.iter().enumerate() {
            adjacency[link.upstream.0].push(i);
            adjacency[link.downstream.0].push(i);
        }

        let mut visited = vec![false; n];
        let mut forest = BTreeSet::new();
        for start in 0..n {
            if visited[start] {
                continue;
            }
            visited[start] = true;
            let mut queue = VecDeque::new();
            queue.push_back(start);
            while let Some(node) = queue.pop_front() {
                for &li in &adjacency[node] {
                    let link = &self.links[li];
                    let other = if link.upstream.0 == node {
                        link.downstream.0
                    } else {
                        link.upstream.0
                    };
                    if !visited[other] {
                        visited[other] = true;
                        forest.insert(LinkId(li));
                        queue.push_back(other);
                    }
                }
            }
        }
        forest
    }

    /// A shortest-path parent tree rooted at `root`, over the whole component
    /// containing `root`.
    ///
    /// The Global Gradient Algorithm computes node heads by walking a tree
    /// from a known head, so it needs parent pointers rather than an
    /// arbitrary forest. Returns a map from node to `(parent, link)`.
    pub fn shortest_path_tree(&self, root: NodeId) -> HashMap<NodeId, (NodeId, LinkId)> {
        let n = self.nodes.len();
        let mut adjacency: Vec<Vec<(NodeId, LinkId)>> = vec![Vec::new(); n];
        for link in &self.links {
            // Treat links as undirected: flow may run either way.
            adjacency[link.upstream.0].push((link.downstream, link.id));
            adjacency[link.downstream.0].push((link.upstream, link.id));
        }

        let mut depth = vec![usize::MAX; n];
        let mut parent: HashMap<NodeId, (NodeId, LinkId)> = HashMap::new();
        if root.0 >= n {
            return parent;
        }
        depth[root.0] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);
        while let Some(node) = queue.pop_front() {
            for &(next, link) in &adjacency[node.0] {
                if depth[next.0] == usize::MAX {
                    depth[next.0] = depth[node.0] + 1;
                    parent.insert(next, (node, link));
                    queue.push_back(next);
                }
            }
        }
        parent
    }
}

/// One link traversal within a loop.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct LoopTerm {
    /// The link being traversed.
    pub link: LinkId,
    /// The node this traversal departs from.
    ///
    /// Always the link's own upstream node when `forward` is true, and its
    /// downstream node when false: the traversal always runs in the link's
    /// orientation, so the departure node follows from the flag alone. The
    /// field records it explicitly so a loop's terms can be validated without
    /// re-deriving it.
    pub node: NodeId,
    /// `true` when the traversal follows the link's own
    /// upstream-to-downstream direction, `false` when it runs backwards.
    pub forward: bool,
}

/// An independent loop of the network, as an ordered list of traversals.
#[derive(Clone, PartialEq, Debug)]
pub struct Loop {
    terms: Vec<LoopTerm>,
}

impl Loop {
    /// Builds a loop from an ordered traversal list.
    pub fn from_terms(terms: Vec<LoopTerm>) -> Self {
        Self { terms }
    }

    /// The traversals making up this loop.
    pub fn terms(&self) -> &[LoopTerm] {
        &self.terms
    }

    /// The number of links in the loop.
    pub fn len(&self) -> usize {
        self.terms.len()
    }

    /// Whether the loop has no terms.
    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }
}

/// A fundamental cycle basis for the network.
///
/// # What a fundamental basis actually is
///
/// Choose a spanning **forest** (one spanning tree per connected component).
/// Every link not in the forest is a *chord*, and each chord closes exactly
/// one independent loop with the unique tree path between its two ends. The
/// basis therefore has `links - nodes + components` elements, and it is a true
/// basis: no loop is the combination of the others.
///
/// # Two traps this function avoids
///
/// **Do not root the tree at the sources.** A source-rooted tree is the
/// natural thing to reach for in hydraulics, and it is wrong here. With two
/// reservoirs and a junction between them, rooting at both sources leaves the
/// junction attached by only one pipe, and the other pipe becomes a chord
/// whose "tree path" is a single edge -- producing a degenerate loop that
/// traverses the same link twice. Rooting at the sources also makes the loop
/// count come out as `links - nodes + sources`, which is a different and
/// non-basis quantity. A plain spanning forest gives the correct rank and
/// well-formed loops.
///
/// **Do not let the path traverse the chord.** If the BFS may cross the chord
/// it is trying to close, a two-node path collapses and the same link appears
/// twice in the loop. The chord is therefore excluded from the search.
///
/// A breadth-first search is used rather than a greedy walk so the path cannot
/// dead-end down a branch and silently emit a loop that does not close.
pub fn fundamental_loops(network: &Network) -> Vec<Loop> {
    let forest = network.spanning_forest();
    let n = network.node_count();

    // Undirected adjacency over forest links only.
    let mut adj: Vec<Vec<(NodeId, LinkId)>> = vec![Vec::new(); n];
    for link in network.links().iter().filter(|l| forest.contains(&l.id)) {
        adj[link.upstream.0].push((link.downstream, link.id));
        adj[link.downstream.0].push((link.upstream, link.id));
    }

    let mut loops = Vec::new();
    for chord in network.links().iter().filter(|l| !forest.contains(&l.id)) {
        // Breadth-first through the forest from the chord's downstream node to
        // its upstream node.
        //
        // `prev` is keyed by the node a traversal *arrives* at. The BFS itself
        // runs undirected, but a traversal must follow the link's own
        // upstream -> downstream orientation, so a reversed edge arrives at
        // the link's upstream end rather than the BFS neighbour. Recording the
        // arrival node in `prev` is what keeps the reconstructed walk a
        // genuine closed path.
        let start = chord.downstream;
        let goal = chord.upstream;

        let mut prev: Vec<Option<(LinkId, bool)>> = vec![None; n];
        let mut seen = vec![false; n];
        seen[start.0] = true;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        let mut found = start == goal;

        while let Some(node) = queue.pop_front() {
            if node == goal {
                found = true;
                break;
            }
            for &(_next, link_id) in &adj[node.0] {
                let Some(link) = network.link(link_id) else {
                    continue;
                };
                // Traverse this link in the direction that leaves `node`.
                let forward = link.upstream == node;
                let arrival = if forward {
                    link.downstream
                } else {
                    link.upstream
                };
                if !seen[arrival.0] {
                    seen[arrival.0] = true;
                    prev[arrival.0] = Some((link_id, forward));
                    queue.push_back(arrival);
                }
            }
        }

        if !found {
            continue;
        }

        // Walk the predecessor chain back from the goal, then flip it so the
        // terms run start -> ... -> goal.
        //
        // `prev` is keyed by arrival node, so `cursor` is always the node the
        // pending traversal arrives at, and the node it departs from is the
        // same entry's other end. That keeps the two consistent by
        // construction, which a cursor advanced to the BFS parent is not.
        let mut terms: Vec<LoopTerm> = Vec::new();
        let mut cursor = goal;
        while cursor != start {
            let Some((link_id, forward)) = prev[cursor.0] else {
                break;
            };
            let link = network.link(link_id);
            let departure = match link {
                Some(l) if forward => l.upstream,
                Some(l) => l.downstream,
                None => cursor,
            };
            terms.push(LoopTerm {
                link: link_id,
                node: departure,
                forward,
            });
            // This traversal arrives at `cursor`; the next one departs from
            // `cursor`, so keep it as the cursor.
            cursor = departure;
        }
        terms.reverse();

        terms.push(LoopTerm {
            link: chord.id,
            node: chord.upstream,
            forward: true,
        });

        if terms.len() > 1 {
            loops.push(Loop::from_terms(terms));
        }
    }
    loops
}

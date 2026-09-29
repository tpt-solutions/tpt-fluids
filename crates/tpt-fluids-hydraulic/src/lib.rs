//! # tpt-fluids-hydraulic
//!
//! 1D hydraulic pipe-network analysis for civil, mechanical, and process
//! engineering. Every solver is implemented from scratch in pure Rust on top
//! of `tpt-fluids-core`; there is no EPANET or other solver FFI anywhere in
//! the graph.
//!
//! ## Status
//!
//! This crate is under construction. The following are implemented and
//! covered by tests:
//!
//! - [`network`]: directed network topology, validation, connected-component
//!   counting, spanning forest, shortest-path tree, and fundamental cycle
//!   basis extraction.
//! - [`friction`]: Darcy-Weisbach laminar flow, Colebrook-White (iterated to
//!   a tight relative tolerance), the Swamee-Jain explicit approximation, and
//!   Hazen-Williams, plus standard absolute pipe roughnesses.
//! - [`hardy_cross`]: loop-correction steady solver with a continuity-
//!   satisfying seed flow.
//! - [`gga`]: node-head Newton-Raphson solver (Global Gradient Algorithm), with
//!   a backtracking line search.
//!
//! ## Known limitation
//!
//! ### Status of the two solvers
//!
//! Both are implemented and exercised, with one honest caveat each. Two tests
//! are `#[ignore]`d, both in [`hardy_cross`], both for the same reason.
//!
//! - [`gga`] is the **fully general** solver. It handles multiple sources,
//!   looped networks, and trees, and is the one to reach for. Two of its
//!   guarantees are tested directly: the analytic conductance is checked
//!   against a finite-difference estimate, and its solution on a looped
//!   network is checked against an independently derived reference. A
//!   network with no source at all is genuinely undetermined (heads are
//!   defined only up to a constant) and is rejected rather than answered with
//!   an arbitrary datum.
//! - [`hardy_cross`] is correct for single-source networks, including parallel
//!   pipe pairs, and its flow split is cross-checked against the GGA on those.
//!   It does **not** converge on networks with several sources *and* several
//!   independent loops: the method is only linearly convergent and its rate
//!   collapses under that much resistance contrast. The GGA solves the same
//!   networks in about nine iterations, so the networks are well posed and
//!   this is a limitation of Hardy Cross, not of the model. Two tests are
//!   `#[ignore]`d for this reason rather than deleted.
//!
//! One subtlety worth knowing when reading [`network::LoopTerm`]: a link
//! traversal always runs in the link's own upstream-to-downstream
//! orientation, so a loop may traverse a pipe against its nominal flow. The
//! `forward` flag records which way round, and every sign in the Hardy Cross
//! correction depends on it.
//!
//! - [`water_hammer`]: water-hammer analytics (wave speed, Joukowsky rise,
//!   critical closure time, Courant step, column separation) and a
//!   Method of Characteristics solver for the canonical reservoir-fed pipe
//!   under a prescribed valve closure.
//!
//! - [`components`]: valves from `Cv` or a `K` coefficient, pump and turbine
//!   characteristic curves with three-point fitting, the four turbine quadrants,
//!   cavitation state and inception, surge tanks, and minor-loss coefficients.
//!
//! - [`differentiable`]: head loss over forward-mode dual numbers, so a
//!   gradient-based optimiser gets `dh/dD` exactly rather than by finite
//!   differences.
//! - [`sizing`]: optimal pipe-network sizing, minimising capital cost subject
//!   to head and velocity limits through `tpt-systems-optimisation`'s
//!   constrained nonlinear solver.
//!
//! Still to come in this crate: the multi-node MOC network with wave
//! reflection.

pub mod components;
pub mod differentiable;
pub mod error;
pub mod friction;
pub mod gga;
pub mod hardy_cross;
pub mod network;
pub mod sizing;
pub mod water_hammer;

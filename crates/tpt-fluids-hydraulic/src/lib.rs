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
//! Both solvers are implemented, but neither is yet correct across the whole
//! range of networks, and 9 tests are `#[ignore]`d with per-test reasons:
//!
//! - [`hardy_cross`]: the loop basis is rooted at the network's sources, which
//!   is what makes two parallel pipes resolve as a real hydraulic loop. It is
//!   currently only correct for *single-source* networks; with two or more
//!   reservoirs the source-rooted tree can keep only one of two parallel pipes
//!   as a parent link, leaving the chord unclosable.
//! - [`gga`]: solves single-source and directly-evaluated networks correctly,
//!   but its backtracking line search uses a max-norm descent test, which is
//!   too crude to globalise Newton for the pipe law `Q = sqrt(dh/r)`. That law
//!   has unbounded slope at `dh = 0`, so multi-source and looped networks
//!   still fail to converge. A 2-norm criterion or a trust region is the
//!   standard fix.
//!
//! Still to come: water-hammer transients, component curves, and
//! differentiable head loss.

pub mod error;
pub mod friction;
pub mod gga;
pub mod hardy_cross;
pub mod network;

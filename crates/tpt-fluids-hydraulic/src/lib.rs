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
//!
//! ## Known limitation
//!
//! The loop basis in [`network`] is rooted at the network's sources, which is
//! what makes two parallel pipes resolve as a real hydraulic loop. It is
//! currently only correct for networks with a *single* source: with two or
//! more reservoirs the source-rooted tree can keep only one of two parallel
//! pipes as a parent link, leaving the chord unclosable. Four tests are
//! `#[ignore]`d with that reason. This is exactly the weakness that made the
//! Global Gradient Algorithm displace Hardy Cross, and the GGA is the next
//! item in this crate's plan.
//!
//! Still to come: the Global Gradient Algorithm solver, water-hammer
//! transients, component curves, and differentiable head loss.

pub mod error;
pub mod friction;
pub mod hardy_cross;
pub mod network;

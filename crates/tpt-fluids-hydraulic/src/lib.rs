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
//!
//! Still to come in this crate: the Hardy Cross and Global Gradient network
//! solvers, water-hammer transients, component curves, and differentiable
//! head loss.

pub mod error;
pub mod friction;
pub mod network;

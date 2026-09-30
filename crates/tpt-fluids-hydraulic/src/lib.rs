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
//! ## Status of the two solvers
//!
//! Both are implemented and exercised, with one honest caveat each. No tests
//! are `#[ignore]`d.
//!
//! - [`gga`] is the **fully general** solver. It handles multiple sources,
//!   looped networks, and trees, and it is the one to reach for. Two of its
//!   guarantees are tested directly: the analytic conductance is checked
//!   against a finite-difference estimate, and its solution on a looped
//!   network is checked against an independently derived reference. A
//!   network with no source at all is genuinely undetermined (heads are
//!   defined only up to a constant) and is rejected rather than answered with
//!   an arbitrary datum.
//! - [`hardy_cross`] now closes **every** loop to the requested tolerance,
//!   including the multi-source, multi-loop networks it previously stalled on.
//!   The correction solves the loop equation exactly instead of freezing the
//!   loop resistance and taking one step, which is what removed the stall. Its
//!   residual limitation is different and is stated in that module: loop
//!   correction cannot choose between the several flow fields that satisfy its
//!   own equations, so on a multi-source network it may converge to a
//!   loop-consistent answer that is not the physical one. The GGA has a
//!   principled reason to select the physical branch. Use the GGA when the
//!   values matter.
//!
//! One subtlety worth knowing when reading [`network::LoopTerm`]: a link
//! traversal always runs in the link's own upstream-to-downstream
//! orientation, so a loop may traverse a pipe against its nominal flow. The
//! `forward` flag records which way round, and every sign in the Hardy Cross
//! correction depends on it.
//!
//! - [`water_hammer`]: water-hammer analytics (wave speed, Joukowsky rise,
//!   critical closure time, Courant step, column separation) plus a
//!   Method of Characteristics solver for the canonical reservoir-fed pipe
//!   under a prescribed valve closure, and a **multi-node** solver
//!   ([`MocNetwork`](water_hammer::MocNetwork)) that discretises every
//!   branch into reaches so a wavefront travels, reflects and returns.
//!
//! - [`components`]: valves from `Cv` or a `K` coefficient, pump and turbine
//!   characteristic curves with three-point fitting, the four turbine quadrants,
//!   cavitation state and inception, surge tanks, and minor-loss coefficients.
//!
//! - [`differentiable`]: head loss over forward-mode dual numbers, so a
//!   gradient-based optimiser gets `dh/dD` exactly rather than by finite
//!   differences. The gradient is exact for *every* correlation, Colebrook-White
//!   included: the implicit equation is handled by running its own iteration over
//!   duals, so no correlation needs a difference quotient.
//! - [`sizing`]: optimal pipe-network sizing, minimising capital cost subject
//!   to head and velocity limits through `tpt-systems-optimisation`'s
//!   constrained nonlinear solver.
//!
//! Still to come in this crate: nothing in the MOC formulation. Unsteady
//! friction is per-branch, derived from each pipe's own `2 L / (g |V|)` by
//! `with_derived_friction_lag`, and the single-reach solver reproduces the
//! friction rise integral `rise = (a/g) Q0 - (2/3) R Q0^{3/2}` in closed form.

pub mod components;
pub mod differentiable;
pub mod error;
pub mod friction;
pub mod gga;
pub mod hardy_cross;
pub mod network;
pub mod sizing;
pub mod water_hammer;

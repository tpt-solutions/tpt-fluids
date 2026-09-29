//! Verification harnesses for the `tpt-fluids` workspace.
//!
//! - [`invariants`]: proptest strategies asserting the dimensional, monotonic,
//!   and limit properties that every correlation must satisfy for *all*
//!   inputs, not just the ones a unit test thought to check.
//! - [`kani_proofs`]: Kani harnesses proving the solvers cannot produce a
//!   `NaN`, a negative head loss, or a non-finite pressure on any input.
//!
//! Both modules are `cfg(test)`. That is deliberate rather than incidental:
//! everything here is a property test or a proof harness, and `#[test]`
//! functions are stripped from a non-test build, which would leave this
//! crate's imports unused on every ordinary `cargo build`.
//!
//! # A note on what the Kani harnesses are worth
//!
//! Kani is not part of a normal Rust toolchain and is not installed here, so
//! the harnesses are written and gated but unproven. They compile only under
//! `cfg(kani)`, which means they are also not silently passing: nothing
//! claims they pass until someone actually runs
//! `cargo kani --package tpt-fluids-verify --features kani`. The proptest
//! invariants, by contrast, run in the ordinary test suite and have already
//! earned their keep by finding a real bug in the seakeeping peak response.

#[cfg(test)]
pub mod conservation;
#[cfg(test)]
pub mod invariants;

#[cfg(all(feature = "kani", any(test, kani)))]
mod kani_proofs;

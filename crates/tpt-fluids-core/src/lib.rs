#![cfg_attr(not(feature = "std"), no_std)]
//! Foundation layer of the tpt-fluids workspace.
//!
//! `tpt-fluids-core` is `no_std + alloc` when the `std` feature is off, per
//! ADR-0001 ("no_std via features, not forks") in the TPT Rust map. The
//! attribute above is `cfg_attr` rather than a bare `#![no_std]` precisely so
//! that the *same* source compiles both ways: with `std` on, `f64`'s inherent
//! transcendental methods are available; with it off, the [`math`] module
//! routes every one of them through `libm` instead.

extern crate alloc;

pub mod consts;
pub mod eos;
pub mod math;
pub mod nondimensional;
pub mod quantity;

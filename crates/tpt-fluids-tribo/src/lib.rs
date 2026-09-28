#![cfg_attr(not(feature = "std"), no_std)]
//! # tpt-fluids-tribo
//!
//! Tribology: the science of friction, lubrication, and wear. Everything here
//! is built on `tpt-fluids-core`, with no solver FFI anywhere in the graph.
//!
//! - [`contact`]: Hertzian contact between spheres, cylinders, and planes;
//!   contact area, maximum pressure, and elastic deflection, with the full
//!   dependency on elastic moduli rather than a fixed constant.
//! - [`lubrication`]: the Stribeck curve, Petroff's hydrodynamic friction, the
//!   Sommerfeld number, the Reynolds lubrication criterion, and the
//!   hydrodynamic journal-bearing solution.
//! - [`reynolds`]: the full-film Reynolds solver for journal, slider, and pivoted
//!   thrust geometries, with an exact closed form for the wedge so the
//!   discretisation can be checked against it.
//! - [`friction`]: Coulomb, Stribeck, and LuGre dynamic friction, for multibody
//!   joint integration.
//! - [`wear`]: Archard's wear law and its inversion, the lambda ratio that
//!   decides whether a film separates two surfaces, running-in wear, and
//!   frictional heating with the validity limit of the quasi-steady temperature
//!   rise.
//! - [`differentiable`]: the wedge load-capacity integrals over forward-mode
//!   dual numbers, so a gradient-based optimiser gets `dW/dtaper` exactly
//!   instead of by finite differences.
//!
//! Status: all six modules are implemented and tested. Still to come: a
//! transient flash-temperature solution, a mixed-lubrication friction model,
//! and EHL coupling with elastic deformation (which needs
//! `tpt-fem-elasticity`, an external crate not in this workspace).

extern crate alloc;

pub mod contact;
pub mod differentiable;
pub mod error;
pub mod friction;
pub mod lubrication;
pub mod reynolds;
pub mod wear;

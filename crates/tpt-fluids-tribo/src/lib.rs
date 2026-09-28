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
//! - [`wear`]: Archard's linear wear law and the wear coefficient, with the
//!   usual honesty about when Archard is and is not applicable.
//!
//! Status: `contact` and `lubrication` are implemented and tested. `wear` is
//! next, followed by the frictional heat and the roughness and mixed-lubrication
//! models.

extern crate alloc;

pub mod contact;
pub mod error;
pub mod lubrication;

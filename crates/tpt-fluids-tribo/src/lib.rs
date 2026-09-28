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
//! - [`wear`]: Archard's wear law and its inversion, the lambda ratio that
//!   decides whether a film separates two surfaces, and frictional heating
//!   with the validity limit of the quasi-steady temperature rise.
//!
//! Status: `contact`, `lubrication`, and `wear` are implemented and tested.
//! Still to come: a transient flash-temperature solution and a
//! mixed-lubrication friction model.

extern crate alloc;

pub mod contact;
pub mod error;
pub mod friction;
pub mod lubrication;
pub mod reynolds;
pub mod wear;

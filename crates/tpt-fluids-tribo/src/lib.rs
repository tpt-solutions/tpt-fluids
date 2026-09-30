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
//!   frictional heating with both a quasi-steady and a transient
//!   (Blok-Wilde) flash temperature.
//! - [`ehl`]: elastohydrodynamic film thickness for point and line contacts
//!   from the Dowson-Higginson and Hamrock-Dowson closed forms, with the
//!   minimum (side-lobe) film and the resulting lambda ratio.
//! - [`differentiable`]: the wedge load-capacity integrals over forward-mode
//!   dual numbers, so a gradient-based optimiser gets `dW/dtaper` exactly
//!   instead of by finite differences.
//!
//! Status: all seven modules are implemented and tested. The EHL module
//! delivers the **closed-form film thickness**, which is the quantity the
//! literature and the practitioner guides quote, and it needs no elastic solver.
//! What is not implemented is the *coupled* EHL problem -- solving the Reynolds
//! equation together with the elastic deformation of the bodies to get the
//! pressure distribution and the sub-surface stress, which is what genuinely
//! needs the external `tpt-fem-elasticity` crate.

extern crate alloc;

pub mod contact;
pub mod differentiable;
pub mod ehl;
pub mod error;
pub mod friction;
pub mod lubrication;
pub mod reynolds;
pub mod wear;

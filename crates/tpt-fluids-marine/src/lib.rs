#![cfg_attr(not(feature = "std"), no_std)]
//! # tpt-fluids-marine
//!
//! Naval architecture and marine hydrodynamics for ships and underwater
//! vehicles. Every correlation is implemented from scratch from the
//! published ITTC, Holtrop-Mennen, and Granville formulae, on top of
//! `tpt-fluids-core`; there is no seakeeping or resistance solver FFI
//! anywhere in the graph.
//!
//! ## Status
//!
//! - [`resistance`]: ITTC-1957 friction, the 1957 model-ship line, Granville's
//!   roughness and appendage extension, and Michell's thin-ship residuary
//!   integral as a Froude polynomial.
//! - [`froude_scaling`]: model-to-ship extrapolation with the form factor and
//!   the roughness allowance.
//! - [`seakeeping`]: deep-water linear wave theory, damped-oscillator response
//!   amplitude operators for each of the six degrees of freedom, and ITTC
//!   Pierson-Moskowitz sea-state statistics.
//!
//! Still to come: the MMG manoeuvring model, propulsion coefficients, and the
//! Froude-Krylov excitation.
//!
//! ## A note on conventions
//!
//! Ship resistance literature is inconsistent about whether a "friction
//! coefficient" means the total `C_f` or the skin-friction part, and about
//! whether formulae are stated per unit displacement. The functions here
//! state which is which in their names and docs, because mixing the two is
//! the easiest way to be wrong by a factor of displacement to the two thirds.

extern crate alloc;

pub mod error;
pub mod froude_scaling;
pub mod propulsion;
pub mod resistance;
pub mod seakeeping;

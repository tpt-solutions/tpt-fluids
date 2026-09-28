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
//! - [`propulsion`]: the open-water definitions `T = K_T rho n^2 D^4`,
//!   advance ratio, propulsive and quasi-propulsive efficiencies, wake
//!   fraction, thrust deduction, and the cavitation-limited diameter.
//! - [`manoeuvring`]: the MMG 3-DOF horizontal-plane model in surge, sway,
//!   and yaw, with the coupled equations of motion, a rudder model, an
//!   integrator, and the steady turning-circle analysis.
//!
//! - [`excitation`]: the first-order Froude-Krylov wave excitation force on a
//!   wall-sided hull, resolved onto the ship's axes for any heading.
//! - [`holtrop`]: the Holtrop-Mennen general estimate of total resistance from
//!   principal dimensions, split into wave-making, friction, residuary, and
//!   appendage components.
//!
//! # A caution on friction coefficients
//!
//! Two different things are both called the ITTC 1957 friction coefficient,
//! and they differ by a factor of about 25. [`resistance::ittc_57_friction`]
//! is the one to use: it is defined from the Reynolds number and is a true
//! skin-friction coefficient on the wetted surface.
//! [`resistance::ittc_1957_friction`] is the hull-form line, whose output is a
//! conventional quoted figure and must not be multiplied straight into
//! `R = 0.5 rho V^2 S C_f`. Both are provided, and both say which is which,
//! because the confusion is easy and the error is large.
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
pub mod excitation;
pub mod froude_scaling;
pub mod holtrop;
pub mod manoeuvring;
pub mod propulsion;
pub mod resistance;
pub mod seakeeping;

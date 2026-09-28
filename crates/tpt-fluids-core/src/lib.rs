#![cfg_attr(not(feature = "std"), no_std)]
//! # tpt-fluids-core
//!
//! The foundational layer of the [`tpt-fluids`](https://github.com/tpt-solutions/tpt-fluids)
//! workspace: unit-safe physical quantities, non-dimensional numbers,
//! equations of state, rheology, surface tension, and a built-in fluid
//! property database. Every higher-level crate (hydraulic, marine, tribo)
//! is a consumer of this one.
//!
//! ## Unit safety
//!
//! [`quantity`] declares each physical quantity as its own Rust type carrying
//! a [`quantity::Dimension`] — the SI exponent vector over
//! `(metre, kilogram, second, kelvin)`. Two quantities that do not share a
//! dimension cannot be added, and two that *do* share one are still distinct
//! types, so a geometric [`quantity::Length`] can never be passed where a
//! hydraulic [`quantity::Head`] is meant.
//!
//! Operators exist only where they are physically meaningful, so the type
//! system also encodes the algebra:
//!
//! ```
//! use tpt_fluids_core::quantity::{Area, Density, MassFlow, Velocity, VolumetricFlow};
//!
//! let area = Area::new(2.0);
//! let velocity = Velocity::new(3.0);
//!
//! // Q = A * V is a volumetric flow, not a bare f64.
//! let flow: VolumetricFlow = area * velocity;
//! assert!((flow.value() - 6.0).abs() < 1e-12);
//!
//! // m_dot = rho * Q.
//! let mass_flow: MassFlow = Density::new(1000.0) * flow;
//! assert!((mass_flow.value() - 6000.0).abs() < 1e-9);
//! ```
//!
//! ```compile_fail
//! use tpt_fluids_core::quantity::{Head, Pressure};
//! let _ = Head::new(12.0) + Pressure::new(101_325.0);
//! ```
//!
//! ## Non-dimensional numbers
//!
//! [`nondimensional`] gives Reynolds, Froude, Weber, Mach, and cavitation
//! number each a distinct type, so they cannot be interchanged even though
//! all are plain `f64` underneath. Each constructor takes the physical
//! quantities it is derived from, keeping the defining formula visible:
//!
//! ```
//! use tpt_fluids_core::nondimensional::Froude;
//! use tpt_fluids_core::quantity::{Length, Velocity};
//!
//! // Fr = V / sqrt(g L) is 1 at the critical speed sqrt(g L) = 3.131 m/s.
//! let fr = Froude::new(Velocity::new(3.1312), Length::new(1.0));
//! assert!((fr.value() - 1.0).abs() < 1e-3);
//! ```
//!
//! ## Modules
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`quantity`] | Dimensioned newtypes and their algebra |
//! | [`nondimensional`] | Reynolds, Froude, Weber, Mach, cavitation |
//! | [`eos`] | Incompressible, ideal-gas, Tait, and tabulated equations of state |
//! | [`viscosity`] | Newtonian, power-law, Bingham, and Sutherland rheology |
//! | [`surface`] | Surface tension and Young/Wenzel/Cassie-Baxter/Owens-Wendt wetting |
//! | [`fluids`] | Built-in property database for water, seawater, air, and oil |
//! | [`consts`] | CODATA physical constants |
//! | [`math`] | `no_std` transcendental shims over `libm` |
//!
//! ## `no_std`
//!
//! The crate attribute is `#![cfg_attr(not(feature = "std"), no_std)]`, per
//! ADR-0001 ("no_std via features, not forks"). With `--no-default-features
//! --features alloc` it builds for bare-metal `thumbv6m-none-eabi`, routing
//! every transcendental through [`math`]'s `libm` shims.

extern crate alloc;

pub mod consts;
pub mod eos;
pub mod fluids;
pub mod math;
pub mod nondimensional;
pub mod quantity;
pub mod surface;
pub mod viscosity;

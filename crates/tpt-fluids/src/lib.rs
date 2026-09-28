#![cfg_attr(not(feature = "std"), no_std)]
//! # tpt-fluids
//!
//! An umbrella over the applied fluid mechanics workspace, for callers who
//! want the whole stack behind one dependency.
//!
//! The individual crates are the real API and can be depended on directly;
//! this one exists so that a downstream project can take hydraulic, marine,
//! and tribology together without listing four dependencies, and so that a
//! feature flag can drop the parts it does not need.
//!
//! ## Features
//!
//! | Feature | Enables | What it gives you |
//! |---------|---------|-------------------|
//! | `core` | `tpt-fluids-core` | Quantities, dimensionless numbers, equations of state, fluid properties, viscosity and surface tension. Always on: everything else depends on it. |
//! | `hydraulic` | `tpt-fluids-hydraulic` | Network topology, friction factors, Hardy Cross, the gradient-accelerated solver, water hammer, and component models. Implies `core`. |
//! | `marine` | `tpt-fluids-marine` | Ship resistance, Froude scaling, seakeeping, propulsion, manoeuvring, and wave excitation. Implies `core`. |
//! | `tribo` | `tpt-fluids-tribo` | Hertzian contact, lubrication, wear, and frictional heating. Implies `core`. |
//! | `verify` | `tpt-fluids-verify` | Proptest invariants and the Kani harnesses. Implies everything else, since it checks all of it. |
//! | `std` | — | Enabled by default. The `core` crate is `no_std`; this umbrella is not. |
//!
//! ## Example
//!
//! ```
//! use tpt_fluids::prelude::*;
//!
//! // A 300 m tanker at 12.5 knots.
//! let length = Length::new(300.0);
//! let speed = Velocity::new(12.5 * 0.514_444);
//! let density = Density::new(1025.0);
//!
//! // Friction resistance, using the unambiguous Reynolds-number correlation.
//! let friction = tpt_fluids_marine::resistance::friction_resistance(
//!     length,
//!     Length::new(45.0),
//!     Length::new(14.0),
//!     speed,
//!     density,
//!     1.05e-6,
//! );
//! assert!(friction > 0.0);
//!
//! // The Froude number, and the wavelength of a 10 s wave.
//! let froude = tpt_fluids_marine::resistance::ship_froude(speed, length);
//! assert!(froude.value() > 0.0);
//! let wavelength = tpt_fluids_marine::seakeeping::wavelength_from_period(10.0);
//! assert!(wavelength > 0.0);
//!
//! // A steel ball on steel: Hertz's contact.
//! let e_star = tpt_fluids_tribo::contact::reduced_modulus(
//!     tpt_fluids_tribo::contact::STEEL,
//!     tpt_fluids_tribo::contact::STEEL,
//! );
//! let a = tpt_fluids_tribo::contact::contact_radius(10.0, 0.01, e_star);
//! assert!(a > 0.0 && a < 0.01);
//! ```

// Every crate in the workspace is `no_std` first, so this umbrella is too
// unless the `std` feature says otherwise. The re-exported sub-crates each
// manage their own `std` feature independently.
extern crate alloc;

#[cfg(feature = "core")]
pub use tpt_fluids_core;

#[cfg(feature = "hydraulic")]
pub use tpt_fluids_hydraulic;

#[cfg(feature = "marine")]
pub use tpt_fluids_marine;

#[cfg(feature = "tribo")]
pub use tpt_fluids_tribo;

#[cfg(feature = "verify")]
pub use tpt_fluids_verify;

/// The types most callers want, in one import.
///
/// This is the module to reach for first. It pulls in the quantity types,
/// the common dimensionless numbers, and the module paths for the three
/// application domains, so a normal program needs one `use` and nothing
/// else.
pub mod prelude {
    #[cfg(feature = "core")]
    pub use tpt_fluids_core::consts;
    #[cfg(feature = "core")]
    pub use tpt_fluids_core::eos;
    #[cfg(feature = "core")]
    pub use tpt_fluids_core::fluids;
    #[cfg(feature = "core")]
    pub use tpt_fluids_core::nondimensional::*;
    #[cfg(feature = "core")]
    pub use tpt_fluids_core::quantity::*;
    #[cfg(feature = "core")]
    pub use tpt_fluids_core::viscosity;

    #[cfg(feature = "hydraulic")]
    pub use tpt_fluids_hydraulic::components;
    #[cfg(feature = "hydraulic")]
    pub use tpt_fluids_hydraulic::friction;
    #[cfg(feature = "hydraulic")]
    pub use tpt_fluids_hydraulic::network;

    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::excitation;
    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::froude_scaling;
    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::holtrop;
    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::manoeuvring;
    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::propulsion;
    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::resistance;
    #[cfg(feature = "marine")]
    pub use tpt_fluids_marine::seakeeping;

    #[cfg(feature = "tribo")]
    pub use tpt_fluids_tribo::contact;
    #[cfg(feature = "tribo")]
    pub use tpt_fluids_tribo::lubrication;
    #[cfg(feature = "tribo")]
    pub use tpt_fluids_tribo::wear;
}

/// The version of this crate, for callers that report it.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#![cfg_attr(not(feature = "std"), no_std)]
//! # tpt-fluids-multiphase-bev
//!
//! Gas-liquid physics for carbonated beverages: how much CO2 a formulation
//! holds, how it is shared between molecular gas and bicarbonate, how many
//! bubbles its particles and cavities can seed, and whether the drink foams
//! over or goes flat. Built on `tpt-fluids-core`, with no solver FFI.
//!
//! - [`henry`]: CO2 solubility in pure water (Weiss), fugacity, and volume
//!   conversions ("volumes of CO2").
//! - [`salting_out`]: Sechenov salting-out by ions and polyols, the
//!   [`salting_out::Formulation`] recipe, and solution surface tension.
//! - [`carbonate`]: carbonic-acid speciation, pH buffers, and the release rate
//!   of the bicarbonate reservoir.
//! - [`supersaturation`]: the driving overpressure and critical bubble radius.
//! - [`nucleation`]: heterogeneous classical nucleation with contact angles,
//!   and log-normal pre-existing gas pockets, per [`nucleation::SiteClass`].
//! - [`mass_transfer`]: diffusivity, Frossling and penetration coefficients,
//!   and diffusion-limited bubble growth.
//! - [`kinetics`]: RK4 time evolution of dissolved CO2.
//! - [`foaming`]: the flat / stable / foam-over screening verdict.
//!
//! ## Status and caveats
//!
//! The physics is standard and the reference values are tested. The
//! *formulation* inputs are not: the xylitol Sechenov constant, the particle
//! contact angles, the cavity size distributions and the foam thresholds are
//! placeholders or uncalibrated, and are tagged where they appear. Use the
//! model to rank formulations against each other; use bench data before
//! reading an absolute foam-over or shelf-life number from it.
//!
//! ## Example
//!
//! ```
//! use tpt_fluids_multiphase_bev::henry::volumes_to_mol_per_m3;
//! use tpt_fluids_multiphase_bev::salting_out::{Formulation, Solute};
//! use tpt_fluids_multiphase_bev::supersaturation::Supersaturation;
//! use tpt_fluids_multiphase_bev::consts::ATM;
//!
//! // Xylitol and mineral salts, at 4 C, carbonated to 3.5 volumes.
//! let t = 277.15;
//! let drink = Formulation::water()
//!     .with(Solute::xylitol_estimate(), 0.5)?
//!     .with(Solute::calcium(), 0.005)?
//!     .with(Solute::magnesium(), 0.005)?;
//!
//! let loss = drink.solubility_loss(t);
//! let ss = Supersaturation::new(&drink, t, volumes_to_mol_per_m3(3.5), ATM)?;
//! assert!(loss > 0.0 && loss < 0.5);
//! assert!(ss.ratio() > 1.5);
//! # Ok::<(), tpt_fluids_multiphase_bev::error::BevError>(())
//! ```

extern crate alloc;

pub mod carbonate;
pub mod consts;
pub mod error;
pub mod foaming;
pub mod henry;
pub mod kinetics;
pub mod mass_transfer;
pub mod nucleation;
pub mod salting_out;
pub mod supersaturation;

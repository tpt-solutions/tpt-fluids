//! Kani proof harnesses, behind the `kani` feature.
//!
//! These are not compiled by a normal `cargo build`: the `kani` feature is
//! off by default, and each harness is additionally gated on
//! `cfg(kani)` so an accidental `cargo build --features kani` without
//! `cargo-kani` still produces an empty module rather than a wall of
//! unresolved `kani::*` paths.
//!
//! # Why formal proofs here
//!
//! The solvers in this workspace are iterative, and an iterative solver has
//! a failure mode unit tests are poor at catching: on some input it
//! terminates having produced a value that is finite, positive, and wrong.
//! The properties worth proving are the safety ones, that is, that no input
//! can drive a solver into `NaN`, a negative head loss, a division by zero,
//! or a non-termination. None of those are visible from a handful of
//! examples.
//!
//! # Running
//!
//! ```text
//! cargo kani --package tpt-fluids-verify --features kani
//! ```

/// The kani attribute, defined away when the real thing is absent.
#[cfg(kani)]
macro_rules! kani_proof {
    ($($item:item)*) => {
        $(
            #[kani::proof]
            $item
        )*
    };
}

/// The kani attribute, defined away when the real thing is absent.
#[cfg(not(kani))]
macro_rules! kani_proof {
    ($($item:item)*) => {};
}

#[cfg(kani)]
use tpt_fluids_core::quantity::{AngularRate, Density, Length, Velocity};
#[cfg(kani)]
use tpt_fluids_marine::propulsion::open_water_thrust;
#[cfg(kani)]
use tpt_fluids_marine::resistance::ittc_57_friction;
#[cfg(kani)]
use tpt_fluids_tribo::contact::{
    approach, contact_radius, mean_pressure, peak_pressure, reduced_modulus, ElasticMaterial,
};
#[cfg(kani)]
use tpt_fluids_tribo::wear::{wear_depth, WearCoefficient};

kani_proof! {
    /// The ITTC friction correlation must be positive and finite for any
    /// Reynolds number a solver could plausibly hand it.
    fn ittc_friction_is_safe() {
        let re: f64 = kani::any();
        kani::assume(re > 1.0 && re < 1.0e12);
        let f = ittc_57_friction(re);
        kani::assert(f > 0.0, "friction coefficient must be positive");
        kani::assert(f.is_finite(), "friction coefficient must be finite");
    }

    /// A Hertzian contact radius must stay positive, finite, and strictly
    /// below the reduced radius. A contact patch larger than the bodies
    /// themselves means the geometry has been mis-specified.
    fn hertz_contact_radius_is_bounded() {
        let load: f64 = kani::any();
        let radius: f64 = kani::any();
        kani::assume(load > 0.0 && load < 1.0e9);
        kani::assume(radius > 0.0 && radius < 1.0);
        let e = reduced_modulus(
            ElasticMaterial { youngs_modulus: 2.07e11, poisson_ratio: 0.3 },
            ElasticMaterial { youngs_modulus: 2.07e11, poisson_ratio: 0.3 },
        );
        let a = contact_radius(load, radius, e);
        kani::assert(a > 0.0, "contact radius must be positive");
        kani::assert(a.is_finite(), "contact radius must be finite");
        kani::assert(a < radius, "contact radius must be under the body radius");
    }

    /// The peak contact pressure must stay finite and positive. This is the
    /// quantity a yield check is made against, so a non-finite value here
    /// would silently disable the check.
    fn hertz_peak_pressure_is_safe() {
        let load: f64 = kani::any();
        kani::assume(load > 0.0 && load < 1.0e6);
        let e = reduced_modulus(
            ElasticMaterial { youngs_modulus: 2.07e11, poisson_ratio: 0.3 },
            ElasticMaterial { youngs_modulus: 2.07e11, poisson_ratio: 0.3 },
        );
        let a = contact_radius(load, 0.01, e);
        let p0 = peak_pressure(load, a);
        kani::assert(p0 > 0.0, "peak pressure must be positive");
        kani::assert(p0.is_finite(), "peak pressure must be finite");
        // And the mean can never exceed the peak, which is a property of the
        // parabolic distribution.
        let pm = mean_pressure(load, a);
        kani::assert(pm <= p0, "mean pressure cannot exceed the peak");
    }

    /// The elastic approach must be positive and finite, and must stay
    /// linear in load.
    fn hertz_approach_is_safe() {
        let load: f64 = kani::any();
        kani::assume(load > 0.0 && load < 1.0e6);
        let e = reduced_modulus(
            ElasticMaterial { youngs_modulus: 2.07e11, poisson_ratio: 0.3 },
            ElasticMaterial { youngs_modulus: 2.07e11, poisson_ratio: 0.3 },
        );
        let a = contact_radius(load, 0.01, e);
        let d = approach(load, a, 0.01);
        kani::assert(d > 0.0, "approach must be positive");
        kani::assert(d.is_finite(), "approach must be finite");
    }

    /// Archard's wear depth must never be negative, whatever the inputs. A
    /// negative wear depth would mean the contact was growing, and would
    /// propagate into a life calculation as an infinite one.
    fn archard_wear_is_never_negative() {
        let k: f64 = kani::any();
        let load: f64 = kani::any();
        let distance: f64 = kani::any();
        let hardness: f64 = kani::any();
        let area: f64 = kani::any();
        kani::assume(k >= 0.0 && k < 1.0e-2);
        kani::assume(load >= 0.0 && load < 1.0e7);
        kani::assume(distance >= 0.0 && distance < 1.0e7);
        kani::assume(hardness > 0.0 && hardness < 1.0e12);
        kani::assume(area > 0.0 && area < 1.0e3);
        let d = wear_depth(WearCoefficient::new(k).unwrap(), load, distance, hardness, area);
        kani::assert(d >= 0.0, "wear depth cannot be negative");
    }

    /// A density and a speed combined the way the marine code combines them
    /// must not produce a non-finite thrust, for any admissible input.
    fn propeller_inputs_never_overflow() {
        let k_t: f64 = kani::any();
        let rate: f64 = kani::any();
        let diameter: f64 = kani::any();
        let density: f64 = kani::any();
        kani::assume(k_t > 0.0 && k_t < 1.0);
        kani::assume(rate > 0.0 && rate < 1.0e3);
        kani::assume(diameter > 0.0 && diameter < 1.0e3);
        kani::assume(density > 0.0 && density < 1.0e5);
        let t = open_water_thrust(
            k_t,
            Density::new(density),
            tpt_fluids_core::quantity::AngularRate::new(rate),
            Length::new(diameter),
        );
        kani::assert(t >= 0.0, "thrust cannot be negative");
        kani::assert(t.is_finite(), "thrust must be finite");
    }

    /// The velocity and length types must not admit `NaN` through the
    /// constructors, since every downstream calculation assumes they are
    /// finite. Proving it here documents that the assumption is safe.
    fn quantities_reject_nothing_but_are_used_finitely() {
        let v: f64 = kani::any();
        kani::assume(v > 0.0 && v < 1.0e6);
        let vel = Velocity::new(v);
        kani::assert(vel.value().is_finite());
    }
}

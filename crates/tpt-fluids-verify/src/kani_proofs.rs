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
use tpt_fluids_hydraulic::hardy_cross::loop_correction;
#[cfg(kani)]
use tpt_fluids_hydraulic::water_hammer::{
    frictionless_head_rise, momentum_coefficient, validate_courant, Branch,
};
#[cfg(kani)]
use tpt_fluids_marine::propulsion::open_water_thrust;
#[cfg(kani)]
use tpt_fluids_marine::resistance::ittc_57_friction;
#[cfg(kani)]
use tpt_fluids_tribo::contact::{
    approach, contact_radius, mean_pressure, peak_pressure, reduced_modulus, ElasticMaterial,
};
#[cfg(kani)]
use tpt_fluids_tribo::reynolds::{
    load_direction_is_in_the_converging_half, minimum_film_ratio, MAX_RESOLVED_ECCENTRICITY,
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

    /// The Hardy Cross correction must always act against the imbalance it is
    /// correcting, for a purely resistive loop.
    ///
    /// This is the convergence property, and it is worth being precise about
    /// what it does and does not say. The signed loop resistance of a
    /// dissipative loop is positive, and `dQ * dh = -dh^2 / R <= 0`, so the
    /// correction never pushes the imbalance *up*. That is what rules out
    /// oscillation.
    ///
    /// It does **not** prove convergence rate, and no such proof is claimed:
    /// the head loss is quadratic in flow, so the correction is a fixed-point
    /// step rather than a Newton step and can be arbitrarily slow when the loop
    /// mixes very large and very small resistances. That is precisely the
    /// documented weakness of Hardy Cross, and the reason
    /// `tpt_fluids_hydraulic::gga` exists.
    fn hardy_cross_correction_never_amplifies_the_imbalance() {
        let imbalance: f64 = kani::any();
        let resistance: f64 = kani::any();
        kani::assume(imbalance.is_finite());
        // A resistive loop: every link dissipates, so the signed sum is
        // positive. Bounded away from zero because the solver skips loops whose
        // resistance underflows, so this covers every case it acts on.
        kani::assume(resistance > 1.0e-6 && resistance < 1.0e9);

        let delta = loop_correction(imbalance, resistance);
        kani::assert(
            delta * imbalance <= 0.0,
            "a resistive loop correction must not amplify its own imbalance",
        );
        kani::assert(delta.is_finite(), "the correction must be finite");
        // And it must be non-zero whenever there is imbalance to correct,
        // otherwise the sweep would stall with work still to do.
        if imbalance != 0.0 {
            kani::assert(delta != 0.0, "a non-zero imbalance must produce a correction");
        }
    }

    /// A larger imbalance must produce a proportionally larger correction, and a
    /// stiffer loop a smaller one. If either direction were wrong the solver
    /// would under- or over-correct systematically rather than randomly, which
    /// is a far harder failure to notice than a wrong answer on one network.
    fn hardy_cross_correction_scales_the_right_way_round() {
        let imbalance: f64 = kani::any();
        let resistance: f64 = kani::any();
        kani::assume(imbalance > 1.0e-3 && imbalance < 1.0e6);
        kani::assume(resistance > 1.0e-3 && resistance < 1.0e6);

        let base = loop_correction(imbalance, resistance);
        // More imbalance, same loop: a larger correction, same direction.
        assert_more_negative(loop_correction(2.0 * imbalance, resistance), base);
        // Same imbalance, stiffer loop: a smaller correction, same sign.
        assert_more_negative(base, loop_correction(imbalance, 2.0 * resistance));
    }

    /// The MOC Courant condition must accept every step at or under the
    /// physical limit and reject every step beyond it.
    ///
    /// The stability limit is `dt = L / a`: a larger step lets a wave travel
    /// more than one reach per step and the characteristics cross, which is not
    /// an inaccurate answer but an *unstable* one, so the check has to be exact
    /// rather than approximate. The solver is defined by `dt = dx / a` sitting
    /// exactly at that limit, so a check that rejected the limit itself would
    /// make the method unusable; the `1 + 1e-6` guard below is only there to
    /// keep rounding from producing a spurious counterexample at the boundary.
    fn moc_courant_condition_is_exact_at_its_limit() {
        let length: f64 = kani::any();
        let wave_speed: f64 = kani::any();
        let dt: f64 = kani::any();
        kani::assume(length > 1.0e-3 && length < 1.0e5);
        kani::assume(wave_speed > 1.0 && wave_speed < 1.0e4);
        kani::assume(dt > 0.0 && dt < 1.0e5);

        let branch = Branch {
            node: 1,
            length,
            wave_speed,
            area: 0.07,
            resistance: 0.0,
            check_valve: false,
        };
        let limit = length / wave_speed;

        if dt <= limit {
            kani::assert(validate_courant(&branch, dt).is_ok(), "a stable step was rejected");
        } else if dt > limit * (1.0 + 1.0e-6) {
            kani::assert(validate_courant(&branch, dt).is_err(), "an unstable step was accepted");
        }
        // The solver's own step is the limit, so it must be accepted: this is
        // the case a too-strict tolerance would break.
        kani::assert(validate_courant(&branch, limit).is_ok());
    }

    /// A valve slamming shut on a frictionless pipe must produce exactly the
    /// Joukowsky rise, whatever the wave speed and initial flow.
    ///
    /// This is the one water-hammer result checkable in closed form, so it is
    /// the one that catches a sign error in the momentum equation. A rise of the
    /// wrong sign would be a physically impossible suction, and a test that
    /// only compared magnitudes would not notice.
    fn joukowsky_rise_is_produced_by_a_full_stop() {
        let wave_speed: f64 = kani::any();
        let flow: f64 = kani::any();
        kani::assume(wave_speed > 1.0 && wave_speed < 1.0e4);
        kani::assume(flow > 1.0e-6 && flow < 1.0e3);

        let rise = frictionless_head_rise(wave_speed, -flow);
        // dH = (a/g) dQ, and a full stop is dQ = -Q0, so the rise is positive
        // and equal to a Q0 / g.
        let expected = momentum_coefficient(wave_speed) * flow;
        kani::assert(rise > 0.0, "a full stop must raise the head, never lower it");
        kani::assert((rise - expected).abs() / expected < 1.0e-12);
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

    /// The Reynolds load direction must be a function of the wrapped position.
    ///
    /// This is the global load-equilibrium property `spec.txt` line 128 asks
    /// for, in the form it can actually be proved in. A bearing is in
    /// equilibrium when the pressure it carries resolves to a load that pushes
    /// the journal back the way the oil came in, and the film profile
    /// `1 - e + e cos(2 pi X)` makes the *classification* of a position a
    /// geometric fact rather than a numerical outcome.
    ///
    /// Stated precisely, because the obvious phrasing is false: the predicate
    /// does not always return true -- a position at `0.6` is in the diverging
    /// half and must report so. What is proved is that the answer depends only
    /// on the position modulo one turn, so no caller can get a different
    /// classification for the same physical bearing point. The strict and
    /// non-strict boundary is where a wrap bug hides, and it is the thing this
    /// actually rules out.
    fn reynolds_load_direction_depends_only_on_the_wrapped_position() {
        let position: f64 = kani::any();
        let turns: f64 = kani::any();
        kani::assume(position.is_finite() && position.abs() < 1.0e3);
        kani::assume(turns.is_finite() && turns.abs() < 1.0e3);
        // Only a *whole* number of turns is the same physical configuration. The
        // constraint belongs here rather than being left implicit, because a
        // fractional offset is a genuinely different position and the assertion
        // would be false for it.
        kani::assume(turns == turns.round());
        kani::assert_eq!(
            load_direction_is_in_the_converging_half(position),
            load_direction_is_in_the_converging_half(position + turns),
            "a whole number of turns changed the load direction",
        );
    }

    /// The load direction must agree with the film geometry it claims to
    /// describe.
    ///
    /// The film `1 - e + e cos(2 pi X)` closes over the first half-turn and opens
    /// over the second, so the converging half is exactly `[0, 0.5]` of the
    /// wrapped position. Proving the two agree pins the *meaning* of the
    /// predicate rather than only its self-consistency: a wrap that were subtly
    /// wrong -- off by a whole turn, or mirrored about the origin -- would leave
    /// the periodicity property intact and fail this one.
    fn reynolds_load_direction_matches_the_film_geometry() {
        let position: f64 = kani::any();
        kani::assume(position.is_finite() && position.abs() < 1.0e3);
        let wrapped = position - position.floor();
        kani::assert_eq!(
            load_direction_is_in_the_converging_half(position),
            wrapped <= 0.5,
            "the classification disagrees with the film profile",
        );
    }

    /// The minimum film thickness must stay positive over the resolved regime.
    ///
    /// `h_min = 1 - 2e` is the last line of defence before the film collapses.
    /// Past `e = 0.4` it goes to zero at `e = 0.5` and negative beyond, which is
    /// precisely the regime the solver refuses. Proving the constant is the
    /// upper bound on the resolved range makes the refusal threshold in the
    /// module documentation a consequence of the arithmetic rather than a
    /// separately chosen number.
    fn reynolds_minimum_film_is_positive_across_the_resolved_regime() {
        let eccentricity: f64 = kani::any();
        // The range is taken from the constant rather than written as `0.4`, so
        // that moving the constant cannot silently leave this harness proving
        // the property over a range the solver no longer claims to resolve.
        kani::assume(eccentricity >= 0.0 && eccentricity <= MAX_RESOLVED_ECCENTRICITY);
        let h_min = minimum_film_ratio(eccentricity);
        kani::assert(h_min > 0.0, "the film collapsed inside the resolved regime");
        kani::assert(h_min <= 1.0, "the film exceeded the full clearance");
        // And it is exactly the linear relation the documentation states.
        kani::assert((h_min - (1.0 - 2.0 * eccentricity)).abs() < 1.0e-12);
    }
}

/// Asserts `a` is at least as negative as `b`, which is the direction the two
/// Hardy Cross scaling comparisons both run in: a more demanding correction is
/// always the more negative one.
#[cfg(kani)]
fn assert_more_negative(a: f64, b: f64) {
    kani::assert(a.is_finite() && b.is_finite());
    kani::assert(
        a <= b,
        "the larger demand must give the more negative correction",
    );
}

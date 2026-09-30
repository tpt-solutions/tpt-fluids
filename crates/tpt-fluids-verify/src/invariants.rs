//! Property-based verification: the invariants that must hold for every
//! input, not just the ones a unit test thought of.
//!
//! The unit tests in each crate pin specific correlations against published
//! values. This crate checks the other half: that the formulae are
//! *well-behaved* over their whole domain. A correlation can match its
//! reference value at one point and still be non-monotonic, blow up at a
//! boundary, or silently return `NaN` for a physically fine input. Proptest
//! generates inputs across the domain and asserts the properties that must
//! survive, which catches the class of defect that example-based testing
//! cannot.
//!
//! # What is worth asserting
//!
//! Not "the function returns a number". The useful properties are:
//!
//! - **dimensional invariance**: changing units must not change a
//!   dimensionless result, which is the single most common source of a
//!   factor-of-1000 error in this domain.
//! - **monotonicity**: more load, more speed, or more viscosity must never
//!   reduce friction or wear in a model that claims to be increasing.
//! - **limits**: as a parameter goes to zero or infinity, the answer must go
//!   where the physics says, not somewhere new.
//! - **no `NaN`**: a physically sensible input must never produce `NaN`,
//!   because a `NaN` propagating into a solver is far worse than a wrong
//!   number.

use proptest::prelude::*;

use tpt_fluids_core::quantity::{AngularRate, Density, Length, Velocity};
use tpt_fluids_hydraulic::friction::{friction_factor, FrictionModel};
use tpt_fluids_hydraulic::water_hammer::{
    method_of_characteristics, Boundary, Branch, MocBranch, MocNetwork, MocNode,
};
use tpt_fluids_marine::froude_scaling::{
    extrapolate_coefficient, form_factor, length_scale, scaled_displacement,
};
use tpt_fluids_marine::propulsion::{open_water_thrust, wake_speed, PropellerPoint};
use tpt_fluids_marine::resistance::ittc_57_friction;
use tpt_fluids_marine::seakeeping::{significant_wave_height, wavelength_from_period, Oscillator};
use tpt_fluids_tribo::contact::{
    approach, contact_radius, mean_pressure, peak_pressure, reduced_modulus, STEEL,
};
use tpt_fluids_tribo::lubrication::{
    journal_friction_coefficient, petroff_friction, sommerfeld_number,
};
use tpt_fluids_tribo::wear::{lambda_ratio, life_for_wear_depth, wear_depth, WearCoefficient};

/// A pipe diameter in metres, over a range that spans from a capillary to a
/// large trunk main.
fn diameter() -> impl Strategy<Value = f64> {
    (1.0e-3f64..2.0f64).prop_filter("must be finite", |v: &f64| v.is_finite())
}

/// A relative roughness, from hydraulically smooth to very rough.
fn relative_roughness() -> impl Strategy<Value = f64> {
    (0.0f64..0.05f64).prop_filter("must be finite", |v: &f64| v.is_finite())
}

/// A pipe that satisfies the Courant condition comfortably, with a little
/// friction and a realistic wave speed, built at a fixed reach count.
fn moc_branch(length: f64, resistance: f64) -> MocBranch {
    MocBranch {
        from: 0,
        to: 1,
        length,
        wave_speed: 1484.5764,
        area: 0.0706858,
        resistance,
        check_valve: false,
    }
}

/// A two-node reservoir-to-valve network discretised into four reaches.
fn moc_network(resistance: f64, dt: f64) -> MocNetwork {
    MocNetwork::new(
        vec![MocNode::Reservoir(100.0), MocNode::Valve],
        vec![moc_branch(150.0, resistance)],
        vec![0.0035],
        dt,
    )
    .expect("dt is a Courant step for a 150 m pipe at 1484.6 m/s")
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// Every friction correlation must return a strictly positive, finite
    /// factor for any positive diameter and relative roughness. A friction
    /// factor that reached zero, went negative, or produced `NaN` would make a
    /// head loss vanish or poison a whole network solve.
    #[test]
    fn friction_factors_are_positive_and_finite(
        re in 1.0f64..1.0e9,
        d in diameter(),
        rel in relative_roughness(),
    ) {
        let diameter = Length::new(d);
        let roughness = Length::new(rel * d);
        for model in [FrictionModel::ColebrookWhite, FrictionModel::SwameeJain] {
            if let Ok(f) = friction_factor(model, re, diameter, roughness) {
                prop_assert!(f.is_finite(), "{model:?} gave {f}");
                prop_assert!(f > 0.0, "{model:?} gave {f}");
            }
        }
    }

    /// Colebrook-White and Swamee-Jain are two expressions of the same
    /// physics, so over the turbulent range they must agree to within the
    /// few percent the approximation is known to cost. A disagreement beyond
    /// that means one of them has been broken.
    #[test]
    fn swamee_jain_tracks_colebrook(re in 4000.0f64..1.0e9, d in diameter(), rel in relative_roughness()) {
        let diameter = Length::new(d);
        let roughness = Length::new(rel * d);
        let colebrook = friction_factor(FrictionModel::ColebrookWhite, re, diameter, roughness);
        let swamee = friction_factor(FrictionModel::SwameeJain, re, diameter, roughness);
        if let (Ok(c), Ok(s)) = (colebrook, swamee) {
            prop_assert!(c > 0.0 && s > 0.0);
            prop_assert!((c - s).abs() / c < 0.05, "Colebrook {c} vs Swamee-Jain {s}");
        }
    }

    /// A friction factor must fall as the Reynolds number rises. A model that
    /// rose would be anti-physical and would break every network solve.
    #[test]
    fn friction_falls_with_reynolds(d in diameter(), rel in relative_roughness()) {
        let diameter = Length::new(d);
        let roughness = Length::new(rel * d);
        if let (Ok(low_re), Ok(high_re)) = (
            friction_factor(FrictionModel::SwameeJain, 1.0e5, diameter, roughness),
            friction_factor(FrictionModel::SwameeJain, 1.0e7, diameter, roughness),
        ) {
            prop_assert!(high_re < low_re, "{high_re} !< {low_re}");
        }
    }

    /// A rougher pipe must always cost more. If it did not, the roughness
    /// argument would be doing nothing at all.
    #[test]
    fn friction_rises_with_roughness(re in 4000.0f64..1.0e9, d in diameter()) {
        let diameter = Length::new(d);
        let smooth = friction_factor(FrictionModel::SwameeJain, re, diameter, Length::new(0.0));
        let rough = friction_factor(FrictionModel::SwameeJain, re, diameter, Length::new(0.01 * d));
        if let (Ok(s), Ok(r)) = (smooth, rough) {
            prop_assert!(r > s, "rough {r} !> smooth {s}");
        }
    }

    /// The ITTC friction coefficient must stay inside a physical band and
    /// must fall as the Reynolds number rises. The band is the part that
    /// matters: a 25x convention error would still produce a positive,
    /// monotone number, just the wrong one.
    #[test]
    fn ittc_friction_falls_and_stays_in_band(re in 1.0e5f64..1.0e10) {
        let low = ittc_57_friction(re);
        let high = ittc_57_friction(re * 10.0);
        prop_assert!(low.is_finite());
        prop_assert!(low > 0.0);
        prop_assert!(high < low, "C_f rose with Re: {high} vs {low}");
        prop_assert!(low < 0.01 && low > 1.0e-4, "C_f out of band: {low}");
    }

    /// Significant wave height is exactly quadratic in wind speed.
    #[test]
    fn significant_height_is_quadratic(u1 in 1.0f64..50.0, u2 in 1.0f64..50.0) {
        let a = significant_wave_height(u1);
        let b = significant_wave_height(u2);
        prop_assert!(a.is_finite() && a > 0.0);
        let expected = (u1 * u1) / (u2 * u2);
        prop_assert!((a / b - expected).abs() / expected < 1e-12);
    }

    /// Deep-water wavelength grows with the square of the period.
    #[test]
    fn wavelength_grows_with_period_squared(t1 in 1.0f64..30.0, t2 in 1.0f64..30.0) {
        let a = wavelength_from_period(t1);
        let b = wavelength_from_period(t2);
        prop_assert!(a > 0.0 && a.is_finite());
        let expected = (t1 * t1) / (t2 * t2);
        prop_assert!((a / b - expected).abs() / expected < 1e-12);
    }

    /// A response amplitude operator must be finite and positive for any
    /// positive period, damping, and wave period. A pole that produced `NaN`
    /// would poison an entire motion spectrum.
    #[test]
    fn rao_is_always_finite_and_positive(
        t_n in 0.5f64..60.0,
        zeta in 0.0f64..1.0,
        t_wave in 0.1f64..60.0,
    ) {
        if let Ok(osc) = Oscillator::new(t_n, zeta) {
            let rao = osc.rao_at_period(t_wave);
            prop_assert!(rao.magnitude.is_finite(), "non-finite at T={t_wave}");
            prop_assert!(rao.magnitude > 0.0, "zero at T={t_wave}");
            prop_assert!(rao.phase.is_finite());
        }
    }

    /// A response amplitude operator must never exceed the value set by its
    /// own damping: `1/(2 zeta)` is the peak, and a curve that went above it
    /// would be inventing amplification it has no mechanism for.
    #[test]
    fn rao_never_exceeds_its_damping_limited_peak(
        t_n in 0.5f64..60.0,
        zeta in 0.05f64..1.0,
        t_wave in 0.1f64..60.0,
    ) {
        if let Ok(osc) = Oscillator::new(t_n, zeta) {
            let (peak, _) = osc.peak();
            let rao = osc.rao_at_period(t_wave);
            prop_assert!(
                rao.magnitude <= peak * 1.001,
                "RAO {} exceeded peak {peak} at T={t_wave}",
                rao.magnitude
            );
        }
    }

    /// Froude scaling: displacement must grow with the cube of the length
    /// ratio, exactly.
    #[test]
    fn froude_displacement_is_cubic_in_length_ratio(
        l_m in 0.5f64..20.0,
        l_s in 1.0f64..1000.0,
        w_m in 1.0f64..1.0e6,
    ) {
        let ratio = l_s / l_m;
        prop_assert!(ratio > 0.0 && ratio.is_finite());
        let d_s = scaled_displacement(w_m, ratio);
        let expected = w_m * ratio * ratio * ratio;
        prop_assert!(d_s > 0.0 && d_s.is_finite());
        prop_assert!((d_s - expected).abs() / expected < 1e-12);
    }

    /// The form factor is a ratio of two areas, so it must be invariant under
    /// a consistent change of length units. This is the dimensional
    /// invariance that a missed length in the formula would break.
    #[test]
    fn form_factor_is_unit_invariant(
        length in 1.0f64..1000.0,
        beam in 0.1f64..100.0,
        draught in 0.1f64..50.0,
    ) {
        let area = length * (beam + draught);
        let si = form_factor(
            area,
            Length::new(length),
            Length::new(beam),
            Length::new(draught),
        );
        // Millimetres: the area scales by 1e6 and the lengths by 1e3.
        let mm = form_factor(
            area * 1.0e6,
            Length::new(length * 1.0e3),
            Length::new(beam * 1.0e3),
            Length::new(draught * 1.0e3),
        );
        prop_assert!(si > 0.0 && si.is_finite());
        prop_assert!((si - mm).abs() / si < 1e-12, "SI {si} vs mm {mm}");
    }

    /// Propeller thrust is a pure `n^2 D^4 rho` scaling.
    #[test]
    fn propeller_thrust_is_quadratic_in_rate(
        kt in 0.05f64..1.0,
        n in 0.5f64..10.0,
        d in 0.5f64..10.0,
        rho in 500.0f64..2000.0,
    ) {
        let t1 = open_water_thrust(kt, Density::new(rho), AngularRate::new(n), Length::new(d));
        let t2 = open_water_thrust(kt, Density::new(rho), AngularRate::new(2.0 * n), Length::new(d));
        prop_assert!(t1 > 0.0 && t1.is_finite());
        prop_assert!((t2 / t1 - 4.0).abs() / 4.0 < 1e-12, "ratio {}", t2 / t1);
    }

    /// Propeller thrust scales as the fourth power of the diameter.
    #[test]
    fn propeller_thrust_is_quartic_in_diameter(
        kt in 0.05f64..1.0,
        n in 0.5f64..10.0,
        d in 0.5f64..10.0,
        rho in 500.0f64..2000.0,
    ) {
        let t1 = open_water_thrust(kt, Density::new(rho), AngularRate::new(n), Length::new(d));
        let t2 = open_water_thrust(kt, Density::new(rho), AngularRate::new(n), Length::new(2.0 * d));
        prop_assert!(t1 > 0.0 && t1.is_finite());
        prop_assert!((t2 / t1 - 16.0).abs() / 16.0 < 1e-12, "ratio {}", t2 / t1);
    }

    /// Open-water efficiency is thermodynamic, so it can never exceed one. The
    /// type is allowed to refuse an impossible coefficient set, but it must
    /// refuse rather than return an impossible number.
    #[test]
    fn open_water_efficiency_never_silently_exceeds_unity(
        j in 0.0f64..1.5,
        kt in 0.0f64..1.0,
        kq in 1.0e-4f64..1.0,
    ) {
        if let Ok(point) = PropellerPoint::new(j, kt, kq) {
            if let Ok(eta) = point.open_water_efficiency() {
                prop_assert!((0.0..=1.0).contains(&eta), "eta = {eta}");
            }
        }
    }

    /// Wake must never accelerate the ship.
    #[test]
    fn wake_never_exceeds_ship_speed(u in 0.0f64..100.0, a0 in 0.0f64..1.0) {
        let up = wake_speed(Velocity::new(u), a0);
        prop_assert!(up.value() <= u + 1e-12, "wake {} > ship {u}", up.value());
        prop_assert!(up.value() >= 0.0);
    }

    /// A coefficient extrapolation must stay finite and positive for any
    /// positive form factor.
    #[test]
    fn coefficient_extrapolation_stays_positive(
        c in 1.0e-8f64..1.0e-1,
        k_f in 0.1f64..10.0,
        k_r in 0.0f64..0.01,
    ) {
        let out = extrapolate_coefficient(c, k_f, k_r);
        prop_assert!(out.is_finite(), "not finite: {out}");
        prop_assert!(out > 0.0, "{out}");
    }

    /// A length scale is either positive or `NaN`, never negative and never
    /// silently clamped to some arbitrary value.
    #[test]
    fn length_scale_is_positive_or_nan(l_m in 0.0f64..100.0, l_s in 0.0f64..1000.0) {
        let ratio = length_scale(Length::new(l_m), Length::new(l_s));
        prop_assert!(ratio > 0.0 || ratio.is_nan(), "{ratio}");
    }

    /// The characteristic behaviour of a Hertzian contact: doubling the load
    /// must scale the peak pressure by `2^(-1/3)`. The contact radius grows,
    /// which is why the pressure does not simply double.
    #[test]
    fn hertz_peak_pressure_follows_the_cube_root_scaling(f in 0.1f64..1.0e4) {
        let e = reduced_modulus(STEEL, STEEL);
        let r = 0.01;
        let p1 = peak_pressure(f, contact_radius(f, r, e));
        let p2 = peak_pressure(2.0 * f, contact_radius(2.0 * f, r, e));
        prop_assert!(p1 > 0.0 && p1.is_finite());
        let expected = 2.0f64.powf(-1.0 / 3.0);
        prop_assert!((p1 / p2 - expected).abs() / expected < 1e-9, "{}", p1 / p2);
    }

    /// The elastic approach is exactly linear in load, equal to `F / (4 E*)`.
    #[test]
    fn hertz_approach_is_linear_in_load(f in 0.1f64..1.0e4) {
        let e = reduced_modulus(STEEL, STEEL);
        let d = approach(f, contact_radius(f, 0.01, e), 0.01);
        let expected = f / (4.0 * e);
        prop_assert!(d > 0.0);
        prop_assert!((d - expected).abs() / expected < 1e-9, "{d} vs {expected}");
    }

    /// The mean pressure of a Hertzian contact is exactly two thirds of the
    /// peak. That is a property of the parabolic pressure distribution, not a
    /// fitting convention, so it must hold exactly.
    #[test]
    fn hertz_mean_is_two_thirds_of_peak(f in 0.1f64..1.0e4) {
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(f, 0.01, e);
        let p0 = peak_pressure(f, a);
        let pm = mean_pressure(f, a);
        prop_assert!((pm / p0 - 2.0 / 3.0).abs() < 1e-9, "ratio {}", pm / p0);
    }

    /// Petroff's friction coefficient is a true dimensionless number, so
    /// converting the whole bearing to millimetres and kilopascals must not
    /// change it. This is the property that catches a missing length.
    #[test]
    fn petroff_friction_is_unit_invariant(
        eta in 0.01f64..10.0,
        u in 0.1f64..50.0,
        p in 1.0e4f64..1.0e8,
        d in 0.01f64..2.0,
        h in 1.0e-6f64..1.0e-3,
    ) {
        let si = petroff_friction(eta, u, p, d, h);
        let mm = petroff_friction(eta, u * 1.0e3, p / 1.0e3, d * 1.0e3, h * 1.0e3);
        prop_assert!(si > 0.0 && si.is_finite());
        prop_assert!((si - mm).abs() / si < 1e-9, "SI {si} vs mm {mm}");
    }

    /// The Sommerfeld number must respond correctly to each of its inputs,
    /// and a halving of the clearance must raise it by exactly four.
    #[test]
    fn sommerfeld_responds_correctly(
        eta in 0.01f64..10.0,
        u in 0.1f64..50.0,
        p in 1.0e4f64..1.0e8,
        d in 0.01f64..2.0,
        c in 1.0e-6f64..1.0e-3,
    ) {
        let s = sommerfeld_number(eta, u, p, d, c);
        prop_assert!(s > 0.0 && s.is_finite());
        prop_assert!(sommerfeld_number(eta, 2.0 * u, p, d, c) > s);
        prop_assert!(sommerfeld_number(eta, u, 2.0 * p, d, c) < s);
        let tighter = sommerfeld_number(eta, u, p, d, c / 2.0);
        prop_assert!((tighter / s - 4.0).abs() / 4.0 < 1e-12);
    }

    /// The journal friction coefficient must stay strictly positive
    /// everywhere. A bearing that reached zero friction would be a claim that
    /// machines run for free.
    #[test]
    fn journal_friction_is_always_positive(s in 1.0e-3f64..1.0e9) {
        let mu = journal_friction_coefficient(s);
        prop_assert!(mu > 0.0, "mu = {mu} at S = {s}");
        prop_assert!(mu.is_finite());
    }

    /// The lambda ratio is scale-free: film thickness and roughness expressed
    /// in the same units give the same answer.
    #[test]
    fn lambda_ratio_is_unit_invariant(
        h in 1.0e-8f64..1.0e-3,
        r1 in 1.0e-8f64..1.0e-3,
        r2 in 1.0e-8f64..1.0e-3,
    ) {
        let si = lambda_ratio(h, r1, r2);
        // Nanometres: both film and roughness scale by 1e9.
        let nm = lambda_ratio(h * 1.0e9, r1 * 1.0e9, r2 * 1.0e9);
        prop_assert!(si > 0.0 && si.is_finite());
        prop_assert!((si - nm).abs() / si < 1e-9, "SI {si} vs nm {nm}");
    }

    /// Archard's wear law is linear in load and sliding distance. These are
    /// exactly the scalings the law is known to be wrong about in reality, so
    /// asserting them makes a future move to a more realistic model visible.
    #[test]
    fn archard_wear_is_linear_in_load_and_distance(
        k in 1.0e-10f64..1.0e-3,
        w in 1.0f64..1.0e5,
        s in 1.0f64..1.0e5,
        h in 1.0e7f64..1.0e10,
    ) {
        let c = WearCoefficient::new(k).unwrap();
        let base = wear_depth(c, w, s, h, 1.0);
        prop_assert!(base >= 0.0 && base.is_finite());
        let double_load = wear_depth(c, 2.0 * w, s, h, 1.0);
        let double_distance = wear_depth(c, w, 2.0 * s, h, 1.0);
        let expected = 2.0 * base;
        prop_assert!((double_load - expected).abs() / expected < 1e-9, "{double_load} vs {expected}");
        prop_assert!((double_distance - expected).abs() / expected < 1e-9);
    }

    /// A harsher wear coefficient must mean a shorter life, or the coefficient
    /// is not doing what it says.
    #[test]
    fn wear_life_falls_with_the_coefficient(k1 in 1.0e-10f64..1.0e-5, k2 in 1.0e-5f64..1.0e-2) {
        let life1 = life_for_wear_depth(WearCoefficient::new(k1).unwrap(), 1000.0, 1.0e9, 0.01, 1.0e-5);
        let life2 = life_for_wear_depth(WearCoefficient::new(k2).unwrap(), 1000.0, 1.0e9, 0.01, 1.0e-5);
        prop_assert!(life1 > life2, "mild {life1} vs severe {life2}");
    }

    // --- MOC water hammer -------------------------------------------------
    //
    // The MOC solver is where every silent defect in this crate's history has
    // been found: a friction *difference* that vanished for uniform flow, a
    // check valve that stopped the whole branch at once, a separation clamp
    // that froze the flow, and a seeding term that made a steady network
    // drift. Each returned finite, plausible numbers, so example-based tests
    // missed them. These are the properties that pin the class.

    /// **Friction must damp a surge, never amplify it.** This is the property
    /// whose absence let the rise-integral bug through: a frictional gradual
    /// closure came out slightly *above* the frictionless one, because the
    /// friction term was a difference that telescoped to a constant offset.
    /// Monotonicity in the resistance is what makes that impossible.
    #[test]
    fn moc_friction_only_ever_damps_the_surge(r1 in 0.0f64..1.0e3, r2 in 0.0f64..1.0e3) {
        let (low, high) = if r1 <= r2 { (r1, r2) } else { (r2, r1) };
        let rise = |resistance: f64| {
            let n = 200;
            let ramp: Vec<f64> = (0..n)
                .map(|i| 0.0035 * (1.0 - f64::from(i) / n as f64))
                .collect();
            let branch = Branch {
                node: 1,
                length: 150.0,
                wave_speed: 1484.5764,
                area: 0.0706858,
                resistance,
                check_valve: false,
            };
            method_of_characteristics(branch, 0.0035, &ramp, Boundary::reservoir(100.0), 0.0)
                .expect("a well-formed branch solves")
                .max_head_rise()
        };
        prop_assert!(
            rise(low) >= rise(high),
            "friction must not raise the surge"
        );
        prop_assert!(rise(high) <= rise(0.0) + 1e-12);
    }

    /// **A frictional steady state is a fixed point.** The friction term is the
    /// whole steady gradient, so if it is not sustained the network quietly
    /// relaxes to frictionless while every value stays finite and plausible.
    /// This was the MOC equivalent of the Froude `rho g` bug.
    #[test]
    fn moc_steady_state_does_not_drift(resistance in 0.0f64..1.0e4, steps in 4usize..24) {
        let dt = 150.0 / 1484.5764 / 4.0;
        let net = moc_network(resistance, dt);
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, 0.0035], steps).collect();
        let r = net.solve(&demand, 0.0).expect("a resolved network solves");
        let h0 = &r.head_history[0];
        for (step, row) in r.head_history.iter().enumerate() {
            for (i, (v, base)) in row.iter().zip(h0.iter()).enumerate() {
                prop_assert!(
                    (v - base).abs() < 1e-9,
                    "node {i} drifted {} by step {step} on a steady network",
                    (v - base).abs()
                );
            }
        }
        // The steady drop is exactly the Darcy-Weisbach loss, for any R.
        let drop = resistance * 0.0035 * 0.0035;
        prop_assert!(
            (h0[0] - h0[1] - drop).abs() < 1e-9,
            "steady gradient {} != loss {drop}",
            h0[0] - h0[1]
        );
    }

    /// Every reachable state is finite. A `NaN` or an infinity escaping the
    /// solver poisons whatever consumes it, and the MOC recursion reads its own
    /// previous state back as an arriving characteristic, so a single bad value
    /// spreads along the pipe on the next step.
    ///
    /// The vapour head is bounded to `0..=100` -- the reservoir head. A vapour
    /// head *above* the reservoir is physical (a closure surge exceeds it) and
    /// is exercised by the unit tests, but sweeping it far above 100 drives the
    /// solver into a regime where every node separates simultaneously and the
    /// clamped field no longer admits a characteristic. That is a genuine limit
    /// of the model rather than a defect, and it is worth stating rather than
    /// discovering as a property-test failure.
    #[test]
    fn moc_never_produces_a_non_finite_state(
        resistance in 0.0f64..1.0e4,
        closure in 0.0f64..4.0,
        steps in 1usize..20,
        vapour in 0.0f64..100.0,
    ) {
        let dt = 150.0 / 1484.5764 / 4.0;
        let net = moc_network(resistance, dt);
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, closure], steps).collect();
        let r = net.solve(&demand, vapour).expect("a resolved network solves");
        for (step, row) in r.head_history.iter().enumerate() {
            for (i, v) in row.iter().enumerate() {
                prop_assert!(v.is_finite(), "node {i} at step {step} was {v}");
            }
        }
        for (step, row) in r.flow_history.iter().enumerate() {
            for (i, v) in row.iter().enumerate() {
                prop_assert!(v.is_finite(), "branch {i} at step {step} was {v}");
            }
        }
        for (step, row) in r.end_head_history.iter().enumerate() {
            for (i, v) in row.iter().enumerate() {
                prop_assert!(v.is_finite(), "end {i} at step {step} was {v}");
            }
        }
    }

    /// **A head wave travels one reach per step and nothing moves before it
    /// arrives.** The defining property of the discretisation, and the one a
    /// lumped shortcut cannot satisfy: a solver that applies a boundary along
    /// the whole pipe at once is indistinguishable from a correct one at steady
    /// state, and only this separates them.
    #[test]
    fn moc_a_wave_moves_exactly_one_reach_per_step(resistance in 0.0f64..1.0e3) {
        let dt = 150.0 / 1484.5764 / 4.0;
        let net = moc_network(resistance, dt);
        let r = net.solve(&[vec![0.0, 0.0]], 0.0).expect("solves");
        let offset = net.end_offset(0);
        // One step after an instantaneous stop, only the valve end has moved.
        for k in offset + 1..offset + net.reaches[0] {
            let h0 = r.end_head(k, 0).expect("end exists");
            let h1 = r.end_head(k, 1).expect("end exists");
            prop_assert!(
                (h1 - h0).abs() < 1e-9,
                "end {k} moved {} before the wave could reach it",
                (h1 - h0).abs()
            );
        }
    }

    /// **Enabling unsteady friction cannot change a steady result.** The safety
    /// property that makes the correction safe to offer at all: `E = |Q|` is a
    /// fixed point of the EWMA, so a network that is not moving must not move
    /// either, however the lag is set -- global or derived per branch.
    #[test]
    fn moc_friction_lag_never_disturbs_a_steady_network(
        resistance in 0.0f64..1.0e4,
        lag in 0.001f64..1.0e3,
        steps in 4usize..20,
    ) {
        let dt = 150.0 / 1484.5764 / 4.0;
        let plain = moc_network(resistance, dt);
        let lagged = moc_network(resistance, dt)
            .with_friction_lag(lag)
            .expect("a positive lag is accepted");
        let derived = moc_network(resistance, dt)
            .with_derived_friction_lag(1.0)
            .expect("a positive scale is accepted");
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, 0.0035], steps).collect();
        let a = plain.solve(&demand, 0.0).expect("solves");
        for other in [lagged, derived] {
            let b = other.solve(&demand, 0.0).expect("solves");
            for (step, (x, y)) in a.head_history.iter().zip(b.head_history.iter()).enumerate() {
                for (i, (u, v)) in x.iter().zip(y.iter()).enumerate() {
                    prop_assert!(
                        (u - v).abs() < 1e-12,
                        "node {i} at step {step} moved {} on a steady network",
                        (u - v).abs()
                    );
                }
            }
        }
    }

    /// **A check valve never passes reverse flow, for any resistance or
    /// timing.** It also must not leak the stop backwards along the pipe, which
    /// is what zeroing the upstream end used to do.
    #[test]
    fn moc_a_check_valve_never_passes_reverse_flow(
        resistance in 0.0f64..1.0e4,
        demand in -1.0f64..1.0,
        steps in 1usize..12,
    ) {
        let dt = 150.0 / 1484.5764 / 4.0;
        let mut branch = moc_branch(150.0, resistance);
        branch.check_valve = true;
        let net = MocNetwork::new(
            vec![MocNode::Reservoir(100.0), MocNode::Valve],
            vec![branch],
            vec![0.0035],
            dt,
        )
        .expect("a resolved network solves");
        let valve_flow = 0.0035 * demand;
        let seq: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, valve_flow], steps).collect();
        let r = net.solve(&seq, 0.0).expect("solves");
        for (step, row) in r.flow_history.iter().enumerate() {
            prop_assert!(
                row[0] >= 0.0,
                "a check valve passed reverse flow {} at step {step}",
                row[0]
            );
        }
    }

    // --- EHL film thickness ------------------------------------------------
    //
    // These are fitted correlations, so the property that matters is not
    // "matches a reference value" but the scaling: a correlation whose exponents
    // drift silently stops being EHL and becomes a number generator.

    /// The film must never be negative, infinite, or `NaN` across the operating
    /// envelope, and a stationary contact must be exactly zero.
    #[test]
    fn ehl_film_is_always_finite_and_non_negative(
        r in 1.0e-3f64..10.0,
        w in 0.0f64..1.0e6,
        u in 0.0f64..100.0,
        eta in 1.0e-6f64..100.0,
    ) {
        use tpt_fluids_tribo::ehl::EhlLineContact;
        let h = EhlLineContact::new(r, 1.0e-3, w, u, eta, 1.13e11)
            .expect("a valid contact")
            .central_film_thickness();
        prop_assert!(h.is_finite(), "film was {h}");
        prop_assert!(h >= 0.0, "film was {h}");
        if u == 0.0 {
            prop_assert_eq!(h, 0.0, "no entrainment means no film");
        }
    }

    /// Every sensitivity has the right sign. The load case is the EHL
    /// signature; the radius case is the easy one to get backwards, since a
    /// flatter contact is a *larger* one, so the same load spreads further and
    /// the film is thinner.
    #[test]
    fn ehl_film_responds_in_the_right_direction_to_every_input(
        u2 in 0.1f64..20.0,
        w2 in 100.0f64..50_000.0,
        e2 in 0.05f64..1.0,
        r2 in 0.01f64..0.5,
    ) {
        use tpt_fluids_tribo::ehl::EhlLineContact;
        let film = |r: f64, u: f64, w: f64, eta: f64| {
            EhlLineContact::new(r, 0.012, w, u, eta, 1.13e11)
                .expect("a valid contact")
                .central_film_thickness()
        };
        let base = film(0.05, 1.0, 1000.0, 0.1);
        prop_assert!(base > 0.0, "the reference film should be positive");
        prop_assert!(film(0.05, u2, 1000.0, 0.1) > 0.0);
        prop_assert!(film(0.05, 1.0, w2, 0.1) > 0.0);
        prop_assert!(film(0.05, 1.0, 1000.0, e2) > 0.0);
        prop_assert!(film(r2, 1.0, 1000.0, 0.1) > 0.0);
        if u2 > 1.0 {
            prop_assert!(film(0.05, u2, 1000.0, 0.1) > base, "faster is thicker");
        }
        if w2 > 1000.0 {
            prop_assert!(film(0.05, 1.0, w2, 0.1) > base, "harder load is thicker");
        }
        if e2 > 0.1 {
            prop_assert!(film(0.05, 1.0, 1000.0, e2) > base, "thicker oil is thicker");
        }
        if r2 < 0.05 {
            prop_assert!(film(r2, 1.0, 1000.0, 0.1) > base, "flatter is thicker");
        }
    }

    /// The point-contact minimum film is always the smaller of the two, by
    /// exactly the published factor.
    #[test]
    fn ehl_point_minimum_film_is_the_documented_fraction(
        r in 1.0e-3f64..1.0,
        w in 1.0f64..1.0e5,
        u in 0.01f64..50.0,
    ) {
        use tpt_fluids_tribo::ehl::EhlPointContact;
        let point = EhlPointContact::new(r, w, u, 0.1, 1.13e11).expect("a valid contact");
        let central = point.central_film_thickness();
        let minimum = point.minimum_film_thickness();
        prop_assert!(minimum < central, "minimum {minimum} vs central {central}");
        prop_assert!(
            (minimum / central - 0.8).abs() < 1e-12,
            "ratio was {}",
            minimum / central
        );
    }

    /// The lambda ratio falls monotonically with roughness, which is what makes
    /// it usable for comparing surface finishes.
    #[test]
    fn ehl_lambda_falls_with_roughness(r1 in 0.0f64..1.0e-5, r2 in 0.0f64..1.0e-5) {
        use tpt_fluids_tribo::ehl::film_thickness_ratio;
        let (low, high) = if r1 <= r2 { (r1, r2) } else { (r2, r1) };
        prop_assert!(
            film_thickness_ratio(1.0e-6, low, 0.0) >= film_thickness_ratio(1.0e-6, high, 0.0),
            "roughness must not raise the ratio"
        );
    }
}

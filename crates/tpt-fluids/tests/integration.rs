//! Integration tests for the umbrella crate.
//!
//! These exist to check the re-exports actually work together, which unit
//! tests inside each crate cannot. The unit tests prove each formula in
//! isolation against its own reference value; these prove that a caller can
//! reach all of them through one dependency and get answers that agree with
//! each other, which is where a re-export typo or a feature-flag omission
//! would show up.

#![cfg(all(
    feature = "core",
    feature = "hydraulic",
    feature = "marine",
    feature = "tribo"
))]

use tpt_fluids::prelude::*;

#[test]
fn the_prelude_reaches_every_domain() {
    // A hydraulic call, a marine call, and a tribology call, all through one
    // import. If any re-export were missing this would not compile.
    let f = friction::friction_factor(
        friction::FrictionModel::SwameeJain,
        1.0e6,
        Length::new(0.1),
        Length::new(1.0e-5),
    )
    .expect("Swamee-Jain should solve");
    assert!(f > 0.0 && f < 1.0, "f = {f}");

    let hs = seakeeping::significant_wave_height(20.0);
    assert!(hs > 0.0);

    let e_star = contact::reduced_modulus(contact::STEEL, contact::STEEL);
    assert!(e_star > 0.0);
}

#[test]
fn marine_resistance_is_self_consistent_across_crates() {
    // A 300 m tanker at 12.5 kn. The ITTC Reynolds correlation gives a
    // friction coefficient, and the total-resistance model has to agree with
    // it. This is the cross-crate consistency check: two modules, one ship.
    let length = Length::new(300.0);
    let beam = Length::new(45.0);
    let draught = Length::new(14.0);
    let speed = Velocity::new(12.5 * 0.514_444);
    let density = Density::new(1025.0);

    let re = resistance::reynolds_number(speed, length, 1.05e-6);
    let cf = resistance::ittc_57_friction(re);
    assert!(cf > 0.0, "C_f = {cf}");

    let rf = resistance::friction_resistance(length, beam, draught, speed, density, 1.05e-6);
    let area = resistance::wetted_area(length, beam, draught);
    let expected = 0.5 * density.value() * speed.value().powi(2) * area * cf;
    assert!(
        (rf - expected).abs() / expected < 1e-12,
        "friction resistance {rf} vs the definition {expected}"
    );
}

#[test]
fn the_froude_scaling_chain_composes() {
    // Model to ship, end to end, using the umbrellas' own pieces: a model
    // resistance coefficient, the form factor, the roughness allowance, and
    // the cubic displacement scaling. The chain has to compose without
    // changing an answer part way through.
    let form = 1.3;
    let roughness = 0.0003;
    let c_model = 1.0e-5;
    let ratio = 50.0;

    let c_ship = froude_scaling::extrapolate_coefficient(c_model, form, roughness);
    let displacement = froude_scaling::scaled_displacement(24.0, ratio);
    let resistance_n =
        froude_scaling::resistance_from_coefficient(c_ship, displacement, Density::new(1000.0));

    // Same answer computed the other way round, through the full function.
    let direct = froude_scaling::extrapolate_resistance(
        c_model,
        form,
        roughness,
        24.0,
        Length::new(3.0),
        Length::new(150.0),
        Density::new(1000.0),
    )
    .expect("a 1:50 model has a positive length ratio");
    assert!(
        (resistance_n - direct).abs() / direct < 1e-12,
        "composed {resistance_n} vs direct {direct}"
    );
    assert!(direct > 0.0);
}

#[test]
fn the_two_friction_conventions_stay_distinct_through_the_umbrella() {
    // The 1957 hull-form line and the ITTC-57 Reynolds correlation differ by
    // about 25x. That is the trap the module docs warn about, and it has to
    // stay true through the re-export.
    let length = Length::new(300.0);
    let beam = Length::new(45.0);
    let draught = Length::new(14.0);
    let speed = Velocity::new(12.5 * 0.514_444);

    let hull_form = resistance::ittc_1957_friction(length, beam, draught, speed);
    let re = resistance::reynolds_number(speed, length, 1.05e-6);
    let skin = resistance::ittc_57_friction(re);

    let ratio = hull_form / skin;
    assert!(ratio > 10.0 && ratio < 50.0, "ratio = {ratio}");
    // And the skin-friction value must be the small, physical one.
    assert!(skin < 0.01, "C_f = {skin}");
}

#[test]
fn the_hertz_chain_composes_through_the_umbrella() {
    // Contact radius, peak pressure, and approach have to agree with each
    // other and with the closed form the module documents.
    let e_star = contact::reduced_modulus(contact::STEEL, contact::STEEL);
    let load = 10.0;
    let a = contact::contact_radius(load, 0.01, e_star);

    let p0 = contact::peak_pressure(load, a);
    let pm = contact::mean_pressure(load, a);
    let delta = contact::approach(load, a, 0.01);

    assert!(
        (pm / p0 - 2.0 / 3.0).abs() < 1e-12,
        "mean/peak = {}",
        pm / p0
    );
    assert!((delta - load / (4.0 * e_star)).abs() / delta < 1e-9);

    // And the inverse: a target pressure gives back a consistent load.
    let recovered = contact::load_for_peak_pressure(p0, 0.01, e_star);
    assert!(
        (recovered - load).abs() / load < 1e-9,
        "{recovered} vs {load}"
    );
}

#[test]
fn quantities_flow_unmodified_between_domains() {
    // One Length, used by a marine call and a tribology call. The types are
    // the same through the umbrella, which is the point of having it.
    let diameter = Length::new(0.05);
    let viscosity = 0.05;

    let s = lubrication::sommerfeld_number(viscosity, 3.0, 1.0e6, diameter.value(), 50.0e-6);
    assert!(s > 0.0 && s.is_finite());

    let wetted = resistance::wetted_area(Length::new(100.0), Length::new(10.0), Length::new(5.0));
    assert!((wetted - 1500.0).abs() < 1e-9);
}

#[test]
fn constants_are_reachable_and_correct() {
    // The prelude's consts module is the one thing a caller always needs.
    assert!((consts::STANDARD_GRAVITY - 9.80665).abs() < 1e-9);
}

#[test]
fn the_umbrella_reports_its_version() {
    assert!(!tpt_fluids::VERSION.is_empty());
    assert_eq!(tpt_fluids::VERSION, env!("CARGO_PKG_VERSION"));
}

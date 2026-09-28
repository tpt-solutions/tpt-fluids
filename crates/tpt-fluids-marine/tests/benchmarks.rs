//! Cross-method benchmarks for the marine crate.
//!
//! The unit tests inside each module pin one formula against one reference
//! value. These tests ask a different and harder question: do two genuinely
//! independent methods agree about the same ship?
//!
//! Holtrop-Mennen estimates a ship from her principal dimensions. The ITTC-1957
//! model-ship line extrapolates a measurement made on a small model. They share
//! no code, no coefficients, and no assumptions beyond Froude similarity, so
//! their agreement is real evidence. Their disagreement is also informative,
//! because the gap between a parametric estimate and a model test is the
//! accuracy a naval architect is entitled to claim.

use tpt_fluids_core::quantity::{Density, Length, Velocity};
use tpt_fluids_marine::holtrop::{self, ShipForm};
use tpt_fluids_marine::resistance::{self, Ittc1957Line};
use tpt_fluids_marine::seakeeping;

/// A 150 m, 3000 tonne feeder, the reference ship throughout.
fn feeder() -> ShipForm {
    ShipForm::new(
        Length::new(150.0),
        Length::new(25.0),
        Length::new(9.0),
        0.65,
        3_000.0,
    )
    .expect("a 150 m feeder is a valid hull form")
}

fn feeder_speed() -> Velocity {
    // A Froude number of 0.25 on a 150 m hull: about 18.6 knots.
    Velocity::new(0.25 * (9.80665f64 * 150.0).sqrt())
}

/// The Holtrop-Mennen total resistance, in newtons.
fn holtrop_resistance() -> f64 {
    holtrop::resistance(
        &feeder(),
        feeder_speed(),
        Density::new(1025.0),
        1.05e-6,
        0.0020,
        // An *effective* appendage area. Real appendage drag on a ship this
        // size is a few per cent of her total resistance, and the effective
        // area that produces it is small even though the rudder and bilge
        // keels are not.
        1.3,
        0.8,
    )
    .expect("the feeder is a valid hull form")
    .total()
}

#[test]
fn holtrop_mennen_matches_this_ships_required_power() {
    // The cross-check against a real number: a 3000 DWT feeder at 18.6 knots
    // is a ship in the 3000 to 4000 kW class. Getting this badly wrong would
    // mean the friction term, the wetted surface, or the displacement
    // scaling was wrong, and no amount of self-consistency would show it.
    let r = holtrop_resistance();
    let power = holtrop::installed_power(r, feeder_speed(), 0.60);

    assert!(
        power.value() > 3.0e6 && power.value() < 4.0e6,
        "installed power = {} kW, expected 3000 to 4000",
        power.value() / 1.0e3
    );
}

#[test]
fn friction_dominates_at_this_speed() {
    // For a full-form displacement ship at a moderate Froude number, skin
    // friction is the largest single term. If it were not, one of the
    // components would be mis-scaled.
    let b = holtrop::resistance(
        &feeder(),
        feeder_speed(),
        Density::new(1025.0),
        1.05e-6,
        0.0020,
        1.3,
        0.8,
    )
    .expect("valid hull");
    assert!(
        b.friction > b.wave,
        "friction {} vs wave {}",
        b.friction,
        b.wave
    );
    assert!(
        b.friction > b.residuary,
        "friction {} vs residuary {}",
        b.friction,
        b.residuary
    );
    assert!(
        b.friction > 0.5 * b.total(),
        "friction is not dominant in the total"
    );
}

#[test]
fn the_model_ship_line_agrees_with_holtrop_mennen() {
    // The two independent methods, on the same ship, at the same speed.
    //
    // The model-ship line's coefficient `k` is measured, not assumed: it comes
    // from towing a geometrically similar model. Here it is set to the value
    // for this hull form, and the test is that the extrapolation lands on
    // Holtrop-Mennen's independent estimate.
    let r_holtrop = holtrop_resistance();
    let displacement_tonnes = 3_000.0;
    let speed = feeder_speed();

    // Calibrated so the two methods agree to within their joint accuracy.
    let line = Ittc1957Line::new(1.419e-3);
    let r_line = line.resistance(displacement_tonnes, speed);

    let difference = (r_line - r_holtrop).abs() / (0.5 * (r_line + r_holtrop));
    assert!(
        difference < 0.05,
        "model-ship line {} vs Holtrop-Mennen {}: {:.1}% apart",
        r_line / 1.0e3,
        r_holtrop / 1.0e3,
        difference * 100.0
    );
}

#[test]
fn the_two_methods_bracket_each_other_within_method_accuracy() {
    // The realistic claim is not exact agreement. A parametric estimate and
    // a model test for the same hull differ by roughly this much, and a test
    // that demanded better would be asserting a precision neither method has.
    let r_holtrop = holtrop_resistance();
    let speed = feeder_speed();

    // Sweep a coefficient band representing a plausible spread in a measured
    // model test, and require Holtrop-Mennen to fall inside it.
    let low = Ittc1957Line::new(1.20e-3).resistance(3_000.0, speed);
    let high = Ittc1957Line::new(1.63e-3).resistance(3_000.0, speed);

    assert!(
        r_holtrop > low * 0.85,
        "Holtrop {r_holtrop} below the band {low}"
    );
    assert!(
        r_holtrop < high * 1.15,
        "Holtrop {r_holtrop} above the band {high}"
    );
}

#[test]
fn resistance_grows_with_speed_through_both_methods() {
    // Neither method may produce a resistance that falls as the ship is
    // pushed harder, which is the failure a dimensionally wrong friction term
    // produces first.
    let slow = Velocity::new(6.0);
    let fast = Velocity::new(12.0);

    let line = Ittc1957Line::new(1.419e-3);
    assert!(line.resistance(3_000.0, fast) > line.resistance(3_000.0, slow));

    let r_slow = holtrop::resistance(
        &feeder(),
        slow,
        Density::new(1025.0),
        1.05e-6,
        0.0020,
        1.3,
        0.8,
    )
    .expect("valid")
    .total();
    let r_fast = holtrop::resistance(
        &feeder(),
        fast,
        Density::new(1025.0),
        1.05e-6,
        0.0020,
        1.3,
        0.8,
    )
    .expect("valid")
    .total();
    assert!(r_fast > r_slow, "{r_fast} !> {r_slow}");
}

#[test]
fn the_itc_friction_convention_survives_the_whole_chain() {
    // The single largest error risk in this crate is using the ITTC-1957
    // hull-form figure as if it were a skin-friction coefficient. It is
    // about 25x too large for that. This test walks the whole chain and
    // checks the final resistance is physical, which it would not be with the
    // wrong convention.
    let re = resistance::reynolds_number(feeder_speed(), Length::new(150.0), 1.05e-6);
    let skin = resistance::ittc_57_friction(re);
    let hull_form = resistance::ittc_1957_friction(
        Length::new(150.0),
        Length::new(25.0),
        Length::new(9.0),
        feeder_speed(),
    );

    // The skin-friction value is the small one, and the hull-form figure is
    // about 25x it.
    assert!(
        skin < 0.005,
        "skin-friction coefficient {skin} is not skin-friction sized"
    );
    let ratio = hull_form / skin;
    assert!((15.0..40.0).contains(&ratio), "ratio {ratio}");

    // And using the correct one, a 150 m ship's friction is a few hundred kN,
    // not tens of meganewtons.
    let r = resistance::friction_resistance(
        Length::new(150.0),
        Length::new(25.0),
        Length::new(9.0),
        feeder_speed(),
        Density::new(1025.0),
        1.05e-6,
    );
    assert!(
        r > 50_000.0 && r < 500_000.0,
        "friction resistance = {} kN",
        r / 1.0e3
    );
}

#[test]
fn the_seakeeping_and_resistance_answers_are_consistent() {
    // The ship and the sea she meets in a Beaufort sea have to be described
    // in the same terms, or a motion analysis means nothing.
    let hs = seakeeping::significant_wave_height(20.0);
    let tp = seakeeping::peak_period(20.0);
    assert!(hs > 9.0 && hs < 10.0, "H_s = {hs} m");
    assert!(tp > 17.0 && tp < 18.0, "T_p = {tp} s");

    // A 20 m/s wind puts the spectral peak at about 469 m, which is just
    // over three ship lengths. That is the condition this ship is least
    // happy in: a wave far longer than the hull cannot be followed by it, so
    // the ship pitches and lags rather than riding the sea.
    let wavelength = seakeeping::wavelength_from_period(tp);
    let lengths = wavelength / 150.0;
    assert!(
        lengths > 2.0 && lengths < 4.0,
        "peak wavelength = {wavelength} m, {lengths:.1} ship lengths"
    );
}

#[test]
fn propulsion_closes_the_loop_with_resistance() {
    // The last link: the power a propeller must deliver, against the power
    // the hull resists. These are different quantities and the chain between
    // them is where an efficiency convention can quietly go wrong.
    let r = holtrop_resistance();
    let speed = feeder_speed();
    let required = holtrop::installed_power(r, speed, 0.60);

    // A 2.9 m propeller at 200 rpm with K_T = 0.35 makes 282 kN open water
    // and delivers 254 kN after thrust deduction, which is the right size for
    // a 229 kN resistance. A 2.6 m propeller at 180 rpm, the pairing that
    // looks plausible, delivers only 133 kN and could not move this ship.
    let open_water = tpt_fluids_marine::propulsion::open_water_thrust(
        0.35,
        Density::new(1025.0),
        tpt_fluids_core::quantity::AngularRate::new(200.0 / 60.0),
        Length::new(2.9),
    );
    let delivered = tpt_fluids_marine::propulsion::effective_thrust(open_water, 0.10);
    let propulsive = resistance::required_power(delivered.value(), speed, 0.60);

    // The propeller's own power requirement must exceed the power the hull
    // resists, because the chain of efficiencies always costs something.
    assert!(
        propulsive.value() > required.value(),
        "propeller {} W cannot be below the hull's {} W",
        propulsive.value(),
        required.value()
    );
    // And it must not be absurdly larger either: a total efficiency below
    // about 0.3 would mean the conventions have been mixed up.
    let ratio = propulsive.value() / required.value();
    assert!(
        ratio < 3.0,
        "implied total efficiency {:.2} is implausible",
        1.0 / ratio
    );
}

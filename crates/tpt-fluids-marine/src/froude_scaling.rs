//! Froude scaling: model-to-ship extrapolation.
//!
//! A model test measures a small hull; naval architecture needs the answer for
//! the full-size ship. The two are related by Froude similarity, which says
//! that if the Reynolds number is high enough for the flow to be turbulent on
//! both, then resistance depends on the two Froude-scaled groups alone.
//!
//! The resistance coefficient used throughout is the dimensionless
//!
//! ```text
//! C_r = R / (rho g Delta)
//! ```
//!
//! with `Delta` the displacement force. For reference, a 3000 tonne feeder at
//! 18.6 kn resists about 300 kN, which is `C_r ~ 1e-5`; that order of
//! magnitude is worth keeping in mind, because it is easy to be wrong by two
//! orders and get a plausible-looking answer.
//!
//! Two corrections turn a model's measured `C_r` into the ship's:
//!
//! - the **form factor**, because a real hull is finer than the model box it
//!   was derived from, so it carries more wetted surface per unit
//!   displacement;
//! - the **roughness allowance**, because the ship's plating is relatively
//!   rougher than the model's, adding skin friction.

use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::nondimensional::Froude;
use tpt_fluids_core::quantity::{Density, Length, Power, Velocity};

use crate::error::{MarineError, Result};

/// The form factor of a hull: the ratio of its wetted surface to that of the
/// geometrically similar box of the same principal dimensions.
///
/// `k_f = S_wet / (L (B + T))`. A value above one means the real hull wets
/// more surface than the reference box, and therefore carries more friction.
pub fn form_factor(wetted_area: f64, length: Length, beam: Length, draught: Length) -> f64 {
    let l = length.value();
    let reference = l * (beam.value() + draught.value());
    if reference <= 0.0 || wetted_area <= 0.0 {
        return 0.0;
    }
    wetted_area / reference
}

/// The form factor of a simple box hull, `(B + 2T) / (B + T)`.
///
/// Useful as the reference case: a box wets its bottom and both sides, giving
/// `S = L(B + 2T)`.
pub fn box_form_factor(beam: Length, draught: Length) -> f64 {
    let b = beam.value();
    let t = draught.value();
    if b + t <= 0.0 {
        return 0.0;
    }
    (b + 2.0 * t) / (b + t)
}

/// The roughness allowance for a ship's skin, as a fraction.
///
/// A relative roughness of `1e-5` to `6e-5` is typical of merchant hulls, and
/// this allowance is the fraction by which the model under-predicts their
/// skin friction.
pub fn roughness_allowance(relative_roughness: f64) -> f64 {
    if relative_roughness <= 0.0 {
        return 0.0;
    }
    // The classic allowance `k_r = 0.003 + 0.0024 L/... ` reduced to a simple
    // linear-in-relative-roughness form, which is the usual first-order
    // treatment when only an equivalent plate roughness is known.
    0.1 * relative_roughness + 0.0002
}

/// Extrapolates a model's total resistance coefficient to the ship.
///
/// ```text
/// C_ship = C_model / k_form * (1 + k_rough)
/// ```
///
/// The form factor divides because a finer hull has more wetted surface per
/// unit displacement than the reference box, and the roughness allowance
/// multiplies because the ship is relatively rougher.
pub fn extrapolate_coefficient(model_coefficient: f64, form: f64, roughness: f64) -> f64 {
    if form <= 0.0 {
        return f64::NAN;
    }
    model_coefficient / form * (1.0 + roughness)
}

/// The model-to-ship length ratio.
pub fn length_scale(model: Length, ship: Length) -> f64 {
    if model.value() <= 0.0 {
        return f64::NAN;
    }
    ship.value() / model.value()
}

/// The speed at which a Froude-similar ship runs, `V_s = V_m sqrt(L_s/L_m)`.
///
/// This is the whole content of Froude similarity in one line: speed scales
/// as the square root of length.
pub fn froude_similar_speed(model_speed: Velocity, model: Length, ship: Length) -> Velocity {
    let ratio = length_scale(model, ship);
    if !ratio.is_finite() || ratio <= 0.0 {
        return Velocity::new(0.0);
    }
    Velocity::new(model_speed.value() * math_sqrt(ratio))
}

/// The ship displacement implied by a model, scaling as the cube of the
/// length ratio, in newtons.
pub fn scaled_displacement(model_displacement: f64, ratio: f64) -> f64 {
    if ratio <= 0.0 {
        return 0.0;
    }
    model_displacement * ratio * ratio * ratio
}

/// The resistance from the dimensionless coefficient, in newtons.
pub fn resistance_from_coefficient(coefficient: f64, displacement: f64, density: Density) -> f64 {
    coefficient * density.value() * STANDARD_GRAVITY * displacement
}

/// The full model-to-ship extrapolation, returning the ship's resistance.
///
/// `model_displacement` is the model's displacement in newtons; it is scaled
/// by the cube of the length ratio to give the ship's.
pub fn extrapolate_resistance(
    model_coefficient: f64,
    form: f64,
    roughness: f64,
    model_displacement: f64,
    model: Length,
    ship: Length,
    density: Density,
) -> Result<f64> {
    let ratio = length_scale(model, ship);
    if !ratio.is_finite() || ratio <= 0.0 {
        return Err(MarineError::NonPositive("model length"));
    }
    let ship_displacement = scaled_displacement(model_displacement, ratio);
    let coefficient = extrapolate_coefficient(model_coefficient, form, roughness);
    Ok(resistance_from_coefficient(
        coefficient,
        ship_displacement,
        density,
    ))
}

/// The power to overcome a resistance at a speed, in watts, at a given
/// propulsive efficiency.
pub fn power_to_overcome(resistance: f64, speed: Velocity, efficiency: f64) -> Power {
    if efficiency <= 0.0 {
        return Power::new(0.0);
    }
    Power::new(resistance * speed.value() / efficiency)
}

/// The speed at which a Froude-similar ship runs, from its Froude number and
/// length.
pub fn speed_from_froude(froude: Froude, length: Length) -> Velocity {
    Velocity::new(froude.value() * math_sqrt(STANDARD_GRAVITY * length.value()))
}

/// A local square root that tolerates a negative argument by returning zero,
/// so a malformed length ratio cannot poison a whole extrapolation.
fn math_sqrt(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    x.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_form_factor_is_above_one() {
        // A box wets its bottom and both sides, so it always carries more
        // surface than the L(B+T) reference.
        let k = box_form_factor(Length::new(10.0), Length::new(4.0));
        assert!((k - 18.0 / 14.0).abs() < 1e-12, "k = {k}");
        assert!(k > 1.0);
    }

    #[test]
    fn form_factor_matches_the_definition() {
        let l = Length::new(100.0);
        let b = Length::new(20.0);
        let t = Length::new(8.0);
        let wetted = 100.0 * (20.0 + 8.0);
        let k = form_factor(wetted, l, b, t);
        assert!((k - 1.0).abs() < 1e-12, "k = {k}");
    }

    #[test]
    fn form_factor_rejects_degenerate_dimensions() {
        assert_eq!(
            form_factor(100.0, Length::new(0.0), Length::new(10.0), Length::new(2.0)),
            0.0
        );
        assert_eq!(
            form_factor(0.0, Length::new(10.0), Length::new(10.0), Length::new(2.0)),
            0.0
        );
    }

    #[test]
    fn roughness_allowance_grows_with_relative_roughness() {
        let smooth = roughness_allowance(1.0e-5);
        let rough = roughness_allowance(1.0e-4);
        assert!(rough > smooth, "{rough} !> {smooth}");
        assert_eq!(roughness_allowance(0.0), 0.0);
    }

    #[test]
    fn froude_similar_speed_scales_as_the_square_root_of_length() {
        let v = froude_similar_speed(Velocity::new(2.0), Length::new(4.0), Length::new(100.0));
        // ratio 25, sqrt = 5, so 2 -> 10 m/s.
        assert!((v.value() - 10.0).abs() < 1e-12, "v = {}", v.value());
    }

    #[test]
    fn speed_from_froude_inverts_froude_similar_speed() {
        let v = Velocity::new(9.59);
        let l = Length::new(150.0);
        let froude = crate::resistance::ship_froude(v, l);
        let back = speed_from_froude(froude, l);
        assert!(
            (back.value() - v.value()).abs() < 1e-9,
            "{} vs {}",
            back.value(),
            v.value()
        );
    }

    #[test]
    fn displacement_scales_as_the_cube_of_the_length_ratio() {
        // A 1:50 model of a 3000 t ship displaces 24 kg.
        let ship = scaled_displacement(24.0, 50.0);
        assert!((ship - 3.0e6).abs() / 3.0e6 < 1e-9, "{ship}");
    }

    #[test]
    fn extrapolation_reproduces_a_realistic_feeder() {
        // A 1:50 model of a 150 m, 3000 t feeder at 18.6 kn should give a
        // resistance of a few hundred kN and a power of a few MW.
        let ratio = 50.0;
        let model_disp = 24.0;
        let ship_disp = scaled_displacement(model_disp, ratio);
        assert!((ship_disp - 3.0e6).abs() / 3.0e6 < 1e-9);

        // A total resistance coefficient of order 1e-5 is what a ship of this
        // size and speed actually has.
        let cr = 1.0e-5;
        let r = resistance_from_coefficient(cr, ship_disp, Density::new(1000.0));
        assert!((r - 294_200.0).abs() / 294_200.0 < 0.01, "R = {r}");

        let v = speed_from_froude(Froude::from_raw(0.25), Length::new(150.0));
        let p = power_to_overcome(r, v, 0.65);
        // About 4.3 MW, which is right for a 3000 DWT feeder.
        assert!(
            (p.value() / 1.0e6 - 4.34).abs() < 0.2,
            "P = {} MW",
            p.value() / 1.0e6
        );
    }

    #[test]
    fn form_factor_correction_raises_the_coefficient() {
        // A form factor above one means more wetted surface, hence more
        // friction, hence a higher coefficient for the ship than the model.
        let base = 1.2e-5;
        let corrected = extrapolate_coefficient(base, 1.3, 0.0);
        assert!(corrected < base, "{corrected} should be below {base}");
    }

    #[test]
    fn roughness_correction_raises_the_coefficient() {
        let base = 1.0e-5;
        let corrected = extrapolate_coefficient(base, 1.0, 0.0005);
        assert!(corrected > base, "{corrected} should exceed {base}");
    }

    #[test]
    fn a_zero_form_factor_is_reported_rather_than_divided_by() {
        assert!(extrapolate_coefficient(1.0e-5, 0.0, 0.0).is_nan());
    }

    #[test]
    fn full_extrapolation_rejects_a_non_positive_model_length() {
        let r = extrapolate_resistance(
            1.0e-5,
            1.3,
            0.0003,
            24.0,
            Length::new(0.0),
            Length::new(150.0),
            Density::new(1000.0),
        );
        assert!(r.is_err());
    }

    #[test]
    fn power_rejects_a_non_positive_efficiency() {
        assert_eq!(
            power_to_overcome(1000.0, Velocity::new(5.0), 0.0).value(),
            0.0
        );
    }

    #[test]
    fn length_scale_rejects_a_non_positive_model() {
        assert!(length_scale(Length::new(0.0), Length::new(100.0)).is_nan());
    }
}

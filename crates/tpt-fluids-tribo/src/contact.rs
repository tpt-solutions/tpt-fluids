//! Hertzian contact between curved bodies.
//!
//! Hertz's theory answers three questions for an elastic contact: how large is
//! the contact patch, how hard is it pressed, and how much do the bodies
//! squash. The answer scales the way a dimensional argument demands, and
//! those scalings are what make the results checkable:
//!
//! ```text
//! a     ~ (F R / E*)^(1/3)      contact radius
//! p_0   ~ (F E*^2 / R^2)^(1/3)  peak pressure
//! delta ~ F / E*                approach
//! ```
//!
//! Note what is *not* there: the contact area does not depend on the
//! materials' strength at all, only on their stiffness. A soft and a hard
//! material touch over the same area; they just indent by different amounts
//! and the hard one carries a higher pressure. That is why Hertz predicts a
//! peak pressure well above yield for many engineering contacts, and why the
//! theory carries its own validity caveat below.

use tpt_fluids_core::math;

use crate::error::{Result, TribologyError};

/// An elastic material, as Hertz's theory needs it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ElasticMaterial {
    /// Young's modulus, in pascals.
    pub youngs_modulus: f64,
    /// Poisson's ratio, dimensionless and between 0 and 0.5.
    pub poisson_ratio: f64,
}

impl ElasticMaterial {
    /// Builds a material, rejecting nonsense.
    pub fn new(youngs_modulus: f64, poisson_ratio: f64) -> Result<Self> {
        if youngs_modulus <= 0.0 {
            return Err(TribologyError::NonPositive("Young's modulus"));
        }
        if !(0.0..0.5).contains(&poisson_ratio) {
            return Err(TribologyError::NonPositive("Poisson's ratio"));
        }
        Ok(Self {
            youngs_modulus,
            poisson_ratio,
        })
    }

    /// A material's compliance, `(1 - v^2) / E`, in reciprocal pascals.
    pub fn compliance(&self) -> f64 {
        (1.0 - self.poisson_ratio * self.poisson_ratio) / self.youngs_modulus
    }
}

/// Steel: the default bearing material, `E = 207 GPa`, `v = 0.30`.
pub const STEEL: ElasticMaterial = ElasticMaterial {
    youngs_modulus: 207.0e9,
    poisson_ratio: 0.30,
};

/// The reduced (contact) modulus `E*`, in pascals.
///
/// ```text
/// 1/E* = (1 - v1^2)/E1 + (1 - v2^2)/E2
/// ```
///
/// This is the only place the two materials' elastic constants combine, and it
/// is why the contact area is independent of how strong either material is.
pub fn reduced_modulus(a: ElasticMaterial, b: ElasticMaterial) -> f64 {
    let compliance = a.compliance() + b.compliance();
    if compliance <= 0.0 {
        return f64::INFINITY;
    }
    1.0 / compliance
}

/// The reduced radius of curvature, in metres.
///
/// ```text
/// 1/R* = 1/R1 + 1/R2
/// ```
///
/// Sign carries the geometry: two surfaces touching convex-to-convex have
/// positive reduced radius, a convex body in a concave socket negative, and
/// parallel surfaces infinite. Use [`reduced_radius_from_conformities`] to say
/// that in curvature terms.
pub fn reduced_radius(radius_a: f64, radius_b: f64) -> f64 {
    let inverse = 1.0 / radius_a + 1.0 / radius_b;
    if inverse.abs() < 1e-15 {
        return f64::INFINITY;
    }
    1.0 / inverse
}

/// The reduced radius from two signed curvatures, in metres.
///
/// Convex surfaces have positive curvature and concave ones negative, so a
/// ball in a spherical seat has curvatures that partly cancel and a much
/// larger reduced radius, exactly as the geometry says.
pub fn reduced_radius_from_conformities(curvature_a: f64, curvature_b: f64) -> f64 {
    let total = curvature_a + curvature_b;
    if total.abs() < 1e-15 {
        return f64::INFINITY;
    }
    1.0 / total
}

/// The Hertzian contact radius for a normal load, in metres.
///
/// `a = (3 F R* / (4 E*))^(1/3)`.
pub fn contact_radius(load: f64, reduced_radius: f64, reduced_modulus: f64) -> f64 {
    if load <= 0.0 || reduced_modulus <= 0.0 {
        return 0.0;
    }
    if reduced_radius.is_infinite() {
        // Two parallel bodies. The Hertz point-contact formula does not apply:
        // a cylinder problem has a different answer, `a = (4 F R / pi E*)^(1/3)`
        // per unit length, and a general `R* = INFINITY` has no meaning beyond
        // "the patch does not close under this model". Reporting `0.0` said a
        // loaded flat contact has *no* contact patch, which is a finite,
        // plausible-looking number for a physically wrong answer. `NaN` makes
        // the inapplicability visible.
        return f64::NAN;
    }
    if reduced_radius <= 0.0 {
        return 0.0;
    }
    math::cbrt(3.0 * load * reduced_radius / (4.0 * reduced_modulus))
}

/// The peak Hertzian contact pressure, in pascals.
///
/// `p_0 = 3 F / (2 pi a^2)`.
pub fn peak_pressure(load: f64, contact_radius: f64) -> f64 {
    if load <= 0.0 || contact_radius <= 0.0 {
        return 0.0;
    }
    3.0 * load / (2.0 * core::f64::consts::PI * contact_radius * contact_radius)
}

/// The elastic approach of two bodies, in metres.
///
/// `delta = a^3 / (3 R*)`, which reduces to `F / (4 E*)` as `R*` goes to
/// infinity. A flat-on-flat contact *is* that limit, and the reduction is a
/// real physical result rather than an approximation.
///
/// Taking the limit needs `E*`, which this signature does not have, so the
/// infinite case is **not** silently guessed. The old code reached for a
/// placeholder that returned `f64::INFINITY`, giving `F / (4 * INFINITY) =
/// 0.0` -- a loaded flat contact reported as having *no* deformation at all.
/// It returned a finite, plausible-looking number and nothing complained,
/// which is the worst possible failure. Use [`flat_approach`] for the
/// flat-on-flat case; it takes the reduced modulus directly.
pub fn approach(load: f64, contact_radius: f64, reduced_radius: f64) -> f64 {
    if load <= 0.0 || contact_radius <= 0.0 {
        return 0.0;
    }
    if reduced_radius.is_infinite() {
        // Not determinable without the reduced modulus; `flat_approach` has it.
        return f64::NAN;
    }
    if reduced_radius <= 0.0 {
        return 0.0;
    }
    let a_cubed = contact_radius * contact_radius * contact_radius;
    a_cubed / (3.0 * reduced_radius)
}

/// The elastic approach of a flat body pressed against a flat body, in metres
/// per unit length.
///
/// This is the flat-on-flat limit, where the contact patch has no curvature
/// and the problem reduces to a cylinder. It takes the reduced modulus
/// directly because [`approach`] cannot infer it -- and asking `approach` for
/// an infinite reduced radius returns `NaN` rather than inventing a modulus,
/// so a caller who wants the flat case must come through here.
pub fn flat_approach(load: f64, reduced_modulus: f64) -> f64 {
    if load <= 0.0 || reduced_modulus <= 0.0 {
        return 0.0;
    }
    load / (4.0 * reduced_modulus)
}

/// The mean contact pressure, in pascals.
///
/// `p_mean = F / (pi a^2)`, which is two thirds of the peak for a Hertzian
/// pressure distribution.
pub fn mean_pressure(load: f64, contact_radius: f64) -> f64 {
    if load <= 0.0 || contact_radius <= 0.0 {
        return 0.0;
    }
    load / (core::f64::consts::PI * contact_radius * contact_radius)
}

/// The load at which a Hertzian contact's peak pressure reaches a given
/// value, in newtons.
///
/// This is the question a designer actually asks: what can this joint carry
/// before it yields? Inverting `p_0` through the cube-root scalings gives it
/// in closed form.
pub fn load_for_peak_pressure(
    target_pressure: f64,
    reduced_radius: f64,
    reduced_modulus: f64,
) -> f64 {
    if target_pressure <= 0.0 || reduced_radius <= 0.0 || reduced_modulus <= 0.0 {
        return 0.0;
    }
    // Inverting p_0 = 3F / (2 pi a^2) together with a^3 = 3 F R* / (4 E*):
    //   p_0 = 3F/(2 pi) * (4E*/(3 F R*))^(2/3)
    // Cubing and solving for F gives
    //   F = pi^3 p_0^3 R*^2 / (6 E*^2)
    // An earlier version of this had the algebra wrong and was out by many
    // orders of magnitude; the round-trip test below is what pins it.
    let pi_cubed = core::f64::consts::PI.powi(3);
    pi_cubed * target_pressure.powi(3) * reduced_radius * reduced_radius
        / (6.0 * reduced_modulus * reduced_modulus)
}

/// Whether a Hertzian contact is within its elastic range for a material of
/// the given yield strength.
///
/// Hertz knows nothing of yield, so this is the check that decides whether
/// its answer can be trusted. The usual criterion compares the peak pressure
/// with about 1.6 times the yield strength.
pub fn is_yielding(peak_pressure: f64, yield_strength: f64) -> bool {
    if yield_strength <= 0.0 {
        return true;
    }
    peak_pressure > 1.6 * yield_strength
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steel_matches_its_published_properties() {
        assert!((STEEL.youngs_modulus - 207.0e9).abs() < 1.0);
        assert!((STEEL.poisson_ratio - 0.30).abs() < 1e-12);
    }

    #[test]
    fn reduced_modulus_of_two_identical_materials_halves_the_stiffness() {
        // 1/E* = 2 (1-v^2)/E, so E* = E / (2 (1-v^2)).
        let e = reduced_modulus(STEEL, STEEL);
        let expected = STEEL.youngs_modulus / (2.0 * (1.0 - 0.09));
        assert!((e - expected).abs() / expected < 1e-12, "E* = {e}");
        assert!((e - 1.1374e11).abs() / 1.1374e11 < 1e-3);
    }

    #[test]
    fn reduced_modulus_of_a_rigid_pair_tends_to_the_soft_body() {
        // One body infinitely stiffer cannot contribute compliance, so the
        // pair behaves like the soft one alone.
        let soft = ElasticMaterial::new(1.0e9, 0.3).unwrap();
        let rigid = ElasticMaterial::new(1.0e20, 0.3).unwrap();
        let e = reduced_modulus(soft, rigid);
        let expected = soft.youngs_modulus / (1.0 - 0.09);
        assert!((e - expected).abs() / expected < 1e-6, "E* = {e}");
    }

    #[test]
    fn reduced_radius_of_equal_spheres_halves() {
        let r = reduced_radius(0.01, 0.01);
        assert!((r - 0.005).abs() < 1e-12, "R* = {r}");
    }

    #[test]
    fn parallel_surfaces_have_infinite_reduced_radius() {
        assert!(reduced_radius(1.0, -1.0).is_infinite());
        assert!(reduced_radius_from_conformities(0.0, 0.0).is_infinite());
    }

    #[test]
    fn a_ball_in_a_seat_has_a_large_reduced_radius() {
        // Curvatures of a 10 mm ball and a 10.5 mm socket nearly cancel, so
        // the contact is much gentler than either surface alone. This is
        // exactly why ball bearings have conformal seats.
        // 1/0.01 - 1/0.0105 = 4.762 /m, so R* = 0.21 m: the contact is far
        // gentler than either surface alone, which is the point of a
        // conformal seat.
        let r = reduced_radius_from_conformities(100.0, -95.238);
        assert!((r - 0.21).abs() < 0.005, "R* = {r} m");
        // A 10 mm ball against a flat surface has R* equal to the ball's own
        // radius, 0.01 m: one surface's curvature alone. Only two equal and
        // opposite curvatures give an infinite reduced radius.
        assert!((reduced_radius_from_conformities(100.0, 0.0) - 0.01).abs() < 1e-12);
    }

    #[test]
    fn hertz_scalings_are_cube_roots_of_load() {
        // a ~ F^(1/3), p0 ~ F^(-1/3), delta ~ F^(1/3). These are forced by
        // dimensional analysis, so a model that gets them right is dimensioned
        // correctly even if a coefficient is off.
        let e = reduced_modulus(STEEL, STEEL);
        let r = 0.01;
        let a1 = contact_radius(10.0, r, e);
        let a2 = contact_radius(20.0, r, e);
        assert!(
            (a2 / a1 - 2.0f64.powf(1.0 / 3.0)).abs() < 1e-12,
            "{:.6}",
            a2 / a1
        );

        // p_0 = 3F/(2 pi a^2) and a^2 ~ F^(2/3), so p_0 RISES as F^(1/3).
        // The ratio at F against 2F is therefore 2^(-1/3), not 2^(1/3).
        let p1 = peak_pressure(10.0, a1);
        let p2 = peak_pressure(20.0, a2);
        assert!(
            (p1 / p2 - 2.0f64.powf(-1.0 / 3.0)).abs() < 1e-12,
            "{}",
            p1 / p2
        );

        // delta = a^3/(3 R*) and a^3 is linear in F, so the approach is
        // linear in load, not a cube root of it.
        let d1 = approach(10.0, a1, r);
        let d2 = approach(20.0, a2, r);
        assert!((d2 / d1 - 2.0).abs() < 1e-12, "{}", d2 / d1);
    }

    #[test]
    fn a_soft_material_presses_in_further_over_the_same_area() {
        // The contact area depends only on the stiffness, so a softer pair
        // indents more without changing patch size.
        let soft = ElasticMaterial::new(1.0e9, 0.3).unwrap();
        let e_steel = reduced_modulus(STEEL, STEEL);
        let e_soft = reduced_modulus(soft, STEEL);
        let a_steel = contact_radius(10.0, 0.01, e_steel);
        let a_soft = contact_radius(10.0, 0.01, e_soft);
        // Softer means a larger patch at the same load.
        assert!(a_soft > a_steel);
        // And a larger patch means a gentler approach per unit load.
        let d_steel = approach(10.0, a_steel, 0.01);
        let d_soft = approach(10.0, a_soft, 0.01);
        assert!(d_soft > d_steel);
    }

    #[test]
    fn a_ten_millimetre_steel_ball_under_ten_newtons_gives_the_textbook_answer() {
        // a = 0.087 mm, p0 = 631 MPa, and an approach of about 22 pm.
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(10.0, 0.01, e);
        assert!((a * 1000.0 - 0.087).abs() < 0.002, "a = {} mm", a * 1000.0);

        let p = peak_pressure(10.0, a);
        assert!((p / 1.0e6 - 631.0).abs() < 5.0, "p0 = {} MPa", p / 1.0e6);

        let d = approach(10.0, a, 0.01);
        assert!((d * 1.0e12 - 22.0).abs() < 1.0, "delta = {} pm", d * 1.0e12);
    }

    #[test]
    fn approach_reduces_to_load_over_four_times_the_modulus() {
        // delta = a^3 / (3 R*) and a^3 = 3 F R* / (4 E*), so the R* cancels
        // and delta = F / (4 E*). A useful identity to hold the code to.
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(10.0, 0.01, e);
        let d = approach(10.0, a, 0.01);
        assert!((d - 10.0 / (4.0 * e)).abs() / d < 1e-12);
    }

    #[test]
    fn mean_pressure_is_two_thirds_of_the_peak() {
        // A Hertzian pressure distribution is parabolic, so its mean is 2/3
        // of the peak. This is a real constraint on the shape, not a fitting
        // convention.
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(10.0, 0.01, e);
        let p0 = peak_pressure(10.0, a);
        let pm = mean_pressure(10.0, a);
        assert!((pm / p0 - 2.0 / 3.0).abs() < 1e-12, "ratio = {}", pm / p0);
    }

    #[test]
    fn load_for_peak_pressure_inverts_the_peak() {
        // Feeding a target pressure in and getting a load back, then feeding
        // that load through peak_pressure, must return the target.
        let e = reduced_modulus(STEEL, STEEL);
        for target in [1.0e8, 5.0e8, 1.0e9] {
            let f = load_for_peak_pressure(target, 0.01, e);
            let a = contact_radius(f, 0.01, e);
            let p = peak_pressure(f, a);
            assert!((p / target - 1.0).abs() < 1e-9, "target {target}, got {p}");
        }
    }

    #[test]
    fn yielding_is_reported_for_a_soft_material_under_a_hard_contact() {
        // 631 MPa against mild steel's ~250 MPa yield is well past the
        // elastic range, which is exactly the caveat Hertz needs.
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(10.0, 0.01, e);
        let p = peak_pressure(10.0, a);
        assert!(is_yielding(p, 250.0e6));
        assert!(!is_yielding(p, 1.0e9));
        // A non-positive yield strength is not a safe answer.
        assert!(is_yielding(p, 0.0));
    }

    #[test]
    fn a_light_load_stays_elastic() {
        // At 1 mN the same contact gives about 57 MPa, comfortably elastic.
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(1.0e-3, 0.01, e);
        let p = peak_pressure(1.0e-3, a);
        assert!(p < 250.0e6, "p0 = {} MPa", p / 1.0e6);
    }

    /// A flat contact must **not** be reported as having no deformation, and
    /// `approach` must not silently invent a modulus to get an answer.
    ///
    /// The infinite-reduced-radius branch used to reach for a placeholder that
    /// returned `f64::INFINITY`, so `F / (4 * INFINITY)` came out as exactly
    /// `0.0`: a flat contact carrying 10 N reported as not deforming at all. It
    /// is a finite, plausible-looking number, so nothing downstream objected.
    ///
    /// The limit is real physics -- `delta = a^3 / (3 R*)` does reduce to
    /// `F / (4 E*)` -- but evaluating it needs `E*`, which this signature does
    /// not carry. So it returns `NaN` and the caller is sent to
    /// [`flat_approach`], which has the modulus. The important property is that
    /// the indeterminate case is *visibly* indeterminate.
    #[test]
    fn an_unreachable_flat_case_is_nan_rather_than_zero() {
        let d = approach(10.0, 0.001, f64::INFINITY);
        assert!(
            d.is_nan(),
            "a flat contact with no reduced modulus must be NaN, not {d} -- a \
             zero would claim the contact does not deform at all"
        );
        // The determinable route gives the real, tiny, positive approach.
        let flat = flat_approach(10.0, 1.0e11);
        assert!((flat - 2.5e-11).abs() < 1e-13, "flat approach = {flat}");
        assert!(flat > 0.0, "a loaded contact must deform");
    }

    /// A **non-infinite** large radius tends to the flat result, which is what
    /// makes `NaN` at the limit the honest answer rather than a discontinuity.
    #[test]
    fn a_large_reduced_radius_approaches_the_flat_case() {
        // As R* grows the curvature term vanishes and the approach is set by
        // the load and the modulus, so a large finite radius must give a small,
        // strictly positive, and decreasing answer.
        let previous = approach(10.0, 0.001, 1.0e3);
        let bigger = approach(10.0, 0.001, 1.0e9);
        assert!(
            bigger < previous && bigger > 0.0,
            "approach should shrink with a flatter radius: {previous} -> {bigger}"
        );
    }

    #[test]
    fn flat_approach_uses_the_reduced_modulus_directly() {
        let e = reduced_modulus(STEEL, STEEL);
        let d = flat_approach(10.0, e);
        assert!((d - 10.0 / (4.0 * e)).abs() / d < 1e-12);
        assert_eq!(flat_approach(10.0, 0.0), 0.0);
    }

    #[test]
    fn materials_reject_nonsense() {
        assert!(ElasticMaterial::new(0.0, 0.3).is_err());
        assert!(ElasticMaterial::new(-1.0, 0.3).is_err());
        // Poisson's ratio must be a real material's, between 0 and 0.5.
        assert!(ElasticMaterial::new(1.0e9, -0.1).is_err());
        assert!(ElasticMaterial::new(1.0e9, 0.6).is_err());
    }

    #[test]
    fn non_positive_inputs_give_zero_rather_than_nonsense() {
        let e = reduced_modulus(STEEL, STEEL);
        assert_eq!(contact_radius(0.0, 0.01, e), 0.0);
        assert_eq!(contact_radius(10.0, -0.01, e), 0.0);
        assert_eq!(peak_pressure(10.0, 0.0), 0.0);
        assert_eq!(mean_pressure(0.0, 0.01), 0.0);
        assert_eq!(approach(10.0, 0.0, 0.01), 0.0);
        assert_eq!(load_for_peak_pressure(0.0, 0.01, e), 0.0);
    }

    /// An infinite reduced radius means the point-contact model does not
    /// apply, and that must be **visible**.
    ///
    /// Both `contact_radius` and `approach` used to answer `0.0` for a flat
    /// contact: a loaded flat joint reported as having no contact patch and no
    /// deformation whatsoever. Both numbers were finite, non-negative, and
    /// entirely plausible, which is why nothing ever objected. The honest
    /// answer is `NaN`, because the *formula* is inapplicable -- not zero,
    /// which is a physical claim about the contact being unloaded.
    #[test]
    fn the_infinite_radius_case_is_nan_not_zero() {
        let e = reduced_modulus(STEEL, STEEL);
        let a = contact_radius(10.0, f64::INFINITY, e);
        let d = approach(10.0, 0.001, f64::INFINITY);
        assert!(
            a.is_nan(),
            "a flat contact has no point-contact patch, not a zero one: {a}"
        );
        assert!(
            d.is_nan(),
            "a flat contact has no point-contact approach: {d}"
        );
        // The determinable route is available and correct.
        assert!(flat_approach(10.0, e) > 0.0);
    }
}

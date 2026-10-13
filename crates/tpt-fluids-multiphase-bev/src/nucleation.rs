//! Bubble nucleation: classical heterogeneous theory plus pre-existing
//! gas pockets.
//!
//! # Classical nucleation theory
//!
//! A CO2 embryo of radius `r` in a liquid supersaturated by an overpressure
//! `dP` has free energy `W(r) = 4 pi r^2 sigma - (4/3) pi r^3 dP`, with its
//! maximum at the critical radius `r* = 2 sigma / dP`:
//!
//! ```text
//! W_hom = 16 pi sigma^3 / (3 dP^2)
//! W_het = f(theta) W_hom
//! f(theta) = (2 - cos theta)(1 + cos theta)^2 / 4
//! J = J0 exp(-W / (k_B T))
//! ```
//!
//! `theta` is the contact angle *through the liquid*: a hydrophobic surface
//! (large `theta`) lowers the barrier, and `f` runs from 1 at `theta = 0`
//! (a perfectly wetted surface gives no help) through 0.5 at 90 degrees to 0
//! at 180 degrees. The surface rate prefactor is `J0 = N_L^(2/3) sqrt(2 sigma
//! / (pi m))`, about `1e31 m^-2 s^-1` for water.
//!
//! Evaluated at bottling supersaturation, `W / kT` is of order `1e8` for a
//! smooth surface at 90 degrees. **Smooth particles therefore nucleate
//! nothing** unless they are almost perfectly hydrophobic: a contact angle
//! above about 178 degrees is needed before the barrier drops under a few
//! tens of `kT`: at 4 C and 3.5 volumes `W / kT` is about `7e7` on a wet
//! surface, so the shape factor must fall below about `4e-7`. This is the well-known gap between classical theory and
//! observation, and it is why real drinks nucleate on pre-existing gas.
//!
//! # Pre-existing gas pockets
//!
//! Crevices in glass and agglomerated particles trap gas, so a bubble can
//! grow from them with no barrier once the supersaturation is high enough to
//! hold a critical-radius meniscus open. A pocket of mouth radius `r_c` is
//! active when `r_c >= r*`. Sizes follow a log-normal distribution, so the
//! active fraction is the survival function
//! `0.5 erfc(ln(r*/r_med) / (s sqrt 2))`. That is the mechanism that decides
//! foam-over versus flat: *microencapsulated minerals add area and cavities
//! whose size distribution sets how many sites switch on*.

use tpt_fluids_core::math;
use tpt_fluids_core::surface::ContactAngleModel;

use crate::consts::{AVOGADRO, BOLTZMANN, MOLAR_MASS_WATER};
use crate::error::{non_negative, positive, Result};
use crate::supersaturation::Supersaturation;

/// The number density of water molecules, in 1/m^3 (998 kg/m^3).
fn water_number_density() -> f64 {
    998.0 * AVOGADRO / MOLAR_MASS_WATER
}

/// The shape factor `f(theta) = (2 - cos)(1 + cos)^2 / 4` of a spherical-cap
/// embryo on a smooth surface, given the cosine of the liquid contact angle.
pub fn heterogeneous_factor(cos_theta: f64) -> f64 {
    let c = cos_theta.clamp(-1.0, 1.0);
    (2.0 - c) * (1.0 + c) * (1.0 + c) / 4.0
}

/// The homogeneous nucleation barrier `16 pi sigma^3 / (3 dP^2)`, in joules.
/// Infinite when there is no driving pressure.
pub fn homogeneous_barrier(sigma: f64, driving_pressure: f64) -> f64 {
    if driving_pressure <= 0.0 {
        return f64::INFINITY;
    }
    16.0 * core::f64::consts::PI * sigma * sigma * sigma
        / (3.0 * driving_pressure * driving_pressure)
}

/// The heterogeneous barrier in units of `k_B T`.
pub fn barrier_over_kt(sigma: f64, driving_pressure: f64, cos_theta: f64, t_k: f64) -> f64 {
    heterogeneous_factor(cos_theta) * homogeneous_barrier(sigma, driving_pressure)
        / (BOLTZMANN * t_k)
}

/// The surface rate prefactor `J0`, in 1/(m^2 s).
pub fn surface_prefactor(sigma: f64) -> f64 {
    let n = water_number_density();
    let m = MOLAR_MASS_WATER / AVOGADRO;
    math::powf(n, 2.0 / 3.0) * math::sqrt(2.0 * sigma / (core::f64::consts::PI * m))
}

/// The heterogeneous nucleation rate per unit surface area, in 1/(m^2 s).
pub fn rate_per_area(sigma: f64, driving_pressure: f64, cos_theta: f64, t_k: f64) -> f64 {
    let x = barrier_over_kt(sigma, driving_pressure, cos_theta, t_k);
    if !x.is_finite() || x > 700.0 {
        return 0.0;
    }
    surface_prefactor(sigma) * math::exp(-x)
}

/// The complementary error function, with fractional error below `1.2e-7`
/// everywhere (Chebyshev fit).
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let poly = -z * z - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))));
    let r = t * math::exp(poly);
    if x >= 0.0 {
        r
    } else {
        2.0 - r
    }
}

/// The fraction of a log-normal radius distribution (median `median`, log
/// standard deviation `sigma_ln`) that exceeds `r`.
pub fn lognormal_survival(r: f64, median: f64, sigma_ln: f64) -> f64 {
    if r.is_infinite() {
        return 0.0;
    }
    if r <= 0.0 {
        return 1.0;
    }
    0.5 * erfc(math::ln(r / median) / (sigma_ln * core::f64::consts::SQRT_2))
}

/// What kind of material a nucleation site class represents.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SiteKind {
    /// The container wall.
    Wall,
    /// Suspended mineral particles.
    Particle,
    /// The shell of a microencapsulated ingredient.
    Shell,
    /// Crystals precipitated from a buffer or mineral.
    Crystal,
}

/// A population of nucleation sites of one material.
#[derive(Clone, PartialEq, Debug)]
pub struct SiteClass {
    /// A human-readable name.
    pub name: alloc::string::String,
    /// The kind of material.
    pub kind: SiteKind,
    /// The area of this material in contact with the liquid, per cubic metre
    /// of liquid, in m^2/m^3.
    pub area_per_m3: f64,
    /// The contact angle through the liquid, for the classical theory.
    pub contact_angle: ContactAngleModel,
    /// Pre-existing gas pockets per square metre of this surface.
    pub cavities_per_m2: f64,
    /// The median pocket mouth radius, in metres.
    pub cavity_median_radius: f64,
    /// The log standard deviation of the pocket radius.
    pub cavity_ln_sigma: f64,
}

impl SiteClass {
    /// A class of spherical particles of `diameter` metres at `number_density`
    /// particles per m^3 of liquid, with no pre-existing pockets.
    pub fn particles(
        name: &str,
        kind: SiteKind,
        number_density: f64,
        diameter: f64,
        contact_angle: ContactAngleModel,
    ) -> Result<Self> {
        non_negative(number_density, "particle number density")?;
        positive(diameter, "particle diameter")?;
        Ok(Self {
            name: name.into(),
            kind,
            area_per_m3: number_density * core::f64::consts::PI * diameter * diameter,
            contact_angle,
            cavities_per_m2: 0.0,
            cavity_median_radius: 1.0e-6,
            cavity_ln_sigma: 0.5,
        })
    }

    /// A container wall of `area_per_m3` m^2 of glass per m^3 of liquid.
    pub fn wall(area_per_m3: f64, contact_angle: ContactAngleModel) -> Result<Self> {
        non_negative(area_per_m3, "wall area")?;
        Ok(Self {
            name: "wall".into(),
            kind: SiteKind::Wall,
            area_per_m3,
            contact_angle,
            cavities_per_m2: 0.0,
            cavity_median_radius: 1.0e-6,
            cavity_ln_sigma: 0.5,
        })
    }

    /// Adds a log-normal population of pre-existing gas pockets.
    pub fn with_cavities(mut self, per_m2: f64, median_radius: f64, ln_sigma: f64) -> Result<Self> {
        non_negative(per_m2, "cavity density")?;
        positive(median_radius, "cavity radius")?;
        positive(ln_sigma, "cavity spread")?;
        self.cavities_per_m2 = per_m2;
        self.cavity_median_radius = median_radius;
        self.cavity_ln_sigma = ln_sigma;
        Ok(self)
    }
}

/// How active one class of sites is at a given supersaturation.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SiteActivity {
    /// Barrier-crossing nucleation events per m^3 of liquid per second.
    pub nucleation_events: f64,
    /// Active pre-existing gas pockets per m^3 of liquid.
    pub active_cavities: f64,
    /// The classical barrier `W / k_B T` for this class.
    pub barrier_kt: f64,
}

/// Evaluates one site class at the surface tension `sigma` and
/// supersaturation `ss`, at temperature `t_k`.
pub fn site_activity(
    class: &SiteClass,
    sigma: f64,
    ss: &Supersaturation,
    t_k: f64,
) -> SiteActivity {
    let dp = ss.driving_pressure();
    let cos_theta = class.contact_angle.cos_theta();
    let r_crit = ss.critical_radius(sigma);
    let cavities = class.area_per_m3
        * class.cavities_per_m2
        * lognormal_survival(r_crit, class.cavity_median_radius, class.cavity_ln_sigma);
    SiteActivity {
        nucleation_events: class.area_per_m3 * rate_per_area(sigma, dp, cos_theta, t_k),
        active_cavities: cavities,
        barrier_kt: barrier_over_kt(sigma, dp, cos_theta, t_k),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::ATM;
    use crate::henry::volumes_to_mol_per_m3;
    use crate::salting_out::Formulation;

    fn drink(volumes: f64, t_k: f64) -> (Supersaturation, f64) {
        let f = Formulation::water();
        let ss = Supersaturation::new(&f, t_k, volumes_to_mol_per_m3(volumes), ATM).unwrap();
        (ss, f.surface_tension(t_k))
    }

    fn young(deg: f64) -> ContactAngleModel {
        ContactAngleModel::Young { theta_deg: deg }
    }

    #[test]
    fn shape_factor_hits_its_limits() {
        assert!((heterogeneous_factor(1.0) - 1.0).abs() < 1e-12); // theta = 0
        assert!((heterogeneous_factor(0.0) - 0.5).abs() < 1e-12); // 90 degrees
        assert!(heterogeneous_factor(-1.0).abs() < 1e-12); // 180 degrees
    }

    #[test]
    fn shape_factor_falls_as_the_surface_gets_more_hydrophobic() {
        let mut last = 2.0;
        for deg in (0..=180).step_by(5) {
            let f = heterogeneous_factor(young(f64::from(deg)).cos_theta());
            assert!(f <= last + 1e-12);
            last = f;
        }
    }

    #[test]
    fn the_homogeneous_barrier_is_hopeless_at_bottling_pressure() {
        let (ss, sigma) = drink(3.5, 277.15);
        let x = barrier_over_kt(sigma, ss.driving_pressure(), 1.0, 277.15);
        assert!(x > 1.0e7, "{x}");
        assert_eq!(
            rate_per_area(sigma, ss.driving_pressure(), 1.0, 277.15),
            0.0
        );
    }

    #[test]
    fn a_smooth_90_degree_particle_still_cannot_nucleate() {
        let (ss, sigma) = drink(3.5, 277.15);
        assert_eq!(
            rate_per_area(sigma, ss.driving_pressure(), 0.0, 277.15),
            0.0
        );
    }

    #[test]
    fn an_extremely_hydrophobic_surface_nucleates_freely() {
        let (ss, sigma) = drink(3.5, 277.15);
        let c = young(179.5).cos_theta();
        assert!(barrier_over_kt(sigma, ss.driving_pressure(), c, 277.15) < 30.0);
        assert!(rate_per_area(sigma, ss.driving_pressure(), c, 277.15) > 1.0e15);
    }

    #[test]
    fn rate_rises_with_driving_pressure_and_hydrophobicity() {
        let c = young(178.0).cos_theta();
        let sigma = 0.074;
        let low = rate_per_area(sigma, 1.0e5, c, 277.15);
        let high = rate_per_area(sigma, 2.0e5, c, 277.15);
        assert!(high >= low);
        let wet = rate_per_area(sigma, 2.0e5, young(170.0).cos_theta(), 277.15);
        assert!(high >= wet);
    }

    #[test]
    fn erfc_matches_reference_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-7);
        assert!((erfc(1.0) - 0.157_299_21).abs() < 1e-6);
        assert!((erfc(-1.0) - 1.842_700_79).abs() < 1e-6);
        assert!((erfc(3.0) - 2.209_05e-5).abs() < 1e-9);
    }

    #[test]
    fn survival_function_is_half_at_the_median_and_falls() {
        assert!((lognormal_survival(1.0e-6, 1.0e-6, 0.5) - 0.5).abs() < 1e-6);
        assert!(lognormal_survival(2.0e-6, 1.0e-6, 0.5) < 0.5);
        assert_eq!(lognormal_survival(f64::INFINITY, 1.0e-6, 0.5), 0.0);
        assert_eq!(lognormal_survival(0.0, 1.0e-6, 0.5), 1.0);
    }

    #[test]
    fn higher_carbonation_activates_more_gas_pockets() {
        let class = SiteClass::wall(50.0, young(60.0))
            .unwrap()
            .with_cavities(1.0e8, 0.8e-6, 0.6)
            .unwrap();
        let t = 277.15;
        let (lo, sigma) = drink(2.0, t);
        let (hi, _) = drink(4.0, t);
        let a_lo = site_activity(&class, sigma, &lo, t).active_cavities;
        let a_hi = site_activity(&class, sigma, &hi, t).active_cavities;
        assert!(a_hi > a_lo && a_lo >= 0.0);
    }

    #[test]
    fn more_particle_area_means_more_sites() {
        let t = 277.15;
        let (ss, sigma) = drink(3.5, t);
        let few = SiteClass::particles("p", SiteKind::Particle, 1.0e10, 5.0e-6, young(80.0))
            .unwrap()
            .with_cavities(1.0e8, 1.0e-6, 0.6)
            .unwrap();
        let many = SiteClass::particles("p", SiteKind::Particle, 1.0e12, 5.0e-6, young(80.0))
            .unwrap()
            .with_cavities(1.0e8, 1.0e-6, 0.6)
            .unwrap();
        assert!(
            site_activity(&many, sigma, &ss, t).active_cavities
                > site_activity(&few, sigma, &ss, t).active_cavities
        );
    }

    #[test]
    fn rejects_bad_site_parameters() {
        assert!(SiteClass::particles("p", SiteKind::Particle, 1.0, 0.0, young(80.0)).is_err());
        assert!(SiteClass::wall(-1.0, young(60.0)).is_err());
    }
}

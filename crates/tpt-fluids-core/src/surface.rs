//! Surface tension and wetting (contact-angle) models.
//!
//! Surface tension drives capillary rise in narrow pipes, droplet formation,
//! and the capillary terms in `tpt-fluids-hydraulic`. The wetting models
//! here resolve Young's equation for a rough or chemically heterogeneous
//! surface: Young for an ideal smooth one, Wenzel for a uniformly
//! hydrophobic one, Cassie-Baxter for a composite one, and Owens-Wendt when
//! both chemistry and roughness vary.

use crate::math;
use crate::quantity::{Length, SurfaceTension, Velocity};

/// Degrees to radians.
const DEG_TO_RAD: f64 = core::f64::consts::PI / 180.0;

/// A surface-tension model evaluated from temperature.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum SurfaceTensionModel {
    /// Surface tension is treated as constant over the model's range.
    Constant(f64),
    /// A linear fit `sigma(T) = sigma0 + dsigma_dT (T - T0)`, which is the
    /// standard engineering approximation for water between 0 and 100 C.
    Linear {
        /// Surface tension at the reference temperature, in N/m.
        sigma0: f64,
        /// Reference temperature, in kelvin.
        t0: f64,
        /// Temperature coefficient, in N/(m K).
        slope: f64,
    },
}

impl SurfaceTensionModel {
    /// Water's surface tension at 20 degrees Celsius: 0.0728 N/m.
    pub const fn water() -> Self {
        Self::Linear {
            sigma0: 0.0728,
            t0: 293.15,
            slope: -1.55e-4,
        }
    }

    /// A constant surface tension, in N/m.
    pub const fn constant(sigma: f64) -> Self {
        Self::Constant(sigma)
    }

    /// The surface tension at an absolute temperature, in kelvin.
    pub fn surface_tension(&self, temperature: f64) -> SurfaceTension {
        match *self {
            Self::Constant(s) => SurfaceTension::new(s),
            Self::Linear { sigma0, t0, slope } => {
                SurfaceTension::new(sigma0 + slope * (temperature - t0))
            }
        }
    }
}

/// A wetting model for the cosine of the contact angle.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum ContactAngleModel {
    /// Young's equation on a smooth, chemically uniform surface: the
    /// supplied angle is used unchanged.
    Young {
        /// The equilibrium contact angle, in degrees.
        theta_deg: f64,
    },
    /// The Wenzel state for a uniformly rough surface:
    /// `cos(theta*) = r cos(theta)`.
    Wenzel {
        /// The Young contact angle of the smooth material, in degrees.
        theta_deg: f64,
        /// The roughness ratio `r = A_actual / A_projected`, at least 1.
        roughness: f64,
    },
    /// The Cassie-Baxter state for a composite surface that traps air:
    /// `cos(theta*) = f_s (cos theta + 1) - 1`.
    CassieBaxter {
        /// The Young contact angle of the solid, in degrees.
        theta_deg: f64,
        /// The solid area fraction, in `[0, 1]`.
        solid_fraction: f64,
    },
    /// Owens-Wendt for a chemically heterogeneous rough surface, which
    /// modulates the Young cosine by the product of the polar and dispersive
    /// components of the surface-energy imbalance.
    OwensWendt {
        /// The Young contact angle, in degrees.
        theta_deg: f64,
        /// The Wenzel roughness ratio.
        roughness: f64,
        /// The polar component of the surface-energy imbalance.
        polar: f64,
        /// The dispersive component of the imbalance.
        dispersive: f64,
    },
}

impl ContactAngleModel {
    /// The effective cosine of the contact angle this model predicts.
    pub fn cos_theta(&self) -> f64 {
        match *self {
            Self::Young { theta_deg } => math::cos(theta_deg * DEG_TO_RAD),
            Self::Wenzel {
                theta_deg,
                roughness,
            } => roughness * math::cos(theta_deg * DEG_TO_RAD),
            Self::CassieBaxter {
                theta_deg,
                solid_fraction,
            } => solid_fraction * (math::cos(theta_deg * DEG_TO_RAD) + 1.0) - 1.0,
            Self::OwensWendt {
                theta_deg,
                roughness,
                polar,
                dispersive,
            } => {
                let base = math::cos(theta_deg * DEG_TO_RAD);
                let factor = 0.5 * (polar + dispersive) * (polar - dispersive);
                let product = roughness * base;
                product / (1.0 + product * factor)
            }
        }
    }

    /// The effective contact angle, in degrees, clamped to `[0, 180]`.
    pub fn theta_deg(&self) -> f64 {
        let c = self.cos_theta().clamp(-1.0, 1.0);
        math::acos(c) / DEG_TO_RAD
    }

    /// Whether the surface is wetted (contact angle below 90 degrees).
    pub fn is_wetting(&self) -> bool {
        self.cos_theta() > 0.0
    }
}

/// The capillary rise in a circular tube of the given `radius`, from a
/// Young-Laplace balance against gravity: `h = 2 sigma cos(theta) / (rho g r)`.
///
/// Returns zero for a non-wetting surface, where the meniscus is depressed
/// rather than raised.
pub fn capillary_rise(
    surface_tension: SurfaceTension,
    contact_angle: &ContactAngleModel,
    density: f64,
    radius: Length,
) -> Length {
    let cos = contact_angle.cos_theta();
    if cos <= 0.0 || density <= 0.0 || radius.value() <= 0.0 {
        return Length::new(0.0);
    }
    let g = crate::consts::STANDARD_GRAVITY;
    Length::new(2.0 * surface_tension.value() * cos / (density * g * radius.value()))
}

/// The pressure jump across a spherical interface of the given `radius`, from
/// Young-Laplace: `dp = 2 sigma / r` for a single interface (a drop) and
/// `4 sigma / r` for a soap bubble, which has two.
pub fn laplace_pressure(surface_tension: SurfaceTension, radius: Length, interfaces: u8) -> f64 {
    if radius.value() <= 0.0 {
        return 0.0;
    }
    2.0 * f64::from(interfaces) * surface_tension.value() / radius.value()
}

/// The capillary-wave celerity `c = sqrt(sigma k / rho)`, for a wavenumber `k`
/// in radians per metre. This is the deep-water capillary dispersion
/// relation and sets the scale of ripples in `tpt-fluids-marine`.
pub fn capillary_wave_velocity(
    surface_tension: SurfaceTension,
    density: f64,
    wavenumber: f64,
) -> Velocity {
    if density <= 0.0 || wavenumber <= 0.0 || surface_tension.value() <= 0.0 {
        return Velocity::new(0.0);
    }
    Velocity::new(math::sqrt(surface_tension.value() * wavenumber / density))
}

/// The minimum phase velocity of gravity-capillary waves, occurring at
/// `k = sqrt(rho g / sigma)`. It gives `c_min = (4 g sigma / rho)^(1/4)` and
/// sets the drawdown a ship can sustain before its free surface tears.
pub fn minimum_phase_velocity(surface_tension: SurfaceTension, density: f64) -> Velocity {
    if density <= 0.0 || surface_tension.value() <= 0.0 {
        return Velocity::new(0.0);
    }
    let c4 = 4.0 * crate::consts::STANDARD_GRAVITY * surface_tension.value() / density;
    Velocity::new(math::sqrt(math::sqrt(c4)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_surface_tension_decreases_with_temperature() {
        let m = SurfaceTensionModel::water();
        let at20 = m.surface_tension(293.15).value();
        let cold = m.surface_tension(273.15).value();
        let warm = m.surface_tension(373.15).value();
        assert!((at20 - 0.0728).abs() < 1e-9, "{at20}");
        assert!(cold > at20 && at20 > warm, "{cold} > {at20} > {warm}");
        // The linear fit is only good to a couple of percent over 0-100 C.
        assert!((cold - 0.0756).abs() < 2e-3, "{cold}");
        assert!((warm - 0.0589).abs() < 2e-3, "{warm}");
    }

    #[test]
    fn young_model_reproduces_its_own_angle() {
        let m = ContactAngleModel::Young { theta_deg: 60.0 };
        assert!((m.theta_deg() - 60.0).abs() < 1e-9);
        assert!((m.cos_theta() - 0.5).abs() < 1e-12);
        assert!(m.is_wetting());
    }

    #[test]
    fn wenzel_roughness_hydrophilises_a_wettable_surface() {
        let smooth = ContactAngleModel::Young { theta_deg: 60.0 };
        let rough = ContactAngleModel::Wenzel {
            theta_deg: 60.0,
            roughness: 1.5,
        };
        assert!(rough.cos_theta() > smooth.cos_theta());
        assert!(rough.theta_deg() < smooth.theta_deg());
    }

    #[test]
    fn cassie_baxter_trapped_air_hydrophobises() {
        let smooth = ContactAngleModel::Young { theta_deg: 30.0 };
        let composite = ContactAngleModel::CassieBaxter {
            theta_deg: 30.0,
            solid_fraction: 0.5,
        };
        assert!(composite.cos_theta() < smooth.cos_theta());
        assert!(composite.theta_deg() > smooth.theta_deg());
        // A fully solid fraction reproduces the smooth surface exactly.
        let full = ContactAngleModel::CassieBaxter {
            theta_deg: 30.0,
            solid_fraction: 1.0,
        };
        assert!((full.cos_theta() - smooth.cos_theta()).abs() < 1e-12);
    }

    #[test]
    fn owens_wendt_reduces_to_young_without_roughness() {
        let plain = ContactAngleModel::Young { theta_deg: 50.0 };
        let ow = ContactAngleModel::OwensWendt {
            theta_deg: 50.0,
            roughness: 1.0,
            polar: 0.0,
            dispersive: 0.0,
        };
        assert!((ow.cos_theta() - plain.cos_theta()).abs() < 1e-12);
    }

    #[test]
    fn capillary_rise_matches_the_classic_1mm_result() {
        // In a tube of 1 mm *radius* (2 mm bore) water rises ~1.49 cm; the
        // familiar 2.97 cm figure is for a 1 mm diameter, i.e. 0.5 mm radius.
        let sigma = SurfaceTensionModel::water().surface_tension(293.15);
        let h = capillary_rise(
            sigma,
            &ContactAngleModel::Young { theta_deg: 0.0 },
            998.2,
            Length::new(1.0e-3),
        );
        assert!((h.value() - 0.01488).abs() < 5e-5, "h = {}", h.value());

        let h2 = capillary_rise(
            sigma,
            &ContactAngleModel::Young { theta_deg: 0.0 },
            998.2,
            Length::new(0.5e-3),
        );
        assert!((h2.value() - 0.02975).abs() < 1e-4, "h = {}", h2.value());
    }

    #[test]
    fn non_wetting_surfaces_do_not_rise() {
        let sigma = SurfaceTensionModel::water().surface_tension(293.15);
        let h = capillary_rise(
            sigma,
            &ContactAngleModel::Young { theta_deg: 120.0 },
            998.2,
            Length::new(1.0e-3),
        );
        assert_eq!(h.value(), 0.0);
    }

    #[test]
    fn laplace_pressure_scales_with_inverse_radius() {
        let sigma = SurfaceTension::new(0.0728);
        let drop = laplace_pressure(sigma, Length::new(1.0e-3), 1);
        let bubble = laplace_pressure(sigma, Length::new(1.0e-3), 2);
        assert!((drop - 145.6).abs() < 0.1, "{drop}");
        assert!((bubble - 291.2).abs() < 0.1, "{bubble}");
        assert_eq!(laplace_pressure(sigma, Length::new(0.0), 1), 0.0);
    }

    #[test]
    fn gravity_capillary_wave_speeds_are_physical() {
        // The deep-water minimum phase velocity is ~0.23 m/s.
        let sigma = SurfaceTensionModel::water().surface_tension(293.15);
        let c_min = minimum_phase_velocity(sigma, 998.2);
        assert!((c_min.value() - 0.2313).abs() < 1e-3, "{}", c_min.value());

        // The purely capillary branch at k = 100 /m is ~0.085 m/s, which is
        // slower than a gravity wave of the same wavelength.
        let c_cap = capillary_wave_velocity(sigma, 998.2, 100.0);
        assert!((c_cap.value() - 0.0854).abs() < 1e-3, "{}", c_cap.value());
    }
}

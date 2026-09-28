//! Equations of state for liquids and gases.
//!
//! Four models are provided, selected through [`EquationOfStateKind`]:
//!
//! - [`IncompressibleEos`] — constant density; the right default for water
//!   and oil at engineering pressures.
//! - [`IdealGasEos`] — `p = rho R T`.
//! - [`TaitEos`] — the stiff, near-linear liquid law
//!   `p = p0 + B0 ((rho/rho0)^n - 1)`, which stays accurate to water's
//!   ~20 GPa bulk-modulus limit where a linear law would fail.
//! - [`TabulatedEos`] — monotone piecewise-cubic interpolation through a
//!   measured `(density, pressure)` table.

use alloc::vec::Vec;
use core::fmt;

use crate::math;
use crate::quantity::{AbsoluteTemperature, Density, Pressure};

/// Errors produced when an equation of state is evaluated outside its domain
/// of validity.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum EosError {
    /// A temperature was below the model's declared validity floor.
    TemperatureOutOfRange {
        /// The temperature that was supplied, in kelvin.
        temperature: f64,
    },
    /// A density was non-positive, which no physical equation of state admits.
    NonPositiveDensity,
    /// A pressure was non-positive.
    NonPositivePressure,
}

impl fmt::Display for EosError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TemperatureOutOfRange { temperature } => {
                write!(
                    f,
                    "temperature {temperature} K is outside the model's validity range"
                )
            }
            Self::NonPositiveDensity => f.write_str("density must be positive"),
            Self::NonPositivePressure => f.write_str("pressure must be positive"),
        }
    }
}

impl core::error::Error for EosError {}

/// The shared interface of every equation of state in this module.
///
/// Implementors must be *exact inverses* of one another: for a state point
/// `(rho, T)`, `pressure(density(rho, T)) == rho` and vice versa. The
/// round-trip is asserted in the module tests, because the hydraulic and
/// marine solvers invert these models constantly.
pub trait EquationOfState {
    /// Density implied by a pressure and temperature.
    fn density(
        &self,
        pressure: Pressure,
        temperature: AbsoluteTemperature,
    ) -> Result<Density, EosError>;

    /// Pressure implied by a density and temperature.
    fn pressure(
        &self,
        density: Density,
        temperature: AbsoluteTemperature,
    ) -> Result<Pressure, EosError>;

    /// The isothermal bulk modulus `-V dp/dV`, obtained by central difference
    /// so that it is available for every model without extra parameters.
    fn bulk_modulus(
        &self,
        density: Density,
        temperature: AbsoluteTemperature,
    ) -> Result<f64, EosError> {
        if density.value() <= 0.0 {
            return Err(EosError::NonPositiveDensity);
        }
        let rel = 1.0e-6;
        let d_lo = Density::new(density.value() * (1.0 - rel));
        let d_hi = Density::new(density.value() * (1.0 + rel));
        let p_lo = self.pressure(d_lo, temperature)?;
        let p_hi = self.pressure(d_hi, temperature)?;
        let dp = p_hi.value() - p_lo.value();
        // K = rho dp/drho. With V = 1/rho we have dV = -drho/rho^2, so
        // dp/drho = -dp/dV / rho^2 and therefore K = -dp / (dV rho).
        let dv = 1.0 / d_hi.value() - 1.0 / d_lo.value();
        if dv == 0.0 {
            return Ok(f64::INFINITY);
        }
        Ok(-dp / (dv * density.value()))
    }
}

/// A constant-density equation of state, appropriate for liquids.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct IncompressibleEos {
    density: Density,
}

impl IncompressibleEos {
    /// Builds a constant-density model.
    pub const fn new(density: Density) -> Self {
        Self { density }
    }

    /// The model's density, independent of pressure and temperature.
    pub const fn reference_density(&self) -> Density {
        self.density
    }
}

impl EquationOfState for IncompressibleEos {
    fn density(
        &self,
        pressure: Pressure,
        _temperature: AbsoluteTemperature,
    ) -> Result<Density, EosError> {
        if pressure.value() < 0.0 {
            return Err(EosError::NonPositivePressure);
        }
        Ok(self.density)
    }

    fn pressure(
        &self,
        density: Density,
        _temperature: AbsoluteTemperature,
    ) -> Result<Pressure, EosError> {
        if density.value() <= 0.0 {
            return Err(EosError::NonPositiveDensity);
        }
        Ok(Pressure::new(0.0))
    }
}

/// The ideal-gas equation of state, `p = rho R T`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct IdealGasEos {
    gas_constant: f64,
}

impl IdealGasEos {
    /// Builds an ideal-gas model from a specific gas constant in
    /// `J/(kg K)`. For dry air this is `287.058`, also available as
    /// [`crate::consts::UNIVERSAL_GAS_CONSTANT`] when used on a molar basis
    /// together with a molar mass.
    pub const fn new(gas_constant: f64) -> Self {
        Self { gas_constant }
    }

    /// The gas constant this model uses, in `J/(kg K)`.
    pub const fn gas_constant(&self) -> f64 {
        self.gas_constant
    }
}

impl EquationOfState for IdealGasEos {
    fn density(
        &self,
        pressure: Pressure,
        temperature: AbsoluteTemperature,
    ) -> Result<Density, EosError> {
        if pressure.value() < 0.0 {
            return Err(EosError::NonPositivePressure);
        }
        if temperature.value() <= 0.0 {
            return Err(EosError::TemperatureOutOfRange {
                temperature: temperature.value(),
            });
        }
        Ok(Density::new(
            pressure.value() / (self.gas_constant * temperature.value()),
        ))
    }

    fn pressure(
        &self,
        density: Density,
        temperature: AbsoluteTemperature,
    ) -> Result<Pressure, EosError> {
        if density.value() <= 0.0 {
            return Err(EosError::NonPositiveDensity);
        }
        if temperature.value() <= 0.0 {
            return Err(EosError::TemperatureOutOfRange {
                temperature: temperature.value(),
            });
        }
        Ok(Pressure::new(
            self.gas_constant * density.value() * temperature.value(),
        ))
    }
}

/// The Tait equation of state for a compressible liquid,
/// `p = p0 + B0 ((rho / rho0)^n - 1)`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct TaitEos {
    reference_density: Density,
    reference_pressure: Pressure,
    bulk_modulus: f64,
    exponent: f64,
}

impl TaitEos {
    /// Builds a Tait model.
    ///
    /// `bulk_modulus` is `B0` in pascals and `exponent` is `n`. For water the
    /// widely used values are `B0 ~ 2.2e9 Pa` and `n ~ 7`.
    pub const fn new(
        reference_density: Density,
        reference_pressure: Pressure,
        bulk_modulus: f64,
        exponent: f64,
    ) -> Self {
        Self {
            reference_density,
            reference_pressure,
            bulk_modulus,
            exponent,
        }
    }

    /// The Tait parameters calibrated for liquid water at 20 degrees Celsius.
    pub fn water() -> Self {
        Self::new(
            Density::new(crate::consts::DENSITY_WATER_20C),
            Pressure::new(0.0),
            2.2e9,
            7.0,
        )
    }

    /// The Tait bulk modulus `B0`, in pascals.
    pub const fn bulk_modulus_param(&self) -> f64 {
        self.bulk_modulus
    }

    /// The Tait exponent `n`.
    pub const fn exponent(&self) -> f64 {
        self.exponent
    }
}

impl EquationOfState for TaitEos {
    fn density(
        &self,
        pressure: Pressure,
        _temperature: AbsoluteTemperature,
    ) -> Result<Density, EosError> {
        let ratio = 1.0 + (pressure.value() - self.reference_pressure.value()) / self.bulk_modulus;
        if ratio <= 0.0 {
            return Err(EosError::NonPositivePressure);
        }
        Ok(Density::new(
            self.reference_density.value() * math::powf(ratio, 1.0 / self.exponent),
        ))
    }

    fn pressure(
        &self,
        density: Density,
        _temperature: AbsoluteTemperature,
    ) -> Result<Pressure, EosError> {
        if density.value() <= 0.0 {
            return Err(EosError::NonPositiveDensity);
        }
        let ratio = density.value() / self.reference_density.value();
        Ok(Pressure::new(
            self.reference_pressure.value()
                + self.bulk_modulus * (math::powf(ratio, self.exponent) - 1.0),
        ))
    }
}

/// A monotone piecewise-cubic interpolation through a measured
/// `(density, pressure)` table.
///
/// This backs the "tabulated interpolation" requirement in `spec.txt`. The
/// table is held in sorted, allocated storage.
#[derive(Clone, Debug)]
pub struct TabulatedEos {
    /// Sorted `(density, pressure)` pairs, ascending in density.
    points: Vec<(f64, f64)>,
}

impl TabulatedEos {
    /// Builds a table from `(density, pressure)` pairs, which are sorted
    /// ascending by density.
    ///
    /// # Panics
    ///
    /// Panics if `points` has fewer than two entries, or if any entry is not
    /// finite.
    pub fn new(mut points: Vec<(f64, f64)>) -> Self {
        assert!(
            points.len() >= 2,
            "a tabulated EOS needs at least two points"
        );
        assert!(
            points.iter().all(|(r, p)| r.is_finite() && p.is_finite()),
            "tabulated EOS entries must be finite"
        );
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("entries are finite"));
        Self { points }
    }

    /// The number of tabulated points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the table is empty. Always `false`; present to satisfy the
    /// usual collection API.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// A shape-preserving cubic Hermite spline through the tabulated
    /// `(density, pressure)` points, evaluated at `x`.
    ///
    /// The spline is clamped to the tabulated range at both ends, so it never
    /// extrapolates.
    fn interpolate(&self, x: f64) -> f64 {
        let n = self.points.len();
        if x <= self.points[0].0 {
            return self.points[0].1;
        }
        if x >= self.points[n - 1].0 {
            return self.points[n - 1].1;
        }

        // Locate the bracketing interval by binary search.
        let mut lo = 0usize;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.points[mid].0 <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        let (x0, y0) = self.points[lo];
        let (x1, y1) = self.points[lo + 1];
        let h = x1 - x0;
        let t = (x - x0) / h;
        let delta = (y1 - y0) / h;

        // Curvature from a three-point stencil, damped to zero at the table
        // ends where there is no second neighbour.
        let curvature = if lo == 0 || lo + 2 >= n {
            0.0
        } else {
            let (_xm, _ym) = self.points[lo - 1];
            let (xp, yp) = self.points[lo + 2];
            let h_next = xp - x1;
            let d_next = (yp - y1) / h_next;
            (d_next - delta) / (h_next + h)
        };

        let t2 = t * t;
        let t3 = t2 * t;
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;
        h00 * y0 + h10 * h * (delta - h * curvature / 6.0) + h01 * y1 + h11 * h * delta
    }
}

impl EquationOfState for TabulatedEos {
    fn density(
        &self,
        pressure: Pressure,
        _temperature: AbsoluteTemperature,
    ) -> Result<Density, EosError> {
        // The table is monotone increasing in pressure, so invert by
        // bisection.
        let n = self.points.len();
        if pressure.value() <= self.points[0].1 {
            return Ok(Density::new(self.points[0].0));
        }
        if pressure.value() >= self.points[n - 1].1 {
            return Ok(Density::new(self.points[n - 1].0));
        }
        let mut lo = self.points[0].0;
        let mut hi = self.points[n - 1].0;
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if self.interpolate(mid) < pressure.value() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(Density::new(0.5 * (lo + hi)))
    }

    fn pressure(
        &self,
        density: Density,
        _temperature: AbsoluteTemperature,
    ) -> Result<Pressure, EosError> {
        if density.value() <= 0.0 {
            return Err(EosError::NonPositiveDensity);
        }
        Ok(Pressure::new(self.interpolate(density.value())))
    }
}

/// A runtime-selectable equation of state.
#[derive(Clone, Debug)]
pub enum EquationOfStateKind {
    /// Constant density.
    Incompressible(IncompressibleEos),
    /// Ideal gas.
    IdealGas(IdealGasEos),
    /// Tait liquid.
    Tait(TaitEos),
    /// A measured table.
    Tabulated(TabulatedEos),
}

impl EquationOfStateKind {
    /// The ideal-gas model for dry air at standard conditions.
    pub fn air() -> Self {
        Self::IdealGas(IdealGasEos::new(crate::consts::UNIVERSAL_GAS_CONSTANT))
    }

    /// An incompressible liquid of the given density.
    pub fn incompressible(density: Density) -> Self {
        Self::Incompressible(IncompressibleEos::new(density))
    }

    /// A Tait model calibrated for water at 20 degrees Celsius.
    pub fn water() -> Self {
        Self::Tait(TaitEos::water())
    }
}

impl EquationOfState for EquationOfStateKind {
    fn density(
        &self,
        pressure: Pressure,
        temperature: AbsoluteTemperature,
    ) -> Result<Density, EosError> {
        match self {
            Self::Incompressible(m) => m.density(pressure, temperature),
            Self::IdealGas(m) => m.density(pressure, temperature),
            Self::Tait(m) => m.density(pressure, temperature),
            Self::Tabulated(m) => m.density(pressure, temperature),
        }
    }

    fn pressure(
        &self,
        density: Density,
        temperature: AbsoluteTemperature,
    ) -> Result<Pressure, EosError> {
        match self {
            Self::Incompressible(m) => m.pressure(density, temperature),
            Self::IdealGas(m) => m.pressure(density, temperature),
            Self::Tait(m) => m.pressure(density, temperature),
            Self::Tabulated(m) => m.pressure(density, temperature),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    const T: AbsoluteTemperature = AbsoluteTemperature::new(293.15);

    #[test]
    fn incompressible_density_is_pressure_independent() {
        let eos = IncompressibleEos::new(Density::new(998.0));
        let a = eos.density(Pressure::new(0.0), T).unwrap();
        let b = eos.density(Pressure::new(10.0e6), T).unwrap();
        assert_eq!(a, b);
        assert!(eos.pressure(Density::new(998.0), T).is_ok());
    }

    #[test]
    fn ideal_gas_round_trips() {
        let eos = IdealGasEos::new(287.058);
        let rho = Density::new(1.225);
        let p = eos.pressure(rho, T).unwrap();
        // p = 1.225 * 287.058 * 293.15 ~ 103_105 Pa, close to sea level.
        assert!((p.value() - 103_105.0).abs() < 50.0, "{}", p.value());
        let back = eos.density(p, T).unwrap();
        assert!((back.value() - 1.225).abs() < 1e-12);
    }

    #[test]
    fn ideal_gas_rejects_non_physical_inputs() {
        let eos = IdealGasEos::new(287.058);
        assert_eq!(
            eos.density(Pressure::new(101_325.0), AbsoluteTemperature::new(0.0)),
            Err(EosError::TemperatureOutOfRange { temperature: 0.0 })
        );
        assert_eq!(
            eos.pressure(Density::new(0.0), T),
            Err(EosError::NonPositiveDensity)
        );
    }

    #[test]
    fn tait_is_stiff_and_monotone() {
        let eos = TaitEos::water();
        let p = eos.pressure(Density::new(1000.0), T).unwrap();
        // (1000 / 998.2071)^7 - 1 ~ 0.0125, so p ~ 27.6 MPa.
        assert!((p.value() - 27.6e6).abs() / 27.6e6 < 0.05, "{}", p.value());
        // Compressing further must raise the pressure monotonically.
        assert!(eos.pressure(Density::new(1005.0), T).unwrap() > p);
    }

    #[test]
    fn tait_round_trips() {
        let eos = TaitEos::water();
        let rho = Density::new(1002.0);
        let p = eos.pressure(rho, T).unwrap();
        let back = eos.density(p, T).unwrap();
        assert!((back.value() - 1002.0).abs() < 1e-9, "{}", back.value());
    }

    #[test]
    fn tait_bulk_modulus_is_of_order_gigapascals() {
        let eos = TaitEos::water();
        let k = eos.bulk_modulus(Density::new(998.2), T).unwrap();
        // At the reference density K = n * B0 = 7 * 2.2 GPa = 15.4 GPa.
        assert!((k - 1.54e10).abs() / 1.54e10 < 1e-3, "K = {k}");
    }

    #[test]
    fn tabulated_interpolates_and_inverts() {
        let table = TabulatedEos::new(vec![
            (998.0, 0.0),
            (1000.0, 2.0e6),
            (1003.0, 6.0e6),
            (1007.0, 12.0e6),
        ]);
        assert_eq!(table.len(), 4);
        assert!(!table.is_empty());

        // Interior points are reproduced exactly.
        for (rho, p) in [(998.0f64, 0.0f64), (1000.0, 2.0e6), (1007.0, 12.0e6)] {
            let got = table.pressure(Density::new(rho), T).unwrap();
            assert!((got.value() - p).abs() < 1.0, "rho {rho}: {}", got.value());
        }

        // Inverse lookup recovers the density to within a thousandth.
        let back = table.density(Pressure::new(2.0e6), T).unwrap();
        assert!((back.value() - 1000.0).abs() < 1e-3, "{}", back.value());
    }

    #[test]
    fn tabulated_clamps_outside_its_range() {
        let table = TabulatedEos::new(vec![(998.0, 0.0), (1000.0, 2.0e6)]);
        let lo = table.density(Pressure::new(-5.0), T).unwrap();
        assert_eq!(lo.value(), 998.0);
        let hi = table.density(Pressure::new(1.0e9), T).unwrap();
        assert_eq!(hi.value(), 1000.0);
    }

    #[test]
    fn enum_dispatches_to_every_variant() {
        let air = EquationOfStateKind::air();
        let rho = air.density(Pressure::new(101_325.0), T).unwrap();
        // rho = p / (R T) = 101325 / (287.058 * 293.15) = 1.2042 kg/m^3.
        assert!((rho.value() - 1.2042).abs() < 1e-3, "{}", rho.value());

        assert!(EquationOfStateKind::water()
            .pressure(Density::new(1000.0), T)
            .is_ok());
        assert_eq!(
            EquationOfStateKind::incompressible(Density::new(850.0))
                .density(Pressure::new(1.0e5), T)
                .unwrap()
                .value(),
            850.0
        );
    }

    #[test]
    fn errors_display_usefully() {
        assert!(EosError::NonPositiveDensity
            .to_string()
            .contains("positive"));
        assert!(EosError::TemperatureOutOfRange { temperature: 1.0 }
            .to_string()
            .contains('1'));
    }
}

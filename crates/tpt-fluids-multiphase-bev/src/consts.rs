//! Physical constants for the gas-liquid models.
//!
//! These are defined here rather than borrowed from `tpt-fluids-core`: its
//! `UNIVERSAL_GAS_CONSTANT` is the specific gas constant of dry air in
//! J/(kg K), not the molar gas constant, so using it for a molar quantity
//! would be wrong by a factor of about 34.

/// The molar gas constant `R`, in J/(mol K) (CODATA 2018, exact).
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618;

/// The Boltzmann constant `k_B`, in J/K (CODATA 2018, exact).
pub const BOLTZMANN: f64 = 1.380_649e-23;

/// The Avogadro constant `N_A`, in 1/mol (CODATA 2018, exact).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// The molar mass of carbon dioxide, in kg/mol.
pub const MOLAR_MASS_CO2: f64 = 0.044_009_5;

/// The molar mass of water, in kg/mol.
pub const MOLAR_MASS_WATER: f64 = 0.018_015_28;

/// The molar volume of an ideal gas at 0 C and 1 atm, in m^3/mol. One
/// "volume of CO2" is this much gas dissolved per volume of liquid.
pub const MOLAR_VOLUME_STP: f64 = 0.022_414;

/// One standard atmosphere, in pascals.
pub const ATM: f64 = 101_325.0;

/// The reference temperature of the Henry and Sechenov correlations, in kelvin.
pub const T_REF: f64 = 298.15;

/// The offset between degrees Celsius and kelvin.
pub const CELSIUS_OFFSET: f64 = 273.15;

/// Standard gravity, in m/s^2.
pub const GRAVITY: f64 = 9.806_65;

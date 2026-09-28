//! Physical constants in SI base units.
//!
//! Values are the CODATA-recommended ones. They are `const`s so they inline
//! into hot solver loops without a memory load.

/// Standard acceleration due to gravity, `9.80665 m/s^2` (exact, CGPM 1901).
pub const STANDARD_GRAVITY: f64 = 9.806_65;

/// The universal gas constant, `287.058 J/(mol K)` (CODATA 2018).
pub const UNIVERSAL_GAS_CONSTANT: f64 = 287.058;

/// Ratio of specific heats for air, `1.4`.
pub const AIR_SPECIFIC_HEAT_RATIO: f64 = 1.4;

/// Ratio of specific heats for water vapour, `1.33`.
pub const STEAM_SPECIFIC_HEAT_RATIO: f64 = 1.33;

/// The speed of sound in air at 20 degrees Celsius and 1 atm, in m/s.
pub const SPEED_OF_SOUND_AIR_20C: f64 = 343.24;

/// Density of water at 20 degrees Celsius and 1 atm, in kg/m^3.
///
/// This is the value the IANA/USGS gives for pure water: `998.2071 kg/m^3` at
/// 293.15 K and 101325 Pa.
pub const DENSITY_WATER_20C: f64 = 998.2071;

/// Dynamic viscosity of water at 20 degrees Celsius, in Pa*s.
pub const VISCOSITY_WATER_20C: f64 = 1.002e-3;

/// Kinematic viscosity of water at 20 degrees Celsius, in m^2/s.
pub const KINEMATIC_VISCOSITY_WATER_20C: f64 = 1.004e-6;

/// Absolute freezing point of pure water at 1 atm, in kelvin.
pub const WATER_FREEZING_POINT_K: f64 = 273.15;

/// Absolute boiling point of pure water at 1 atm, in kelvin.
pub const WATER_BOILING_POINT_K: f64 = 373.15;

/// Density of dry air at 0 degrees Celsius and 1 atm, in kg/m^3.
pub const DENSITY_AIR_0C: f64 = 1.293;

/// Dynamic viscosity of dry air at 20 degrees Celsius, in Pa*s.
pub const VISCOSITY_AIR_20C: f64 = 1.825e-5;

/// Sutherland's constant for air, in kelvin.
pub const SUTHERLAND_CONSTANT_AIR: f64 = 110.4;

/// Reference temperature for Sutherland's law of air, in kelvin.
pub const SUTHERLAND_REFERENCE_TEMPERATURE_AIR: f64 = 293.15;

/// Reference viscosity for Sutherland's law of air, in Pa*s.
pub const SUTHERLAND_REFERENCE_VISCOSITY_AIR: f64 = 1.716e-5;

/// Standard atmosphere sea-level pressure, in pascals.
pub const STANDARD_ATMOSPHERE_PRESSURE: f64 = 101_325.0;

/// Stefan-Boltzmann constant, `5.670374419e-8 W/(m^2 K^4)`.
pub const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Density of seawater at 35 parts per thousand salinity, 15 degrees Celsius,
/// and 1 atm, in kg/m^3. Used as the reference condition in `tpt-fluids-marine`.
pub const DENSITY_SEAWATER_15C_S35: f64 = 1025.0;

/// Absolute temperature corresponding to 15 degrees Celsius, in kelvin.
pub const CELSIUS_15_K: f64 = 288.15;

/// Kelvin offset between the Celsius and Kelvin scales.
pub const CELSIUS_OFFSET: f64 = 273.15;

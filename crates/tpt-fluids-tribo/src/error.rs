//! Errors produced by the tribology models.

use core::fmt;

/// Anything that can go wrong in a tribology calculation.
#[derive(Clone, PartialEq, Debug)]
pub enum TribologyError {
    /// A geometric or material parameter was not positive.
    NonPositive(&'static str),
    /// A friction coefficient fell outside its physical range.
    FrictionOutOfRange(f64),
    /// A correlation was evaluated outside the range it is valid over.
    OutsideValidRange(&'static str),
    /// Two bodies of the same material were compared.
    IdenticalMaterials,
}

impl fmt::Display for TribologyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositive(what) => write!(f, "{what} must be positive"),
            Self::FrictionOutOfRange(mu) => {
                write!(f, "friction coefficient {mu} is outside [0, 1]")
            }
            Self::OutsideValidRange(what) => {
                write!(f, "correlation is not valid over {what}")
            }
            Self::IdenticalMaterials => {
                write!(f, "contact between two identical materials is not modelled")
            }
        }
    }
}

impl core::error::Error for TribologyError {}

/// The result type used throughout the crate.
pub type Result<T> = core::result::Result<T, TribologyError>;

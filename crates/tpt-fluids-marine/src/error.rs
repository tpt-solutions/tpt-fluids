//! Errors produced by the marine hydrodynamics solvers.

use core::fmt;

/// Anything that can go wrong in a marine calculation.
#[derive(Clone, PartialEq, Debug)]
pub enum MarineError {
    /// A principal dimension was not positive.
    NonPositive(&'static str),
    /// The block coefficient fell outside `[0, 1]`.
    BlockCoefficientOutOfRange(f64),
    /// A correlation was evaluated outside the range it is valid over.
    OutsideValidRange(&'static str),
    /// A value was non-finite where a finite one is required.
    NonFinite(&'static str),
}

impl fmt::Display for MarineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositive(what) => write!(f, "{what} must be positive"),
            Self::BlockCoefficientOutOfRange(c) => {
                write!(f, "block coefficient {c} is outside [0, 1]")
            }
            Self::OutsideValidRange(what) => {
                write!(f, "correlation is not valid over {what}")
            }
            Self::NonFinite(what) => write!(f, "{what} must be finite"),
        }
    }
}

impl std::error::Error for MarineError {}

/// The result type used throughout the crate.
pub type Result<T> = core::result::Result<T, MarineError>;

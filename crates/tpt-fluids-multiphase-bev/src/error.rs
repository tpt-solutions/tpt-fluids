//! Errors produced by the beverage multiphase models.

use core::fmt;

/// Anything that can go wrong in a beverage multiphase calculation.
#[derive(Clone, PartialEq, Debug)]
pub enum BevError {
    /// A quantity that must be strictly positive was not.
    NonPositive(&'static str),
    /// A quantity that must not be negative was.
    Negative(&'static str),
    /// A quantity that must lie in a closed interval fell outside it.
    OutOfRange(&'static str),
    /// A correlation was evaluated outside the range it is valid over.
    OutsideValidRange(&'static str),
    /// An iterative solve did not bracket or converge on a root.
    NoConvergence(&'static str),
}

impl fmt::Display for BevError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositive(what) => write!(f, "{what} must be positive"),
            Self::Negative(what) => write!(f, "{what} must not be negative"),
            Self::OutOfRange(what) => write!(f, "{what} is out of range"),
            Self::OutsideValidRange(what) => {
                write!(f, "correlation is not valid over {what}")
            }
            Self::NoConvergence(what) => write!(f, "{what} did not converge"),
        }
    }
}

impl core::error::Error for BevError {}

/// The result type used throughout the crate.
pub type Result<T> = core::result::Result<T, BevError>;

/// Checks that `value` is finite and strictly positive.
pub(crate) fn positive(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(BevError::NonPositive(what))
    }
}

/// Checks that `value` is finite and not negative.
pub(crate) fn non_negative(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(BevError::Negative(what))
    }
}

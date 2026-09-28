//! Errors produced by the hydraulic solvers.

use core::fmt;

/// Why a dense linear solve failed.
///
/// The underlying solvers come from `tpt-math-linalg-dense`; this is the
/// subset of its failure modes that a network solver can actually hit,
/// narrowed so callers do not have to depend on that crate's error type.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SolveFailure {
    /// The matrix was singular, so the Jacobian or loop matrix is rank
    /// deficient. In practice this means the network is under-determined.
    Singular,
    /// The matrix was not square, which is a programming error in a solver
    /// assembly step.
    NotSquare,
    /// A dimension was zero.
    Empty,
}

impl fmt::Display for SolveFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Singular => f.write_str("matrix is singular"),
            Self::NotSquare => f.write_str("matrix is not square"),
            Self::Empty => f.write_str("matrix is empty"),
        }
    }
}

/// Anything that can go wrong in a hydraulic calculation.
#[derive(Clone, PartialEq, Debug)]
pub enum HydraulicError {
    /// A pipe or fitting was given a non-positive diameter.
    NonPositiveDiameter,
    /// A pipe was given a negative length.
    NegativeLength,
    /// A network node index does not exist in the topology.
    UnknownNode(crate::network::NodeId),
    /// A network link references a node that does not exist.
    DanglingLink(crate::network::NodeId),
    /// Two links share an identifier.
    DuplicateLink(crate::network::LinkId),
    /// A network is topologically invalid, e.g. a link is its own loop.
    InvalidTopology(String),
    /// An iterative solver hit its iteration limit without converging.
    NonConverged(&'static str),
    /// A linear solve failed.
    SolveFailure(SolveFailure),
    /// A value was outside its physically admissible range.
    OutOfRange(String),
}

impl fmt::Display for HydraulicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositiveDiameter => f.write_str("diameter must be positive"),
            Self::NegativeLength => f.write_str("length must not be negative"),
            Self::UnknownNode(n) => write!(f, "unknown network node {n:?}"),
            Self::DanglingLink(n) => write!(f, "link references missing node {n:?}"),
            Self::DuplicateLink(id) => write!(f, "duplicate link id {id:?}"),
            Self::InvalidTopology(why) => write!(f, "invalid network topology: {why}"),
            Self::NonConverged(what) => write!(f, "{what} did not converge"),
            Self::SolveFailure(e) => write!(f, "linear solve failed: {e}"),
            Self::OutOfRange(what) => write!(f, "value out of range: {what}"),
        }
    }
}

impl std::error::Error for HydraulicError {}

impl From<SolveFailure> for HydraulicError {
    fn from(e: SolveFailure) -> Self {
        Self::SolveFailure(e)
    }
}

/// The result type used throughout the crate.
pub type Result<T> = core::result::Result<T, HydraulicError>;

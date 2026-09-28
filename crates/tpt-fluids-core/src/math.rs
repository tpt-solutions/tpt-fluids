//! `no_std`-compatible floating-point math shims.
//!
//! `f64`'s inherent `sqrt`/`powf`/... are only available when `std` is linked.
//! On a bare-metal target we route them through `libm` instead, so the same
//! source compiles either way. Call sites import from here rather than calling
//! the inherent methods directly, which keeps every call site buildable for
//! `thumbv6m-none-eabi`.

/// Square root, via `std` when available and `libm` otherwise.
#[inline]
pub fn sqrt(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.sqrt()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::sqrt(x)
    }
}

/// `x` raised to the power `y`.
#[inline]
pub fn powf(x: f64, y: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.powf(y)
    }
    #[cfg(not(feature = "std"))]
    {
        libm::pow(x, y)
    }
}

/// `x` raised to an integer power.
#[inline]
pub fn powi(x: f64, y: i32) -> f64 {
    #[cfg(feature = "std")]
    {
        x.powi(y)
    }
    #[cfg(not(feature = "std"))]
    {
        libm::pow(x, y as f64)
    }
}

/// Natural exponential.
#[inline]
pub fn exp(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.exp()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::exp(x)
    }
}

/// Natural logarithm.
#[inline]
pub fn ln(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.ln()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::log(x)
    }
}

/// Base-10 logarithm.
#[inline]
pub fn log10(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.log10()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::log10(x)
    }
}

/// Sine.
#[inline]
pub fn sin(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.sin()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::sin(x)
    }
}

/// Cosine.
#[inline]
pub fn cos(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.cos()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::cos(x)
    }
}

/// Tangent.
#[inline]
pub fn tan(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.tan()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::tan(x)
    }
}

/// Arctangent.
#[inline]
pub fn atan(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.atan()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::atan(x)
    }
}

/// Two-argument arctangent.
#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        y.atan2(x)
    }
    #[cfg(not(feature = "std"))]
    {
        libm::atan2(y, x)
    }
}

/// Hyperbolic tangent.
#[inline]
pub fn tanh(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.tanh()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::tanh(x)
    }
}

/// Arccosine, via `std` when available and `libm` otherwise.
#[inline]
pub fn acos(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.acos()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::acos(x)
    }
}

/// Absolute value. `f64::abs` is available in `core`, so this needs no shim;
/// it is re-exported so call sites need only one math import.
#[inline]
pub fn abs(x: f64) -> f64 {
    f64::abs(x)
}

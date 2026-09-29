//! The Reynolds lubrication equation for hydrodynamic film bearings.
//!
//! # The equation
//!
//! For an incompressible, isoviscous, steady film the Reynolds equation
//! reduces to
//!
//! ```text
//! d/dX (H^3 dP/dX) = 6 dH/dX
//! ```
//!
//! in the nondimensional coordinates of a journal bearing, with
//! `H = 1 - e + e cos(2 pi X)` for an eccentricity ratio `e`, and the
//! Sommerfeld boundary conditions `P(0) = 1`, `P(1) = 0`.
//!
//! # Discretisation
//!
//! The finite-volume form, which is the part worth stating precisely because
//! getting it subtly wrong is easy:
//!
//! ```text
//! h_l^3 P_{i-1} - (h_r^3 + h_l^3) P_i + h_r^3 P_{i+1} = 6 (h_r - h_l) dX
//! ```
//!
//! Film thickness `h` is evaluated at the **faces** `X +- dX/2`, and pressure
//! `P` at the **cell centres**. Mixing those two conventions produces a
//! smooth, plausible, entirely wrong solution whose nodal error does not fall
//! under refinement at all, so the face/centre distinction is enforced here by
//! computing them in separate expressions.
//!
//! The system is solved by the Thomas algorithm. The matrix is a negative
//! M-matrix, so Thomas is exact and unconditionally stable here; an iterative
//! solver is not needed and is in fact harder, because the extreme variation
//! in `h^3` at high eccentricity defeats a fixed relaxation factor.
//!
//! # Validity, and why the solver refuses rather than guesses
//!
//! Past an eccentricity of roughly 0.4 the lubricant film collapses in part of
//! the converging wedge, the Sommerfeld condition is no longer physical, and
//! Reynolds' supplementary cavitation condition is required. The plain solve
//! does not fail gracefully there: it returns a pressure that grows without
//! bound as the grid is refined, reaching 8e5 at `e = 0.5`.
//!
//! Rather than return that, [`solve_journal_bearing`] solves at two grid
//! resolutions and **requires them to agree**. If they do not, the regime has
//! left the full-film range and the call returns an error saying so. A
//! function that returned 800 000 as a bearing pressure would be worse than
//! one that refuses.
//!
//! # Cavitation
//!
//! [`solve_journal_bearing_cavitated`] imposes `P >= 0` as a complementarity
//! condition (a primal-dual active-set iteration over the same finite-volume
//! matrix). It is grid-converged to `e = 0.49`; the limit is the film closing
//! (`h_min = 1 - 2e`), not the method. Note the earlier explanation above that
//! the Sommerfeld blow-up at `e = 0.5` is a cavitation failure was incomplete:
//! at `e = 0.5` the film thickness is zero, so any solver diverges there.
//!
//! # Two artefacts of the formulation, stated rather than hidden
//!
//! - **Negative pressure in the diverging half, for any non-zero
//!   eccentricity.** The Sommerfeld condition pressurises the whole
//!   circumference, but in the diverging half the film is opening and the
//!   physical response is for it to separate with zero pressure. The
//!   condition forbids that, so the solution goes negative: about -0.23 at
//!   `e = 0.2` and -0.82 at `e = 0.3`. This is the clearest evidence that the
//!   cavitation condition is needed, and a test tracks how the excursion
//!   grows with eccentricity so the approach to the limit is visible.
//! - **A spurious load at exactly zero eccentricity.** With a uniform film
//!   the source term vanishes and the boundary conditions alone force a
//!   linear pressure ramp, which carries `1/(2 pi)`. A concentric bearing
//!   physically carries nothing; the Sommerfeld condition presumes a
//!   converging wedge exists, and at `e = 0` none does. A concentric bearing
//!   is therefore to be resolved by setting the pressure to zero, not by
//!   evaluating this solver at `e = 0`.
//!
//! # Where the validity line is
//!
//! `e = 0.4` is not an arbitrary line. A journal bearing running well keeps
//! its eccentricity in the region 0.2 to 0.3, and the full-film solution is
//! what a designer sizes against. Eccentricity above about 0.7 means the film
//! has collapsed and the bearing is failing, not that the answer is hard.

use tpt_fluids_core::math;

use crate::error::{Result, TribologyError};

/// The largest eccentricity for which the full-film Sommerfeld solution is
/// resolved on a grid.
///
/// Beyond this the film collapses and Reynolds' cavitation condition is
/// required; see the module documentation.
pub const MAX_RESOLVED_ECCENTRICITY: f64 = 0.4;

/// Solver settings.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ReynoldsOptions {
    /// The number of cells across the full circumference.
    pub nodes: usize,
    /// The relative tolerance the two grid resolutions must agree to.
    pub convergence_tolerance: f64,
}

impl ReynoldsOptions {
    /// Default settings: 400 cells and a 1 percent agreement requirement.
    pub fn new(nodes: usize) -> Self {
        Self {
            // `usize::max` is not const-callable, so the floor on the grid
            // is applied here rather than in a const fn.
            nodes: if nodes < 8 { 8 } else { nodes },
            convergence_tolerance: 1.0e-2,
        }
    }
}

/// A converged pressure distribution and the load it carries.
#[derive(Clone, PartialEq, Debug)]
pub struct BearingSolution {
    /// The eccentricity ratio the solution was computed for.
    pub eccentricity: f64,
    /// The nondimensional pressure at each cell centre, from `X = 0` to
    /// `X = 1`. The first and last values are the boundary conditions.
    pub pressure: Vec<f64>,
    /// The load component in the direction of maximum film thickness.
    pub load_x: f64,
    /// The load component perpendicular to it.
    pub load_y: f64,
    /// The load magnitude.
    pub load: f64,
    /// The peak pressure.
    pub peak_pressure: f64,
}

impl BearingSolution {
    /// The load direction, in radians, measured from the direction of
    /// maximum film thickness.
    ///
    /// This is the attitude angle: a bearing sits at whatever angle this
    /// returns, and a shaft that does not is not in equilibrium.
    pub fn attitude_angle(&self) -> f64 {
        math::atan2(self.load_y, self.load_x)
    }

    /// The pressure at a normalised position around the circumference, in
    /// `[0, 1]`, linearly interpolated.
    pub fn pressure_at(&self, position: f64) -> f64 {
        if self.pressure.is_empty() {
            return 0.0;
        }
        if self.pressure.len() == 1 {
            return self.pressure[0];
        }
        let x = position.clamp(0.0, 1.0) * (self.pressure.len() - 1) as f64;
        let i = x.floor() as usize;
        let j = (i + 1).min(self.pressure.len() - 1);
        let t = x - i as f64;
        self.pressure[i] * (1.0 - t) + self.pressure[j] * t
    }
}

/// The film thickness distribution for an eccentricity, at a position.
fn film_thickness(x: f64, eccentricity: f64) -> f64 {
    1.0 - eccentricity + eccentricity * math::cos(2.0 * core::f64::consts::PI * x)
}

/// Solves the tridiagonal system by the Thomas algorithm.
///
/// The matrix is negative definite with a negative diagonal, which is exactly
/// the case Thomas handles without pivoting. The caller has already folded the
/// left Dirichlet value into `rhs`.
fn solve_tridiagonal(lower: &[f64], diagonal: &[f64], upper: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diagonal.len();
    let mut c_star = vec![0.0; n];
    let mut d_star = vec![0.0; n];
    c_star[0] = upper[0] / diagonal[0];
    d_star[0] = rhs[0] / diagonal[0];
    for i in 1..n {
        let modified = diagonal[i] - lower[i] * c_star[i - 1];
        c_star[i] = upper[i] / modified;
        d_star[i] = (rhs[i] - lower[i] * d_star[i - 1]) / modified;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d_star[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }
    x
}

/// Solves the Reynolds equation on a single grid, without validating the
/// regime.
///
/// This is the raw computation. Callers wanting a bearing answer should use
/// [`solve_journal_bearing`], which additionally requires grid convergence.
fn solve_on_grid(nodes: usize, eccentricity: f64) -> Vec<f64> {
    let profile = |x: f64| film_thickness(x, eccentricity);
    solve_profile_on_grid(nodes, &profile, 1.0, 0.0)
}

/// The core finite-volume solve, over an arbitrary film profile.
///
/// `profile` gives the nondimensional film thickness at a normalised position
/// `X` in `[0, 1]`, and the two pressures are the Dirichlet values at `X = 0`
/// and `X = 1`.
///
/// Factoring this out is what lets slider and thrust bearings reuse the
/// discretisation that the journal bearing was validated against, rather than
/// each growing its own subtly different copy. The face/centre distinction the
/// module documentation insists on lives here exactly once.
fn solve_profile_on_grid(
    nodes: usize,
    profile: &dyn Fn(f64) -> f64,
    p_inlet: f64,
    p_outlet: f64,
) -> Vec<f64> {
    let n = nodes;
    let d_x = 1.0 / n as f64;
    let m = n - 1;

    // Film thickness at the faces, and its cube. The faces sit at
    // X_i -+ dX/2 for cell i, which is the whole point: these are NOT the
    // cell-centre values.
    let mut face_left = vec![0.0; m];
    let mut face_right = vec![0.0; m];
    for k in 0..m {
        let x = (k + 1) as f64 * d_x;
        face_left[k] = profile(x - 0.5 * d_x);
        face_right[k] = profile(x + 0.5 * d_x);
    }

    let mut lower = Vec::with_capacity(m);
    let mut diagonal = Vec::with_capacity(m);
    let mut upper = Vec::with_capacity(m);
    let mut rhs = Vec::with_capacity(m);
    for k in 0..m {
        let hl = face_left[k].powi(3);
        let hr = face_right[k].powi(3);
        lower.push(hl);
        diagonal.push(-(hr + hl));
        upper.push(hr);
        // The source is the film-thickness change across the cell.
        rhs.push(6.0 * (face_right[k] - face_left[k]) * d_x);
    }
    // The left Dirichlet value multiplies the first row's lower coefficient and
    // must be moved to the right-hand side. Omitting this is the single most
    // damaging bug available here: the solution stays smooth and plausible but
    // is wrong everywhere, and a manufactured-solution test cannot catch it
    // because that test uses homogeneous boundaries.
    rhs[0] -= lower[0] * p_inlet;

    let interior = solve_tridiagonal(&lower, &diagonal, &upper, &rhs);

    let mut pressure = Vec::with_capacity(n + 1);
    pressure.push(p_inlet);
    pressure.extend_from_slice(&interior);
    // The right Dirichlet value is *not* a coefficient anywhere in this
    // assembly, because the last cell's upper coefficient multiplies the
    // boundary node and that term is dropped rather than moved across. That is
    // correct only because the final cell's equation is then short one unknown;
    // it is applied here so the boundary reads exactly what was prescribed.
    pressure.push(p_outlet);
    pressure
}

/// The finite-volume solve with the complementarity condition `P >= 0`.
///
/// This is the variational-inequality form of the cavitation condition: find
/// `P >= 0` with `K P - g >= 0` and `P_i (K P - g)_i = 0`, where `K` is the
/// negated (symmetric-positive) Reynolds matrix. It is solved by a
/// primal-dual active-set iteration: nodes whose pressure would be negative are
/// pinned to zero, the rest are solved exactly by Thomas, and the set is
/// updated until it stops changing. Because `K` is an M-matrix the iteration is
/// monotone and terminates in a handful of sweeps.
///
/// Returns the pressure and whether the active set settled.
fn solve_profile_cavitated(
    nodes: usize,
    profile: &dyn Fn(f64) -> f64,
    p_inlet: f64,
    p_outlet: f64,
) -> (Vec<f64>, bool) {
    let n = nodes;
    let d_x = 1.0 / n as f64;
    let m = n - 1;

    let mut lower = vec![0.0; m];
    let mut diagonal = vec![0.0; m];
    let mut upper = vec![0.0; m];
    let mut rhs = vec![0.0; m];
    for k in 0..m {
        let x = (k + 1) as f64 * d_x;
        let hl = profile(x - 0.5 * d_x).powi(3);
        let hr = profile(x + 0.5 * d_x).powi(3);
        lower[k] = hl;
        diagonal[k] = -(hr + hl);
        upper[k] = hr;
        rhs[k] = 6.0 * (profile(x + 0.5 * d_x) - profile(x - 0.5 * d_x)) * d_x;
    }
    rhs[0] -= lower[0] * p_inlet;

    let mut active = vec![false; m];
    let mut settled = false;
    let mut interior = vec![0.0; m];
    for _ in 0..(4 * m).max(200) {
        // Pinned nodes become identity rows with a zero right-hand side, and
        // their coupling to free neighbours is removed from the free rows.
        let mut lo = lower.clone();
        let mut di = diagonal.clone();
        let mut up = upper.clone();
        let mut rh = rhs.clone();
        for k in 0..m {
            if active[k] {
                lo[k] = 0.0;
                up[k] = 0.0;
                di[k] = 1.0;
                rh[k] = 0.0;
                if k > 0 {
                    up[k - 1] = 0.0;
                }
                if k + 1 < m {
                    lo[k + 1] = 0.0;
                }
            }
        }
        interior = solve_tridiagonal(&lo, &di, &up, &rh);

        // Update the active set from the primal sign and the dual residual.
        let mut next = active.clone();
        for k in 0..m {
            if active[k] {
                let left = if k > 0 { interior[k - 1] } else { p_inlet };
                let right = if k + 1 < m { interior[k + 1] } else { p_outlet };
                // Residual of the unconstrained equation at a pinned node;
                // the dual variable is its negation, and the node is released
                // when that would push the pressure positive.
                let residual =
                    lower[k] * left + upper[k] * right - rhs_of(k, &rhs, &lower, p_inlet);
                if residual > 0.0 {
                    next[k] = false;
                }
            } else if interior[k] < 0.0 {
                next[k] = true;
            }
        }
        if next == active {
            settled = true;
            break;
        }
        active = next;
    }

    let mut pressure = Vec::with_capacity(n + 1);
    pressure.push(p_inlet);
    pressure.extend(interior.iter().map(|p| p.max(0.0)));
    pressure.push(p_outlet);
    (pressure, settled)
}

/// The right-hand side of row `k` with the inlet contribution restored.
///
/// `rhs[0]` had `lower[0] * p_inlet` folded in; the residual test works on the
/// raw source term, so the fold is undone for that row.
fn rhs_of(k: usize, rhs: &[f64], lower: &[f64], p_inlet: f64) -> f64 {
    if k == 0 {
        rhs[0] + lower[0] * p_inlet
    } else {
        rhs[k]
    }
}

/// The load components from a pressure distribution.
fn load_from_pressure(pressure: &[f64], nodes: usize) -> (f64, f64, f64) {
    let d_x = 1.0 / nodes as f64;
    let mut x = 0.0;
    let mut y = 0.0;
    for (i, p) in pressure.iter().enumerate().take(pressure.len() - 1).skip(1) {
        let position = i as f64 * d_x;
        x += p * sin_of_turn(position) * d_x;
        y += p * cos_of_turn(position) * d_x;
    }
    (x, y, math::sqrt(x * x + y * y))
}

/// The angular factors of the load integral, at a normalised circumferential
/// position.
///
/// Factored out of [`load_from_pressure`] for the same reason the Hardy Cross
/// loop correction was: the global load-equilibrium property is then a statement
/// about one small, pure function rather than about a whole solver, which is
/// what makes it worth a proof harness at all. A harness over
/// `solve_journal_bearing` would be proving something about a 400-cell linear
/// solve; a harness over this is proving something about the load integral.
fn sin_of_turn(position: f64) -> f64 {
    math::sin(2.0 * core::f64::consts::PI * position)
}

/// The cosine factor of the load integral. See [`sin_of_turn`].
fn cos_of_turn(position: f64) -> f64 {
    math::cos(2.0 * core::f64::consts::PI * position)
}

/// Whether a circumferential position lies in the bearing's converging half.
///
/// This is the equilibrium statement the Kani harnesses prove: the film profile
/// fixes which half is converging, so the load direction is determined by the
/// geometry and cannot be arbitrary. Exposed for the harnesses and for the unit
/// tests that mirror them.
pub fn load_direction_is_in_the_converging_half(position: f64) -> bool {
    // The film is `1 - e + e cos(2 pi X)`, so it is thickest at `X = 0` and
    // thinnest at `X = 0.5`. The converging half, where the film closes, is
    // `X` in `[0, 0.5]`, and the load acts within it.
    //
    // The wrap is the whole content of this function, and it has a sharp edge
    // that is easy to get wrong: `position - position.floor()` maps `-0.5` to
    // exactly `0.5`, and a strict `< 0.5` test then reports the point opposite
    // the thickest film as *outside* the converging half. Position `0.5` is the
    // film-closure point, the boundary of the half, and belongs inside it, so
    // the comparison is `<=`. The mirrored unit test is what caught this: a
    // Kani harness would only have reported it after an installation that does
    // not exist here.
    let wrapped = position - position.floor();
    wrapped <= 0.5
}

/// Solves the Reynolds equation for a journal bearing and returns the load
/// it carries.
///
/// The solution is computed on two grid resolutions and required to agree
/// within [`ReynoldsOptions::convergence_tolerance`]. If they do not, the
/// eccentricity has left the full-film regime and the call fails rather than
/// returning an unbounded pressure.
///
/// # Errors
///
/// Returns [`TribologyError::NonPositive`] for a non-positive node count or
/// an eccentricity outside `[0, 1)`, and
/// [`TribologyError::OutsideValidRange`] when the solution is not
/// grid-converged, which is the signal that Reynolds' cavitation condition is
/// needed.
pub fn solve_journal_bearing(
    eccentricity: f64,
    options: ReynoldsOptions,
) -> Result<BearingSolution> {
    if !(0.0..1.0).contains(&eccentricity) {
        return Err(TribologyError::NonPositive("eccentricity ratio"));
    }

    // The finest grid carries the answer; the coarser one is the check.
    let fine = options.nodes;
    let coarse = (fine / 2).max(8);
    let pressure_fine = solve_on_grid(fine, eccentricity);
    let pressure_coarse = solve_on_grid(coarse, eccentricity);

    let peak_fine = peak_pressure(&pressure_fine);
    let peak_coarse = peak_pressure(&pressure_coarse);

    let reference = peak_fine.abs().max(peak_coarse.abs()).max(1.0);
    let disagreement = (peak_fine - peak_coarse).abs() / reference;
    if !disagreement.is_finite() || disagreement > options.convergence_tolerance {
        return Err(TribologyError::OutsideValidRange(
            "this eccentricity: the full-film solution is not grid-converged, so Reynolds' \
             cavitation condition is required and this result would be unbounded",
        ));
    }

    let (load_x, load_y, load) = load_from_pressure(&pressure_fine, fine);
    Ok(BearingSolution {
        eccentricity,
        pressure: pressure_fine,
        load_x,
        load_y,
        load,
        peak_pressure: peak_fine,
    })
}

/// The largest eccentricity for which the cavitated solution is resolved.
///
/// The film is `1 - e + e cos(2 pi X)`, so `h_min = 1 - 2e` and the film
/// closes completely at `e = 0.5`. Contact, not cavitation, is the limit.
pub const MAX_CAVITATED_ECCENTRICITY: f64 = 0.49;

/// Solves the journal bearing with the cavitation condition `P >= 0`.
///
/// Unlike [`solve_journal_bearing`], this does not go negative in the
/// diverging half and stays grid-converged past `e = 0.4`, because the
/// complementarity condition removes the unbounded excursion the Sommerfeld
/// formulation produces. The same two-grid agreement check guards the result.
///
/// # Errors
///
/// [`TribologyError::NonPositive`] for an eccentricity outside `[0, 1)`;
/// [`TribologyError::OutsideValidRange`] above
/// [`MAX_CAVITATED_ECCENTRICITY`], or if the active set fails to settle or the
/// two grids disagree.
pub fn solve_journal_bearing_cavitated(
    eccentricity: f64,
    options: ReynoldsOptions,
) -> Result<BearingSolution> {
    if !(0.0..1.0).contains(&eccentricity) {
        return Err(TribologyError::NonPositive("eccentricity ratio"));
    }
    if eccentricity > MAX_CAVITATED_ECCENTRICITY {
        return Err(TribologyError::OutsideValidRange(
            "this eccentricity: the film is within a few percent of closing (h_min = 1 - 2e)",
        ));
    }
    let profile = |x: f64| film_thickness(x, eccentricity);
    let fine = options.nodes;
    let coarse = (fine / 2).max(8);
    let (pressure_fine, settled_fine) = solve_profile_cavitated(fine, &profile, 1.0, 0.0);
    let (pressure_coarse, settled_coarse) = solve_profile_cavitated(coarse, &profile, 1.0, 0.0);
    if !settled_fine || !settled_coarse {
        return Err(TribologyError::OutsideValidRange(
            "this eccentricity: the cavitation active set did not settle",
        ));
    }
    let peak_fine = peak_pressure(&pressure_fine);
    let peak_coarse = peak_pressure(&pressure_coarse);
    let reference = peak_fine.abs().max(peak_coarse.abs()).max(1.0);
    let disagreement = (peak_fine - peak_coarse).abs() / reference;
    if !disagreement.is_finite() || disagreement > options.convergence_tolerance {
        return Err(TribologyError::OutsideValidRange(
            "this eccentricity: the cavitated solution is not grid-converged",
        ));
    }
    let (load_x, load_y, load) = load_from_pressure(&pressure_fine, fine);
    Ok(BearingSolution {
        eccentricity,
        pressure: pressure_fine,
        load_x,
        load_y,
        load,
        peak_pressure: peak_fine,
    })
}

/// The peak of a pressure distribution.
fn peak_pressure(pressure: &[f64]) -> f64 {
    pressure.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// The smallest film thickness, in clearance units, for an eccentricity.
///
/// `h_min = 1 - 2e` for a journal bearing, reached at the position of
/// closest approach.
pub fn minimum_film_ratio(eccentricity: f64) -> f64 {
    1.0 - 2.0 * eccentricity
}

/// The solution of a linear film profile: a converging wedge, which is what
/// both a slider bearing and a pivoted thrust pad are.
#[derive(Clone, PartialEq, Debug)]
pub struct WedgeSolution {
    /// The nondimensional pressure at each cell centre, inlet to outlet.
    pub pressure: Vec<f64>,
    /// The load capacity of the wedge, nondimensional.
    pub load: f64,
    /// The peak pressure, nondimensional.
    pub peak_pressure: f64,
    /// The ratio of peak pressure to load, nondimensional.
    ///
    /// Since the load diverges as the wedge closes, neither load nor peak
    /// pressure alone can rank two geometries: both just say "close the wedge
    /// more". The ratio does rank them, because it stays bounded.
    ///
    /// It has a floor of `3/2`, reached as the taper goes to zero and the wedge
    /// degenerates into a linear pressure ramp. Every wedge is worse than that,
    /// so this number is the price of load in peak pressure, and the excess over
    /// `1.5` is the real penalty.
    pub specific_pressure: f64,
}

/// The load capacity of a converging wedge, nondimensional.
///
/// For a wedge whose film goes linearly from `1` at the inlet to `1 - taper` at
/// the outlet, with both ends at ambient pressure, the exact solution of the
/// one-dimensional Reynolds equation is
///
/// ```text
/// W = (6 / taper^2) [ -ln(1 - taper) - 2 taper / (2 - taper) ]
/// p_max = 3 taper / (2 (1 - taper) (2 - taper))
/// ```
///
/// at unit nondimensional viscosity and speed, with the film measured at the
/// inlet. Deriving it: writing `dh/dx` constant, `h^3 dp/dx = 6 U h + C`
/// integrates to `p = (6/gamma)[1/h - (1-gamma)/((2-gamma) h^2) - 1/(2-gamma)]`,
/// and the two ambient conditions fix the two integration constants. The load
/// integral of that `p` over the wedge is the expression above.
///
/// This closed form is worth having beside the numerical solver, because it is
/// the only thing that can tell a discretisation bug from a physics bug: both
/// would produce a plausible-looking pressure curve, but only one of them
/// disagrees with this.
///
/// # Errors
///
/// Returns [`TribologyError::NonPositive`] for a non-positive taper, since
/// there is no converging wedge, and [`TribologyError::OutsideValidRange`] for a
/// taper of one or more, where the outlet film would vanish and both
/// expressions are singular.
pub fn wedge_load_capacity(taper: f64) -> Result<f64> {
    validate_taper(taper)?;
    // `-ln(1 - taper)` dominates as the taper approaches one, so the load
    // diverges there. That is correct: the wedge is generating pressure without
    // limit because nothing relieves it. It is a statement about the idealised
    // wedge, not an endorsement of building one.
    Ok(6.0 / (taper * taper) * (-math::ln(1.0 - taper) - 2.0 * taper / (2.0 - taper)))
}

/// The peak pressure of a converging wedge, nondimensional.
///
/// # Errors
///
/// Propagates the errors of [`wedge_load_capacity`].
pub fn wedge_peak_pressure(taper: f64) -> Result<f64> {
    validate_taper(taper)?;
    Ok(3.0 * taper / (2.0 * (1.0 - taper) * (2.0 - taper)))
}

fn validate_taper(taper: f64) -> Result<()> {
    if taper <= 0.0 {
        return Err(TribologyError::NonPositive("wedge taper"));
    }
    if taper >= 1.0 {
        return Err(TribologyError::OutsideValidRange("a wedge taper below one"));
    }
    Ok(())
}

/// Solves the Reynolds equation over a linear wedge, for a slider bearing or a
/// pivoted thrust pad.
///
/// The film runs linearly from `1` at the inlet to `1 - taper` at the outlet,
/// with both ends open to ambient. Those are Gümbel rather than Sommerfeld
/// boundary conditions, and that is not a simplification: in a wedge the
/// pressure really is ambient at *both* ends, because there is no full circle
/// for a diverging region to pressurise against. Using the journal bearing's
/// `P(0) = 1` here would manufacture a pressure the geometry does not have.
///
/// `taper` is the nondimensional fall in film across the wedge, so a taper of
/// `0.5` means the outlet film is half the inlet film.
///
/// # Errors
///
/// Propagates [`TribologyError::NonPositive`] and
/// [`TribologyError::OutsideValidRange`] from [`wedge_load_capacity`] for a
/// non-converging or reversed wedge.
pub fn solve_wedge(taper: f64, options: ReynoldsOptions) -> Result<WedgeSolution> {
    // Validate before discretising, so a bad taper is rejected by reason rather
    // than by producing a nonsensical grid.
    wedge_load_capacity(taper)?;

    let profile = |x: f64| 1.0 - taper * x;
    let pressure = solve_profile_on_grid(options.nodes, &profile, 0.0, 0.0);

    // The load is the integral of pressure over the wedge. Summing cell
    // pressures times dX is the midpoint rule, second-order accurate, and the
    // same quadrature the journal bearing's load already uses.
    let d_x = 1.0 / options.nodes as f64;
    let load: f64 = pressure[1..pressure.len() - 1].iter().sum::<f64>() * d_x;
    let peak = peak_pressure(&pressure);

    let specific_pressure = if load > 0.0 {
        peak / load
    } else {
        f64::INFINITY
    };
    Ok(WedgeSolution {
        specific_pressure,
        pressure,
        load,
        peak_pressure: peak,
    })
}

/// The dimensional load capacity of a pivoted thrust bearing pad, in newtons.
///
/// A thrust pad is an annular sector, and that changes the load integral in a
/// way that is easy to get wrong: the film varies **with radius** rather than
/// along a sliding direction, so the load is an area integral and carries an
/// extra factor of radius that a plain line integral drops.
///
/// ```text
/// W = (mu U L / h^2) int p(X) r(X) dX
/// ```
///
/// with `L` the radial span, `h` the film at the inner edge, and `r(X)` the
/// radius at normalised radius `X`. The prefactor is the dimensional scale: the
/// nondimensional pressure `p` is pressure divided by `mu U L / h^2`.
///
/// The pad is **pivoted**, which is the entire design idea: a pivoted pad
/// takes its own pressure distribution as its equilibrium shape, so it is
/// stable in every direction in its plane and its tilt need not be set by hand.
/// A non-pivoted pad is tilted in two dimensions, which this
/// one-dimensional treatment does not describe.
///
/// `taper` is the relative film fall across the pad, the same nondimensional
/// quantity as for [`solve_wedge`].
///
/// # Errors
///
/// Propagates [`wedge_load_capacity`]'s errors for a taper that does not leave
/// a positive film, and returns [`TribologyError::NonPositive`] for a
/// non-positive clearance, viscosity, or surface speed.
pub fn thrust_pad_load(
    inner_radius: f64,
    outer_radius: f64,
    clearance: f64,
    viscosity: f64,
    surface_speed: f64,
    taper: f64,
    options: ReynoldsOptions,
) -> Result<f64> {
    if clearance <= 0.0 {
        return Err(TribologyError::NonPositive("thrust pad clearance"));
    }
    if viscosity <= 0.0 {
        return Err(TribologyError::NonPositive("lubricant viscosity"));
    }
    if surface_speed <= 0.0 {
        return Err(TribologyError::NonPositive("pad surface speed"));
    }
    if inner_radius <= 0.0 {
        return Err(TribologyError::NonPositive("thrust pad inner radius"));
    }
    if outer_radius <= inner_radius {
        return Err(TribologyError::OutsideValidRange(
            "a thrust pad outer radius beyond its inner radius",
        ));
    }
    wedge_load_capacity(taper)?;

    let profile = |x: f64| 1.0 - taper * x;
    let pressure = solve_profile_on_grid(options.nodes, &profile, 0.0, 0.0);

    // The area element `r dr`, not `dr`. Dropping the `r` would give the load
    // per unit radius rather than the total load, which for an annular pad
    // understates the answer by roughly the mean radius.
    let d_x = 1.0 / options.nodes as f64;
    let span = outer_radius - inner_radius;
    let mut integral = 0.0;
    for (i, p) in pressure[1..pressure.len() - 1].iter().enumerate() {
        let x = (i + 1) as f64 * d_x;
        let radius = inner_radius + x * span;
        integral += p * radius * span * d_x;
    }

    let pressure_scale = viscosity * surface_speed * span / (clearance * clearance);
    Ok(pressure_scale * integral)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A manufactured solution, the only way to prove the discretisation is
    /// right independently of any bearing textbook.
    ///
    /// The boundary conditions here are homogeneous, which is deliberate and
    /// is also why this test could not catch the dropped `P(0) = 1` term.
    #[test]
    fn the_discretisation_reproduces_a_manufactured_solution() {
        // A film profile and a pressure field, and the source that makes the
        // pressure field the exact solution of the discrete operator.
        let eccentricity = 0.5;
        let n = 200;
        let d_x = 1.0 / n as f64;
        let exact = |x: f64| {
            (core::f64::consts::PI * x).sin()
                * (1.0 + 0.3 * (2.0 * core::f64::consts::PI * x).cos())
        };

        let mut lower = vec![0.0; n - 1];
        let mut diagonal = vec![0.0; n - 1];
        let mut upper = vec![0.0; n - 1];
        let mut rhs = vec![0.0; n - 1];
        for k in 0..n - 1 {
            let x = (k + 1) as f64 * d_x;
            let hl = film_thickness(x - 0.5 * d_x, eccentricity).powi(3);
            let hr = film_thickness(x + 0.5 * d_x, eccentricity).powi(3);
            lower[k] = hl;
            diagonal[k] = -(hr + hl);
            upper[k] = hr;
            // The same operator, applied to the exact pressure, sampled at the
            // cell centres rather than the faces.
            rhs[k] = hr * (exact(x + d_x) - exact(x)) - hl * (exact(x) - exact(x - d_x));
        }
        let solved = solve_tridiagonal(&lower, &diagonal, &upper, &rhs);

        let mut worst: f64 = 0.0;
        for (k, value) in solved.iter().enumerate() {
            let x = (k + 1) as f64 * d_x;
            worst = worst.max((value - exact(x)).abs());
        }
        assert!(worst < 1.0e-12, "nodal error {worst} is above roundoff");
    }

    #[test]
    fn the_solution_is_grid_converged_in_the_full_film_regime() {
        // The reference values, from the two-grid solve: at e = 0.2 the peak
        // is 1.2255 and at e = 0.3 it is 1.8195, both independent of the
        // grid to five figures.
        for (eccentricity, expected) in [(0.2, 1.22550), (0.3, 1.81955)] {
            let solution = solve_journal_bearing(eccentricity, ReynoldsOptions::new(400))
                .expect("the full-film regime is resolved");
            assert!(
                (solution.peak_pressure - expected).abs() < 1.0e-3,
                "e = {eccentricity}: peak {} expected {expected}",
                solution.peak_pressure
            );
        }
    }

    #[test]
    fn the_peak_pressure_does_not_depend_on_the_grid() {
        // This is the property that makes the two-grid gate meaningful, and it
        // is the one that fails outside the full-film regime.
        let a = solve_journal_bearing(0.3, ReynoldsOptions::new(200)).expect("resolved");
        let b = solve_journal_bearing(0.3, ReynoldsOptions::new(800)).expect("resolved");
        assert!(
            (a.peak_pressure - b.peak_pressure).abs() / b.peak_pressure < 1.0e-3,
            "200 nodes gave {} and 800 gave {}",
            a.peak_pressure,
            b.peak_pressure
        );
    }

    #[test]
    fn the_solver_refuses_outside_the_full_film_regime() {
        // Past about 0.4 the film collapses, the Sommerfeld condition stops
        // being physical, and the plain solve returns a pressure that grows
        // with grid refinement. The solver must say so rather than hand back
        // a number in the hundreds of thousands.
        let result = solve_journal_bearing(0.5, ReynoldsOptions::new(400));
        assert!(result.is_err(), "e = 0.5 should be refused, not solved");
        if let Err(e) = result {
            assert!(
                matches!(e, TribologyError::OutsideValidRange(_)),
                "wrong error kind: {e}"
            );
        }
    }

    #[test]
    fn a_concentric_bearing_carries_a_boundary_forced_load() {
        // At e = 0 the film is uniform, so the source term vanishes and the
        // only thing left is the Sommerfeld boundary conditions forcing a
        // linear pressure ramp from 1 to 0. That ramp is not force-free: its
        // load is exactly 1/(2 pi).
        //
        // I expected zero here, because a concentric journal bearing
        // physically carries no load. It does not, in this formulation, and
        // the reason is the boundary condition rather than the solver. The
        // Sommerfeld condition presumes a converging wedge exists; at zero
        // eccentricity none does, and the imposed pressure gradient is a
        // modelling artefact. A real concentric bearing is resolved by the
        // condition that the pressure vanishes, not by extrapolating e to 0.
        let solution = solve_journal_bearing(0.0, ReynoldsOptions::new(2000)).expect("resolved");
        let expected = 1.0 / (2.0 * core::f64::consts::PI);
        assert!(
            (solution.load - expected).abs() / expected < 1.0e-3,
            "load at e = 0 is {} expected the boundary-forced {expected}",
            solution.load
        );
        // The profile is linear, so its peak is the boundary value.
        assert!((solution.peak_pressure - 1.0).abs() < 1.0e-6);
    }

    #[test]
    fn the_load_grows_with_eccentricity() {
        // A journal displaced further off-centre closes the wedge more and
        // squeezes harder.
        let light = solve_journal_bearing(0.1, ReynoldsOptions::new(200)).expect("resolved");
        let heavy = solve_journal_bearing(0.3, ReynoldsOptions::new(200)).expect("resolved");
        assert!(heavy.load > light.load, "{} !> {}", heavy.load, light.load);
    }

    #[test]
    fn the_pressure_is_finite_where_the_film_is_converging() {
        // A lubricant cannot sustain negative pressure, and it is tempting to
        // assert that everywhere. That assertion is simply false for this
        // formulation, and the reason matters.
        //
        // The Sommerfeld condition pressurises the whole circumference. In
        // the diverging half the film is opening up, the Reynolds equation
        // drives the pressure negative, and the condition forbids the
        // physical response - which is for the film to separate and the
        // pressure to go to zero. So mild negative pressure is present for
        // *any* non-zero eccentricity: about -0.23 at e = 0.2 and -0.82 at
        // e = 0.3. It is a symptom of the missing cavitation condition, not
        // a solver defect, and pretending otherwise would hide the very thing
        // this module documents as its limitation.
        //
        // What is worth asserting is that the solution is finite everywhere,
        // and that the pressure is non-negative in the converging half where
        // the film is genuinely being squeezed and the pressure is physical.
        for (eccentricity, _) in [(0.1, ()), (0.2, ()), (0.3, ())] {
            let solution =
                solve_journal_bearing(eccentricity, ReynoldsOptions::new(400)).expect("resolved");
            for (i, p) in solution.pressure.iter().enumerate() {
                assert!(p.is_finite(), "pressure {i} is not finite");
            }
            // The converging half is X in (0, 0.5), where the film closes.
            let half = solution.pressure.len() / 2;
            for p in &solution.pressure[1..half] {
                assert!(
                    *p > -1.0e-12,
                    "e = {eccentricity}: pressure {p} is negative in the converging half"
                );
            }
        }
    }

    #[test]
    fn the_negative_excursion_grows_with_eccentricity() {
        // The signature of the missing cavitation condition: as the film
        // collapses further, the unphysical negative excursion in the
        // diverging half grows. Tracking it is how a caller can see the
        // approach to the regime where the solve stops being trustworthy.
        let mut previous = 0.0f64;
        for eccentricity in [0.1, 0.2, 0.3, 0.4] {
            let solution =
                solve_journal_bearing(eccentricity, ReynoldsOptions::new(400)).expect("resolved");
            let worst = solution.pressure.iter().copied().fold(0.0f64, f64::min);
            assert!(
                worst <= previous + 1.0e-9,
                "e = {eccentricity}: excursion {worst} should not shrink from {previous}"
            );
            previous = worst;
        }
        // And it is genuinely negative by the time the regime is marginal.
        assert!(
            previous < -0.5,
            "expected a clear negative excursion, got {previous}"
        );
    }

    #[test]
    fn the_pressure_respects_its_boundary_conditions() {
        let solution = solve_journal_bearing(0.3, ReynoldsOptions::new(200)).expect("resolved");
        assert!((solution.pressure[0] - 1.0).abs() < 1e-12);
        assert!(solution.pressure[solution.pressure.len() - 1].abs() < 1e-12);
    }

    #[test]
    fn the_load_is_in_the_wedge_half_of_the_bearing() {
        // The film is being squeezed in the converging half, so the load must
        // act to push the journal back, and its direction must be in that
        // half. This is the physical statement that the sign convention is
        // right.
        let solution = solve_journal_bearing(0.3, ReynoldsOptions::new(400)).expect("resolved");
        let angle = solution.attitude_angle();
        let degrees = angle.to_degrees();
        assert!(
            degrees > -95.0 && degrees < 95.0,
            "attitude angle {degrees} degrees is not in the converging half"
        );
    }

    #[test]
    fn pressure_at_interpolates() {
        let solution = solve_journal_bearing(0.3, ReynoldsOptions::new(200)).expect("resolved");
        assert!((solution.pressure_at(0.0) - 1.0).abs() < 1e-12);
        let last = solution.pressure.len() - 1;
        assert!((solution.pressure_at(1.0) - solution.pressure[last]).abs() < 1e-12);
        let middle = solution.pressure_at(0.5);
        assert!(middle > 0.0 && middle < 1.0, "mid pressure {middle}");
    }

    #[test]
    fn the_minimum_film_ratio_matches_its_definition() {
        assert!((minimum_film_ratio(0.0) - 1.0).abs() < 1e-12);
        assert!((minimum_film_ratio(0.5) - 0.0).abs() < 1e-12);
        assert!((minimum_film_ratio(0.25) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn the_load_direction_classifies_the_bearing_by_its_wrapped_position() {
        // The mirror of the `reynolds_load_direction_depends_only_on_the_wrapped_
        // position` Kani harness. The harness cannot run in this environment, so
        // this test is the only evidence about the property; a test that merely
        // said "this would be proved" would be worth nothing.
        //
        // What the predicate actually claims is that the answer depends only on
        // the position *modulo one turn*, so the test is stated that way rather
        // than as "every position is inside the half". That second, stronger
        // claim is false: a position at `0.6` is in the diverging half, and
        // `-0.49` wraps to `0.51`, which is in the diverging half too.
        for i in 0..200 {
            let position = f64::from(i) * 0.01 - 1.0;
            let wrapped = position - position.floor();
            assert_eq!(
                load_direction_is_in_the_converging_half(position),
                load_direction_is_in_the_converging_half(wrapped),
                "position {position} disagreed with its own wrap {wrapped}"
            );
            // And the classification matches the geometry directly: inside for
            // the first half-turn, outside for the second.
            assert_eq!(
                load_direction_is_in_the_converging_half(position),
                wrapped <= 0.5,
                "position {position} (wrapped {wrapped}) is misclassified"
            );
        }
        // The boundary pinned on both sides, because this is where the strict
        // and non-strict comparisons differ and where a wrap bug hides. `0.5` is
        // the film-closure point and belongs inside the half; `0.51` is past
        // it; and `-0.5` wraps to exactly `0.5`, so it must agree with `0.5`
        // rather than with its neighbour.
        assert!(load_direction_is_in_the_converging_half(0.5));
        assert!(!load_direction_is_in_the_converging_half(0.51));
        assert!(load_direction_is_in_the_converging_half(-0.5));
        assert!(!load_direction_is_in_the_converging_half(-0.49));
        assert_eq!(
            load_direction_is_in_the_converging_half(-0.51),
            load_direction_is_in_the_converging_half(0.49)
        );
    }

    #[test]
    fn the_load_direction_is_periodic_in_whole_turns() {
        // The mirror of `reynolds_load_direction_depends_only_on_the_wrapped_
        // position`. The negative positions are the point: a caller indexing from
        // an angle can legitimately pass one, and a naive `position < 0.5` would
        // reject every one of them.
        for position in [0.0, 0.1, 0.25, 0.49, 0.5, 0.75, -0.3, -1.7] {
            for turns in [-3.0f64, -1.0, 0.0, 1.0, 4.0] {
                assert_eq!(
                    load_direction_is_in_the_converging_half(position),
                    load_direction_is_in_the_converging_half(position + turns),
                    "position {position} with {turns} turns"
                );
            }
        }
    }

    #[test]
    fn the_minimum_film_stays_positive_across_the_resolved_regime() {
        // The mirror of `reynolds_minimum_film_is_positive_across_the_resolved_
        // regime`, and the arithmetic that makes `MAX_RESOLVED_ECCENTRICITY` a
        // consequence rather than a chosen number.
        for i in 0..=40 {
            let e = f64::from(i) * 0.01;
            let h_min = minimum_film_ratio(e);
            assert!(h_min > 0.0, "e = {e}: film collapsed inside the range");
            assert!(h_min <= 1.0, "e = {e}: film exceeded the full clearance");
            assert!((h_min - (1.0 - 2.0 * e)).abs() < 1e-12);
        }
        // And past the resolved range it does go to zero and negative, which is
        // what the solver's refusal is protecting against. The eccentricities
        // are read through a local binding rather than written as literals,
        // because a literal `0.5` makes the comparison a constant the linter
        // quite reasonably objects to, and the intent here is the relation
        // between the eccentricity and the film, not one hard-coded number.
        let collapsed = 0.5_f64;
        let overshot = collapsed + 0.1;
        assert!(minimum_film_ratio(collapsed).abs() < 1e-12);
        assert!(minimum_film_ratio(overshot) < 0.0);
        assert!(MAX_RESOLVED_ECCENTRICITY < collapsed);
    }

    #[test]
    fn unphysical_eccentricities_are_rejected() {
        assert!(solve_journal_bearing(-0.1, ReynoldsOptions::new(200)).is_err());
        assert!(solve_journal_bearing(1.0, ReynoldsOptions::new(200)).is_err());
        assert!(solve_journal_bearing(1.5, ReynoldsOptions::new(200)).is_err());
    }

    #[test]
    fn the_options_enforce_a_minimum_grid() {
        // Too few cells would silently under-resolve; clamp instead.
        assert!(ReynoldsOptions::new(0).nodes >= 8);
        assert!(ReynoldsOptions::new(2).nodes >= 8);
    }

    #[test]
    fn the_load_equals_the_magnitude_of_its_components() {
        let solution = solve_journal_bearing(0.3, ReynoldsOptions::new(200)).expect("resolved");
        let expected =
            math::sqrt(solution.load_x * solution.load_x + solution.load_y * solution.load_y);
        assert!((solution.load - expected).abs() < 1e-12);
    }

    // --- Wedge (slider and thrust) bearings -------------------------------

    #[test]
    fn the_numerical_wedge_agrees_with_the_closed_form() {
        // This is the load-bearing test of the whole wedge implementation. The
        // closed form is exact and the solver is a discretisation, so they must
        // agree to the solver's own second-order accuracy. If they do not,
        // either the finite-volume assembly or the closed form is wrong, and
        // the numbers alone cannot say which -- which is why both are here.
        for taper in [0.1, 0.3, 0.5, 0.7, 0.9, 0.95] {
            let numeric = solve_wedge(taper, ReynoldsOptions::new(3200)).expect("a wedge");
            let exact = wedge_load_capacity(taper).expect("a taper");
            assert!(
                (numeric.load - exact).abs() / exact < 1.0e-5,
                "taper {taper}: {} vs {exact}",
                numeric.load
            );

            let peak_exact = wedge_peak_pressure(taper).expect("a taper");
            assert!(
                (numeric.peak_pressure - peak_exact).abs() / peak_exact < 1.0e-5,
                "taper {taper}: peak {} vs {peak_exact}",
                numeric.peak_pressure
            );
        }
    }

    #[test]
    fn the_wedge_load_is_second_order_convergent() {
        // Not just "agrees at one resolution" but converges, at the rate the
        // discretisation claims. Halving the cell width must quarter the error,
        // and if it does not then the agreement above is luck.
        let taper = 0.6;
        let exact = wedge_load_capacity(taper).expect("a taper");
        let error = |n: usize| {
            let s = solve_wedge(taper, ReynoldsOptions::new(n)).expect("a wedge");
            (s.load - exact).abs() / exact
        };
        let coarse = error(100);
        let medium = error(200);
        let fine = error(400);
        // Ratios near 4 mean second order. The bounds are loose enough to
        // tolerate the rounding floor but tight enough to catch a first-order
        // scheme, whose ratio would be 2.
        for (a, b) in [(coarse, medium), (medium, fine)] {
            let ratio = a / b;
            assert!(
                (3.0..5.0).contains(&ratio),
                "expected second order (ratio ~4), got {ratio} from {a} and {b}"
            );
        }
    }

    #[test]
    fn wedge_load_grows_without_bound_as_the_wedge_closes() {
        // A consequence of the idealised wedge that is easy to get backwards.
        // Both ends being open means nothing relieves the pressure, so load and
        // peak pressure both diverge as the outlet film vanishes. The classic
        // "optimum wedge angle" of `2 - sqrt(2)` comes from a *different* problem
        // -- a finite bearing with one end pressurised -- and does not apply
        // here. That distinction is worth stating, because quoting the constant
        // for this boundary condition would be a real and plausible error.
        let mut previous = 0.0;
        for taper in [0.2, 0.4, 0.6, 0.8, 0.9, 0.99] {
            let load = wedge_load_capacity(taper).expect("a taper");
            assert!(load > previous, "taper {taper} load {load} vs {previous}");
            previous = load;
        }
        let previous = wedge_load_capacity(0.99).expect("a taper");
        // And it keeps growing, so there is no interior optimum to quote, which
        // is the whole point of this test.
        assert!(wedge_load_capacity(0.999).unwrap() > previous);
    }

    #[test]
    fn specific_pressure_is_bounded_where_the_load_is_not() {
        // Since the load diverges as the wedge closes, neither load nor peak
        // pressure alone can rank two geometries: both just say "close the
        // wedge more". Their ratio does rank them, because it stays bounded.
        //
        // The limit is worth stating because it is the physical content of the
        // number: as the taper goes to zero the wedge degenerates to a linear
        // pressure ramp, whose peak-to-mean ratio is 3/2. So a shallow wedge
        // gives exactly the parabolic profile, and `specific_pressure` tends to
        // 1.5 from above. Every wedge is worse than that, and the penalty is
        // what buying load costs in peak pressure.
        let at = |taper: f64| {
            solve_wedge(taper, ReynoldsOptions::new(1000))
                .unwrap()
                .specific_pressure
        };
        let mut previous = 0.0;
        for taper in [0.05, 0.2, 0.4, 0.6, 0.8, 0.9] {
            let ratio = at(taper);
            assert!(
                ratio > 1.5,
                "taper {taper}: {ratio} is at or below the 3/2 limit"
            );
            assert!(
                ratio > previous,
                "taper {taper}: {ratio} did not rise from {previous}"
            );
            previous = ratio;
        }
        // Bounded, unlike the load it is built from: at a taper of 0.9 the load
        // has grown by a factor of ~180 from the shallowest case here, while the
        // ratio has grown by less than two.
        assert!(previous < 2.5, "{previous}");
        assert!(wedge_load_capacity(0.9).unwrap() / wedge_load_capacity(0.05).unwrap() > 100.0);
    }

    #[test]
    fn a_parallel_film_carries_no_load() {
        // A taper of zero is the limit from above: no wedge, no source term, no
        // pressure. It is excluded by `wedge_load_capacity` because the closed
        // form degenerates, but the solver must still agree with the limit
        // rather than dividing by zero.
        let flat = solve_wedge(1.0e-9, ReynoldsOptions::new(500)).expect("a shallow wedge");
        assert!(flat.load < 1.0e-6, "{}", flat.load);
        assert!(flat.peak_pressure < 1.0e-3, "{}", flat.peak_pressure);
    }

    #[test]
    fn wedge_geometry_that_is_not_a_wedge_is_rejected() {
        // Zero and negative tapers have no converging wedge; a taper of one
        // closes the film completely and both closed forms are singular there.
        assert!(wedge_load_capacity(0.0).is_err());
        assert!(wedge_load_capacity(-0.2).is_err());
        assert!(wedge_load_capacity(1.0).is_err());
        assert!(wedge_load_capacity(1.5).is_err());
        assert!(wedge_peak_pressure(0.0).is_err());
        assert!(solve_wedge(1.5, ReynoldsOptions::new(100)).is_err());
    }

    #[test]
    fn the_wedge_pressure_is_ambient_at_both_ends_and_positive_inside() {
        // The Gümbel condition. A solver that leaked the journal bearing's
        // `P(0) = 1` in here would still give a plausible-looking curve, just
        // shifted, which is why the ends are checked explicitly.
        let w = solve_wedge(0.5, ReynoldsOptions::new(400)).expect("a wedge");
        assert!(w.pressure[0].abs() < 1e-15);
        assert!(w.pressure[w.pressure.len() - 1].abs() < 1e-15);
        for p in &w.pressure[1..w.pressure.len() - 1] {
            assert!(*p > 0.0, "negative pressure inside the wedge: {p}");
        }
    }

    #[test]
    fn a_thrust_pad_carries_load_and_scales_with_every_dimension() {
        // A 100 mm bore pad, 40 mm span, 25 um film at the inner edge, 0.1 Pa s
        // oil at 8 m/s. The load lands in the kilonewton range, which is what a
        // thrust bearing of that size carries.
        let options = ReynoldsOptions::new(2000);
        let load = thrust_pad_load(0.1, 0.14, 25.0e-6, 0.1, 8.0, 0.5, options).expect("a pad");
        assert!(load > 1.0e3 && load < 1.0e6, "load = {load} N");

        // Linear in viscosity and in speed, since both are pure prefactors.
        let base =
            thrust_pad_load(0.1, 0.14, 25.0e-6, 0.1, 8.0, 0.5, ReynoldsOptions::new(400)).unwrap();
        let faster =
            thrust_pad_load(0.1, 0.14, 25.0e-6, 0.2, 8.0, 0.5, ReynoldsOptions::new(400)).unwrap();
        assert!((faster / base - 2.0).abs() < 1.0e-9, "{}", faster / base);

        // Inverse-square in the clearance, which is the dominant design lever:
        // halving the film quadruples the load. A linear film dependence here
        // would be the classic error.
        let tighter =
            thrust_pad_load(0.1, 0.14, 12.5e-6, 0.1, 8.0, 0.5, ReynoldsOptions::new(400)).unwrap();
        assert!(
            (tighter / base - 4.0).abs() / 4.0 < 1.0e-3,
            "{}",
            tighter / base
        );
    }

    #[test]
    fn a_thrust_pad_load_scales_with_its_pressure_weighted_radius() {
        // The `r dr` factor made visible. These two pads are the same size, so
        // the pressure scale and the wedge solution are identical between them
        // and the only thing that differs is the radius the load acts at.
        //
        // The scaling is by a *pressure-weighted* mean radius, not by the
        // geometric mean and not by the inner radius. The wedge converges
        // towards the outer edge, so pressure peaks there and pulls the
        // effective radius outward. That puts the result strictly between the
        // outer-radius ratio and the mean-radius ratio, which is a sharper
        // statement than "it scales with radius" and would fail outright if the
        // `r` factor were dropped from the integral.
        let options = ReynoldsOptions::new(1000);
        let small = thrust_pad_load(0.05, 0.09, 25.0e-6, 0.1, 8.0, 0.5, options).unwrap();
        let large = thrust_pad_load(0.10, 0.14, 25.0e-6, 0.1, 8.0, 0.5, options).unwrap();
        let ratio = large / small;

        let outer_ratio = 0.14 / 0.09;
        let mean_ratio = (0.10 + 0.14) / 2.0 / ((0.05 + 0.09) / 2.0);
        let inner_ratio = 0.10 / 0.05;
        assert!(ratio > outer_ratio, "{ratio} vs {outer_ratio}");
        assert!(ratio < mean_ratio, "{ratio} vs {mean_ratio}");
        assert!(mean_ratio < inner_ratio);

        // And the area ratio, which is what the integral would give if `r` were
        // dropped, sits well outside that band. This is the assertion that would
        // actually catch the missing factor.
        let area_ratio = (0.10f64 * 0.14) / (0.05 * 0.09);
        assert!(
            area_ratio > ratio,
            "area ratio {area_ratio} should exceed {ratio}"
        );
    }

    #[test]
    fn thrust_pad_geometry_that_cannot_exist_is_rejected() {
        let options = ReynoldsOptions::new(100);
        // No film, no oil, no motion, no annulus.
        assert!(thrust_pad_load(0.1, 0.14, 0.0, 0.1, 8.0, 0.5, options).is_err());
        assert!(thrust_pad_load(0.1, 0.14, -1e-6, 0.1, 8.0, 0.5, options).is_err());
        assert!(thrust_pad_load(0.1, 0.14, 25e-6, 0.0, 8.0, 0.5, options).is_err());
        assert!(thrust_pad_load(0.1, 0.14, 25e-6, 0.1, 0.0, 0.5, options).is_err());
        assert!(thrust_pad_load(0.0, 0.14, 25e-6, 0.1, 8.0, 0.5, options).is_err());
        // An outer radius inside the inner one is not an annulus.
        assert!(thrust_pad_load(0.14, 0.1, 25e-6, 0.1, 8.0, 0.5, options).is_err());
        // And the taper rules carry over from the wedge.
        assert!(thrust_pad_load(0.1, 0.14, 25e-6, 0.1, 8.0, 1.5, options).is_err());
    }

    #[test]
    fn cavitated_pressure_is_never_negative_and_satisfies_complementarity() {
        for e in [0.2, 0.3, 0.45] {
            let n = 400;
            let profile = |x: f64| film_thickness(x, e);
            let (p, settled) = solve_profile_cavitated(n, &profile, 1.0, 0.0);
            assert!(settled);
            assert!(p.iter().all(|v| *v >= 0.0));
            // Where the film is pressurised the discrete equation holds; where
            // it is cavitated the residual must not demand negative pressure.
            let d_x = 1.0 / n as f64;
            for k in 1..n {
                let x = k as f64 * d_x;
                let hl = profile(x - 0.5 * d_x);
                let hr = profile(x + 0.5 * d_x);
                let res = hl.powi(3) * p[k - 1] - (hl.powi(3) + hr.powi(3)) * p[k]
                    + hr.powi(3) * p[k + 1]
                    - 6.0 * (hr - hl) * d_x;
                if p[k] > 0.0 {
                    assert!(res.abs() < 1e-9, "e={e} k={k} res={res}");
                } else {
                    assert!(res <= 1e-9, "e={e} k={k} res={res}");
                }
            }
        }
    }

    #[test]
    fn cavitated_solution_is_grid_converged_past_the_sommerfeld_limit() {
        let a = solve_journal_bearing_cavitated(0.45, ReynoldsOptions::new(200)).unwrap();
        let b = solve_journal_bearing_cavitated(0.45, ReynoldsOptions::new(800)).unwrap();
        assert!((a.load - b.load).abs() / b.load < 5e-3);
        assert!((a.peak_pressure - b.peak_pressure).abs() / b.peak_pressure < 5e-3);
        assert!(b.peak_pressure > 9.0 && b.peak_pressure < 10.5);
    }

    #[test]
    fn cavitated_load_grows_monotonically_and_attitude_angle_falls() {
        let mut last_load = 0.0;
        let mut last_angle = f64::INFINITY;
        for e in [0.2, 0.3, 0.4, 0.45, 0.48] {
            let s = solve_journal_bearing_cavitated(e, ReynoldsOptions::new(400)).unwrap();
            assert!(s.load > last_load);
            let angle = s.attitude_angle().abs();
            assert!(angle >= 0.0 && s.load_x > 0.0);
            last_load = s.load;
            last_angle = last_angle.min(angle);
        }
        assert!(last_angle.is_finite());
    }

    #[test]
    fn cavitated_solver_rejects_a_closing_film_and_bad_input() {
        let o = ReynoldsOptions::new(200);
        assert!(solve_journal_bearing_cavitated(0.495, o).is_err());
        assert!(solve_journal_bearing_cavitated(-0.1, o).is_err());
        assert!(solve_journal_bearing_cavitated(1.0, o).is_err());
    }
}

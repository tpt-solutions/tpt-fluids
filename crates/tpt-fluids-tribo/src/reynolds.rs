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
        face_left[k] = film_thickness(x - 0.5 * d_x, eccentricity);
        face_right[k] = film_thickness(x + 0.5 * d_x, eccentricity);
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
    // The left Dirichlet value P(0) = 1 multiplies the first row's lower
    // coefficient and must be moved to the right-hand side. Omitting this is
    // the single most damaging bug available here: the solution stays smooth
    // and plausible but is wrong everywhere, and a manufactured-solution
    // test cannot catch it because that test uses homogeneous boundaries.
    rhs[0] -= lower[0];

    let interior = solve_tridiagonal(&lower, &diagonal, &upper, &rhs);

    let mut pressure = Vec::with_capacity(n + 1);
    pressure.push(1.0);
    pressure.extend_from_slice(&interior);
    pressure.push(0.0);
    pressure
}

/// The load components from a pressure distribution.
fn load_from_pressure(pressure: &[f64], nodes: usize) -> (f64, f64, f64) {
    let d_x = 1.0 / nodes as f64;
    let mut x = 0.0;
    let mut y = 0.0;
    for (i, p) in pressure.iter().enumerate().take(pressure.len() - 1).skip(1) {
        let position = i as f64 * d_x;
        let angle = 2.0 * core::f64::consts::PI * position;
        x += p * math::sin(angle) * d_x;
        y += p * math::cos(angle) * d_x;
    }
    (x, y, math::sqrt(x * x + y * y))
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
}

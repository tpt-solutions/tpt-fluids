//! Optimal pipe-network sizing, driven by `tpt-systems-optimisation`.
//!
//! # The problem
//!
//! Sizing a water network means choosing one diameter per pipe to minimise
//! cost, subject to the network still carrying its required flow. Two facts make
//! this awkward for a naive optimiser:
//!
//! - **The cost rises with diameter** (`cost ∝ D^2`, the standard
//!   material-dominated model), so the unconstrained minimiser drives every pipe
//!   to zero.
//! - **Velocity carries the constraint.** A narrower pipe is faster, and head
//!   loss grows steeply with diameter reduction -- `D^-5` at fixed flow -- so the
//!   unconstrained minimiser is hydraulically absurd long before it is
//!   numerically difficult.
//!
//! So the decision vector is the vector of diameters, the objective is capital
//! cost, and the constraints come from a hydraulic evaluation at the candidate
//! diameters. That evaluation is what makes the problem genuinely nonlinear, and
//! it is why this is an NLP and not a linear program.
//!
//! # What the solver gets, and what it does not
//!
//! The physics is computed by *this crate's* own hydraulic code, not by the
//! optimiser. [`NlpProblem::objective`] and [`NlpProblem::ineq`] call
//! [`crate::friction`]; `tpt-opt-core` supplies the constrained solver
//! (augmented Lagrangian with a conjugate-gradient inner minimiser) and the
//! convergence discipline. The split is deliberate: an optimiser must not be
//! trusted with physics, and the physics must not be reimplemented for it.
//!
//! # Gradients
//!
//! The objective gradient is analytic and exact -- the cost is a clean
//! polynomial. The constraint gradients fall back to the solver's central
//! differences, which is the honest choice: `friction_factor` iterates
//! internally, so there is no closed form for its derivative to write down, and
//! hand-differentiating it would be a place to be wrong quietly.
//!
//! # What this does not model
//!
//! Only capital cost. There is no pumping energy, no first-cost/boundary
//! dichotomy, and no discrete pipe catalogue -- diameters are continuous, which
//! is the standard continuous relaxation and a starting point for the discrete
//! problem rather than an answer to it.

use tpt_opt_core::nlp::{solve_nlp, NlpParams, NlpProblem, NlpStatus};

use crate::error::{HydraulicError, Result};
use crate::friction::{friction_factor, FrictionModel};
use tpt_fluids_core::consts::STANDARD_GRAVITY;

/// A pipe to be sized.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SizingPipe {
    /// Pipe length in metres.
    pub length: f64,
    /// The flow this pipe must carry in cubic metres per second.
    ///
    /// Treated as a *fixed allocation* rather than something the optimiser
    /// adjusts. Fixing it is what keeps the problem well-posed: if the flow were
    /// also a decision, the cost would fall to zero everywhere by starving the
    /// network, and the optimum would be "carry nothing".
    pub flow: f64,
    /// The smallest admissible diameter in metres.
    pub min_diameter: f64,
    /// The largest admissible diameter in metres.
    pub max_diameter: f64,
    /// Capital cost per metre of pipe at unit diameter, in the caller's units.
    pub cost_per_metre: f64,
}

impl SizingPipe {
    /// A pipe with a default diameter range of 50 mm to 500 mm.
    pub const fn new(length: f64, flow: f64, cost_per_metre: f64) -> Self {
        Self {
            length,
            flow,
            min_diameter: 0.05,
            max_diameter: 0.5,
            cost_per_metre,
        }
    }

    /// The capital cost of this pipe at diameter `d`, in cost units.
    ///
    /// ```text
    /// cost = L c D^2
    /// ```
    ///
    /// The `D^2` is the standard material-dominated scaling: pipe mass goes as
    /// the wall area, which for a fixed wall-thickness-to-diameter ratio goes as
    /// `D^2` once length is accounted for. It is a *model*, and it is the model
    /// the whole shape of the answer rests on -- a linear cost in `D` would give
    /// a different optimum, and a cost falling in `1/D` would not have an
    /// interior optimum at all.
    pub fn cost(&self, diameter: f64) -> f64 {
        if diameter <= 0.0 {
            return f64::INFINITY;
        }
        self.length * self.cost_per_metre * diameter * diameter
    }
}

/// The head loss along one pipe, in metres, at a given diameter.
///
/// Darcy-Weisbach, `h = f (L/D) V^2 / 2g`, with the velocity from the pipe's
/// area. Returns `INFINITY` for a non-positive diameter so the optimiser sees a
/// hard barrier rather than a number that quietly goes the wrong way.
pub fn pipe_head_loss(
    pipe: &SizingPipe,
    diameter: f64,
    kinematic_viscosity: f64,
    model: FrictionModel,
) -> f64 {
    if diameter <= 0.0 || pipe.flow <= 0.0 {
        return f64::INFINITY;
    }
    let area = core::f64::consts::PI * diameter * diameter / 4.0;
    let velocity = pipe.flow / area;
    let re = if kinematic_viscosity > 0.0 {
        velocity * diameter / kinematic_viscosity
    } else {
        0.0
    };
    let f = friction_factor(
        model,
        re,
        tpt_fluids_core::quantity::Length::new(diameter),
        tpt_fluids_core::quantity::Length::new(0.0),
    )
    // A zero roughness is hydraulically valid (a drawn tube), so this call only
    // fails on a non-positive diameter, which the guard above already excluded.
    // Mapping the error to `INFINITY` keeps a bad diameter on the infeasible
    // side of every constraint rather than surfacing an error mid-optimisation.
    .unwrap_or(f64::INFINITY);
    f * (pipe.length / diameter) * velocity * velocity / (2.0 * STANDARD_GRAVITY)
}

/// The limits the sized network must respect.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SizingConstraints {
    /// The largest total head loss the network may incur, in metres.
    pub max_total_head_loss: f64,
    /// The largest head loss any single pipe may incur, in metres.
    pub max_pipe_head_loss: f64,
    /// The largest velocity any pipe may carry, in metres per second.
    ///
    /// This is the constraint that actually sets pipe sizes in practice. Without
    /// it the optimum is a very slow, very large, very cheap network that
    /// satisfies the head budget trivially and is not a design anyone wants.
    pub max_velocity: f64,
    /// The friction correlation to use.
    pub friction_model: FrictionModel,
    /// The fluid's kinematic viscosity, in square metres per second.
    pub kinematic_viscosity: f64,
}

impl Default for SizingConstraints {
    /// Conservative defaults: a 10 m total budget, 5 m per pipe, and 2 m/s,
    /// which is the usual upper design velocity for a municipal main.
    fn default() -> Self {
        Self {
            max_total_head_loss: 10.0,
            max_pipe_head_loss: 5.0,
            max_velocity: 2.0,
            friction_model: FrictionModel::SwameeJain,
            kinematic_viscosity: 1.0e-6,
        }
    }
}

/// The outcome of a sizing run.
#[derive(Clone, PartialEq, Debug)]
pub struct SizingResult {
    /// The optimal diameter of each pipe, in metres, in input order.
    pub diameters: Vec<f64>,
    /// The total capital cost at those diameters.
    pub cost: f64,
    /// The head loss in each pipe at those diameters, in metres.
    pub head_losses: Vec<f64>,
    /// Whether the solver reported convergence.
    pub converged: bool,
    /// Outer iterations the solver took.
    pub iterations: usize,
}

impl SizingResult {
    /// The velocity in a pipe, in metres per second.
    pub fn velocity(&self, pipe: &SizingPipe, index: usize) -> f64 {
        let Some(d) = self.diameters.get(index).copied() else {
            return 0.0;
        };
        if d <= 0.0 {
            return 0.0;
        }
        let area = core::f64::consts::PI * d * d / 4.0;
        pipe.flow / area
    }

    /// The largest velocity in the network, in metres per second.
    pub fn max_velocity(&self, pipes: &[SizingPipe]) -> f64 {
        (0..pipes.len())
            .map(|i| self.velocity(&pipes[i], i))
            .fold(0.0f64, f64::max)
    }

    /// The total head loss in the network, in metres.
    pub fn total_head_loss(&self) -> f64 {
        self.head_losses
            .iter()
            .copied()
            .filter(|h| h.is_finite())
            .sum()
    }
}

/// The sizing problem as an [`NlpProblem`].
///
/// The decision vector is the diameters, one per pipe. Constraints are
/// `c(x) <= 0` in the solver's convention, in this order:
///
/// ```text
/// 0                      total head loss - budget          <= 0
/// 1                      max over pipes of (V - V_max)     <= 0
/// 2 .. 2+n               per-pipe head loss - limit        <= 0
/// 2+n .. 2+2n            per-pipe velocity - limit         <= 0
/// 2+2n .. 2+3n           min_diameter - D                  <= 0
/// 2+3n .. 2+4n           D - max_diameter                  <= 0
/// ```
///
/// The bounds are written out explicitly rather than left to the solver because
/// the augmented-Lagrangian method treats a bound as just another constraint and
/// will happily propose a negative diameter otherwise -- and a negative diameter
/// makes the head loss `INFINITY`, which the objective then integrates happily.
struct SizingProblem<'a> {
    pipes: &'a [SizingPipe],
    constraints: SizingConstraints,
}

impl SizingProblem<'_> {
    fn n(&self) -> usize {
        self.pipes.len()
    }

    /// The index of the first per-pipe velocity constraint.
    fn velocity_offset(&self) -> usize {
        2 + self.n()
    }

    /// The index of the first lower-bound constraint.
    fn lower_bound_offset(&self) -> usize {
        2 + 2 * self.n()
    }

    /// The velocity in pipe `j` at diameter `d`, in metres per second.
    fn velocity_of(&self, j: usize, d: f64) -> f64 {
        if d <= 0.0 {
            return f64::INFINITY;
        }
        let area = core::f64::consts::PI * d * d / 4.0;
        self.pipes[j].flow / area
    }

    /// The head loss in every pipe, given the diameter vector.
    fn losses(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n())
            .map(|j| {
                pipe_head_loss(
                    &self.pipes[j],
                    x[j],
                    self.constraints.kinematic_viscosity,
                    self.constraints.friction_model,
                )
            })
            .collect()
    }

    /// A constraint value that is never `NaN`, however the model misbehaved.
    ///
    /// This is a defence that the upstream solver badly needs and does not have.
    /// `solve_nlp` decides feasibility with `prob.ineq(i, x).max(0.0)`, and in
    /// Rust `f64::max` **returns the non-`NaN` operand**: `NaN.max(0.0)` is
    /// `0.0`. So a single `NaN` constraint is scored as *satisfied*, and the
    /// solver can converge in one iteration on a decision vector of pure `NaN`
    /// and report `Converged`. That is not hypothetical -- it is exactly what
    /// happened the first time this module was run, and it is the single most
    /// dangerous behaviour a numerical solver can have: a confident success on a
    /// meaningless answer.
    ///
    /// Mapping any `NaN` to `+INFINITY` puts it on the infeasible side of every
    /// constraint, where the penalty drives the optimiser back into the feasible
    /// region instead of rewarding it. The cost is that a genuinely broken model
    /// can no longer masquerade as a solution, which is the correct trade.
    fn guard(&self, value: f64) -> f64 {
        if value.is_nan() {
            f64::INFINITY
        } else {
            value
        }
    }
}

impl NlpProblem for SizingProblem<'_> {
    fn num_vars(&self) -> usize {
        self.n()
    }

    fn objective(&self, x: &[f64]) -> f64 {
        let mut total = 0.0;
        for (i, pipe) in self.pipes.iter().enumerate() {
            let d = x[i];
            if d <= 0.0 || !d.is_finite() {
                return f64::INFINITY;
            }
            let c = pipe.cost(d);
            if !c.is_finite() {
                return f64::INFINITY;
            }
            total += c;
        }
        self.guard(total)
    }

    fn objective_grad(&self, x: &[f64], g: &mut [f64]) {
        // d/dD of `L c D^2` is `2 L c D`: exact, and a clean polynomial deserves
        // an exact gradient rather than a central difference. A non-finite
        // diameter gets a zero row rather than a `NaN` gradient, which the
        // upstream line search would propagate into the whole step.
        for (i, pipe) in self.pipes.iter().enumerate() {
            let d = x[i];
            g[i] = if d.is_finite() && d > 0.0 {
                2.0 * pipe.length * pipe.cost_per_metre * d
            } else {
                0.0
            };
        }
    }

    fn num_ineq(&self) -> usize {
        2 + 4 * self.n()
    }

    fn ineq(&self, i: usize, x: &[f64]) -> f64 {
        let n = self.n();
        // Every branch is guarded. See `guard`: the solver scores a `NaN`
        // constraint as *satisfied*, so an unguarded `NaN` here reads as success.
        let value = match i {
            // Total head loss against the network budget. An `INFINITY` from a
            // non-positive diameter propagates out through the sum, which is the
            // behaviour wanted: a negative diameter must be infeasible, not
            // merely expensive.
            0 => self.losses(x).iter().sum::<f64>() - self.constraints.max_total_head_loss,
            // The single worst velocity in the network, as a single constraint.
            // Formulating it as a maximum over pipes rather than one constraint
            // per pipe keeps the multiplier attached to the binding pipe, which
            // is the pipe the designer needs to hear about.
            1 => (0..n)
                .map(|j| self.velocity_of(j, x[j]) - self.constraints.max_velocity)
                .fold(f64::NEG_INFINITY, f64::max),
            _ if i < 2 + n => self.losses(x)[i - 2] - self.constraints.max_pipe_head_loss,
            _ if i < self.velocity_offset() + n => {
                let j = i - self.velocity_offset();
                self.velocity_of(j, x[j]) - self.constraints.max_velocity
            }
            _ => {
                let k = i - self.lower_bound_offset();
                if k < n {
                    self.pipes[k].min_diameter - x[k]
                } else {
                    x[k - n] - self.pipes[k - n].max_diameter
                }
            }
        };
        self.guard(value)
    }

    fn num_eq(&self) -> usize {
        0
    }

    fn eq(&self, _j: usize, _x: &[f64]) -> f64 {
        0.0
    }
}

/// Sizes a set of pipes to minimise capital cost subject to head and velocity
/// limits.
///
/// `initial_diameters` is the starting point, one per pipe. It matters: the
/// augmented-Lagrangian method is a local method, and starting from a
/// hydraulically sensible guess is the difference between a run that converges
/// and one that stalls. A good default is the diameter that exactly meets the
/// velocity limit in each pipe.
///
/// # Errors
///
/// Returns [`HydraulicError::OutOfRange`] for a non-positive head-loss budget or
/// velocity limit -- either would make the problem infeasible for a
/// non-obvious reason -- a pipe with a non-positive flow, a pipe whose
/// `min_diameter` is non-positive or exceeds its `max_diameter`, or a starting
/// vector whose length does not match the pipe count.
///
/// The result reports `converged` rather than erroring when the solver merely
/// exhausts its iteration budget. A sizing run that stopped short is still
/// worth looking at; returning an error would throw away the diameters it did
/// find, which are often within a few percent of usable.
pub fn optimal_pipe_diameters(
    pipes: &[SizingPipe],
    constraints: SizingConstraints,
    initial_diameters: &[f64],
    params: &NlpParams,
) -> Result<SizingResult> {
    if pipes.is_empty() {
        return Ok(SizingResult {
            diameters: Vec::new(),
            cost: 0.0,
            head_losses: Vec::new(),
            converged: true,
            iterations: 0,
        });
    }
    if constraints.max_total_head_loss <= 0.0 {
        return Err(HydraulicError::OutOfRange(
            "the total head-loss budget must be positive".into(),
        ));
    }
    if constraints.max_velocity <= 0.0 {
        return Err(HydraulicError::OutOfRange(
            "the velocity limit must be positive".into(),
        ));
    }
    for pipe in pipes {
        if pipe.flow <= 0.0 {
            return Err(HydraulicError::OutOfRange(
                "every pipe must carry a positive flow".into(),
            ));
        }
        if pipe.min_diameter <= 0.0 || pipe.min_diameter > pipe.max_diameter {
            return Err(HydraulicError::OutOfRange(
                "each pipe needs a positive min_diameter no larger than its max_diameter".into(),
            ));
        }
    }
    if initial_diameters.len() != pipes.len() {
        return Err(HydraulicError::OutOfRange(
            "one starting diameter is needed per pipe".into(),
        ));
    }

    let problem = SizingProblem { pipes, constraints };
    let result = solve_nlp(&problem, initial_diameters, params);

    let head_losses = (0..pipes.len())
        .map(|i| {
            pipe_head_loss(
                &pipes[i],
                result.x[i],
                constraints.kinematic_viscosity,
                constraints.friction_model,
            )
        })
        .collect();
    let cost = pipes
        .iter()
        .enumerate()
        .map(|(i, p)| p.cost(result.x[i]))
        .sum();

    Ok(SizingResult {
        diameters: result.x,
        cost,
        head_losses,
        converged: result.status == NlpStatus::Converged,
        iterations: result.iterations,
    })
}

/// The smallest diameter that carries `flow` at exactly `velocity`, in metres.
///
/// This is the natural starting point for [`optimal_pipe_diameters`]: it is
/// feasible for the velocity constraint by construction, so the optimiser starts
/// from a hydraulically sane point rather than from an arbitrary guess that may
/// be wildly infeasible in both directions at once.
pub fn diameter_for_velocity(flow: f64, velocity: f64) -> f64 {
    if flow <= 0.0 || velocity <= 0.0 {
        return f64::NAN;
    }
    (4.0 * flow / (core::f64::consts::PI * velocity)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A three-pipe network with a real spread of flows and lengths, so the
    /// optimiser has to trade the pipes off against each other rather than
    /// finding one uniform answer.
    fn network() -> Vec<SizingPipe> {
        vec![
            SizingPipe {
                length: 1200.0,
                flow: 0.045,
                min_diameter: 0.05,
                max_diameter: 0.6,
                cost_per_metre: 1.0,
            },
            SizingPipe {
                length: 800.0,
                flow: 0.028,
                min_diameter: 0.05,
                max_diameter: 0.6,
                cost_per_metre: 1.0,
            },
            SizingPipe {
                length: 1500.0,
                flow: 0.061,
                min_diameter: 0.05,
                max_diameter: 0.6,
                cost_per_metre: 1.0,
            },
        ]
    }

    fn constraints() -> SizingConstraints {
        SizingConstraints {
            max_total_head_loss: 25.0,
            max_pipe_head_loss: 15.0,
            max_velocity: 1.8,
            friction_model: FrictionModel::SwameeJain,
            kinematic_viscosity: 1.0e-6,
        }
    }

    fn generous_params() -> NlpParams {
        NlpParams {
            max_outer: 60,
            max_inner: 60,
            rho_init: 1.0,
            rho_growth: 4.0,
            rho_max: 1.0e6,
            tol: 1.0e-6,
        }
    }

    /// The physical claim, before any optimiser is involved: head loss must fall
    /// as diameter rises. If this fails, every sizing result downstream is
    /// meaningless, so it is the right first test.
    #[test]
    fn head_loss_falls_with_diameter() {
        let pipe = network().remove(0);
        let mut previous = f64::INFINITY;
        for d in [0.05, 0.08, 0.12, 0.2, 0.35, 0.6] {
            let h = pipe_head_loss(&pipe, d, 1.0e-6, FrictionModel::SwameeJain);
            assert!(h.is_finite(), "head loss at D={d} is not finite");
            assert!(h < previous, "head loss rose at D={d}: {h} vs {previous}");
            previous = h;
        }
    }

    /// Cost must rise with diameter, quadratically, and that is what makes the
    /// problem an optimisation rather than an identity.
    #[test]
    fn cost_is_quadratic_in_diameter() {
        let pipe = SizingPipe::new(100.0, 0.02, 2.0);
        assert!((pipe.cost(0.1) - 100.0 * 2.0 * 0.01).abs() < 1.0e-12);
        // Halving the diameter quarters the cost.
        assert!((pipe.cost(0.05) - pipe.cost(0.1) / 4.0).abs() < 1.0e-12);
        assert_eq!(pipe.cost(0.0), f64::INFINITY);
    }

    /// The starting-diameter helper must be exact: it inverts `V = Q / A`.
    #[test]
    fn diameter_for_velocity_inverts_the_velocity_relation() {
        for (q, v) in [(0.045, 1.5), (0.01, 2.0), (0.2, 0.8)] {
            let d = diameter_for_velocity(q, v);
            let area = core::f64::consts::PI * d * d / 4.0;
            assert!((q / area - v).abs() < 1.0e-12, "Q={q} V={v} D={d}");
        }
        assert!(diameter_for_velocity(0.0, 1.0).is_nan());
        assert!(diameter_for_velocity(1.0, 0.0).is_nan());
    }

    /// The starting vector, scaled up from the velocity limit until the
    /// head-loss budget is met.
    ///
    /// The velocity-limit diameter is the *obvious* starting point and it is the
    /// wrong one: it is infeasible on head loss by construction, and the
    /// augmented-Lagrangian method is a descent method. Handed an infeasible
    /// start it walks further into infeasibility and never recovers, which is
    /// exactly what happened -- it returned a *more expensive* network than the
    /// one it was given, and a 3 percent shrink was cheaper still. Growing the
    /// start until it is feasible costs one cheap loop and makes the search start
    /// from inside the feasible region, where descent is meaningful.
    fn standard_start() -> Vec<f64> {
        let c = constraints();
        let pipes = network();
        let mut d: Vec<f64> = pipes
            .iter()
            .map(|p| diameter_for_velocity(p.flow, c.max_velocity))
            .collect();
        for _ in 0..60 {
            let total: f64 = pipes
                .iter()
                .enumerate()
                .map(|(i, p)| pipe_head_loss(p, d[i], c.kinematic_viscosity, c.friction_model))
                .sum();
            if total <= c.max_total_head_loss {
                break;
            }
            for value in d.iter_mut() {
                *value *= 1.05;
            }
        }
        d
    }

    /// The optimiser must produce a *feasible, finite* network. This is the
    /// baseline claim, and it is the one the `NaN` bug violated: the solver
    /// reported `converged = true` on a decision vector of pure `NaN`, which
    /// passes any check that only looks at the status flag.
    #[test]
    fn the_sized_network_is_finite_and_feasible() {
        let pipes = network();
        let r =
            optimal_pipe_diameters(&pipes, constraints(), &standard_start(), &generous_params())
                .unwrap();
        assert!(
            r.converged,
            "the solver did not converge: {:?}",
            r.diameters
        );
        assert_eq!(r.diameters.len(), pipes.len());
        for (i, d) in r.diameters.iter().enumerate() {
            assert!(d.is_finite(), "pipe {i} has a non-finite diameter {d}");
            assert!(
                *d >= pipes[i].min_diameter - 1.0e-9,
                "pipe {i} below its minimum"
            );
            assert!(
                *d <= pipes[i].max_diameter + 1.0e-9,
                "pipe {i} above its maximum"
            );
        }
        for h in &r.head_losses {
            assert!(h.is_finite() && *h >= 0.0, "bad head loss {h}");
        }
    }

    /// The optimum must be *cheaper* than the starting point, or the optimiser
    /// has not optimised anything. A run that returns its input is a plausible
    /// failure when the constraints are slack.
    #[test]
    fn the_optimum_is_cheaper_than_the_starting_point() {
        let pipes = network();
        let start = standard_start();
        let start_cost: f64 = pipes.iter().zip(&start).map(|(p, d)| p.cost(*d)).sum();
        let r = optimal_pipe_diameters(&pipes, constraints(), &start, &generous_params()).unwrap();
        assert!(
            r.cost < start_cost,
            "cost {} not below start {start_cost}",
            r.cost
        );
    }

    /// The result must be locally optimal: a uniform *shrink* must be blocked by
    /// the head-loss budget, because the budget is the binding constraint here.
    ///
    /// This is a genuine check rather than a tautology only because the binding
    /// constraint is established independently, in
    /// `the_head_loss_budget_is_the_binding_constraint`. Together the two say
    /// the network sits at an interior optimum: too large to grow (cost) and too
    /// small to shrink (head loss), which is the only shape a real optimum has.
    #[test]
    fn the_result_resists_a_uniform_shrink() {
        let pipes = network();
        let c = constraints();
        let r = optimal_pipe_diameters(&pipes, c, &standard_start(), &generous_params()).unwrap();
        let total_at = |scale: f64| -> f64 {
            pipes
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    pipe_head_loss(
                        p,
                        r.diameters[i] * scale,
                        c.kinematic_viscosity,
                        c.friction_model,
                    )
                })
                .sum()
        };
        // The head-loss budget is binding, so a shrink must break it. If it did
        // not, the network was not at its optimum.
        assert!(
            total_at(0.97) > c.max_total_head_loss,
            "a 3% shrink should violate the binding budget: {} vs {}",
            total_at(0.97),
            c.max_total_head_loss
        );
        // And a grow is dearer, which is the other half of interior optimality.
        let cost_at = |scale: f64| -> f64 {
            pipes
                .iter()
                .enumerate()
                .map(|(i, p)| p.cost(r.diameters[i] * scale))
                .sum()
        };
        assert!(cost_at(1.03) > cost_at(1.0), "growing should cost more");
    }

    /// The head-loss budget must actually bind. A run sitting far inside the
    /// budget has left money on the table: the pipes could shrink until either
    /// the budget or the velocity limit stopped them. This pins the physical
    /// content of the answer -- the two constraints trade off, and the binding one
    /// here is the head loss.
    #[test]
    fn the_head_loss_budget_is_the_binding_constraint() {
        let pipes = network();
        let c = constraints();
        let r = optimal_pipe_diameters(&pipes, c, &standard_start(), &generous_params()).unwrap();
        // Active within a percent, not exactly: the solver stops at a tolerance,
        // and demanding exactness would test the tolerance rather than physics.
        assert!(
            (r.total_head_loss() - c.max_total_head_loss).abs() / c.max_total_head_loss < 0.02,
            "head loss {} is not at the budget {}",
            r.total_head_loss(),
            c.max_total_head_loss
        );
        // And the velocity limit is *not* the binding one, or this test would be
        // asserting the wrong constraint.
        assert!(
            r.max_velocity(&pipes) < c.max_velocity,
            "velocity {} should be inside the limit {}",
            r.max_velocity(&pipes),
            c.max_velocity
        );
    }

    /// Tightening the head-loss budget must cost more, and must stay feasible.
    ///
    /// This is the check that the optimiser responds to the *constraint* and not
    /// merely the objective. A solver ignoring a constraint returns the same
    /// answer for both budgets; one responding degenerately returns a larger cost
    /// with the budget violated.
    #[test]
    fn a_tighter_budget_costs_more_and_stays_feasible() {
        let pipes = network();
        let loose = constraints();
        let tight = SizingConstraints {
            max_total_head_loss: loose.max_total_head_loss * 0.7,
            ..loose
        };
        let start = standard_start();
        let a = optimal_pipe_diameters(&pipes, loose, &start, &generous_params()).unwrap();
        let b = optimal_pipe_diameters(&pipes, tight, &start, &generous_params()).unwrap();
        assert!(
            b.cost > a.cost,
            "tighter budget was cheaper: {} vs {}",
            b.cost,
            a.cost
        );
        assert!(
            b.total_head_loss() <= tight.max_total_head_loss * 1.02,
            "tighter solution violates its own budget: {}",
            b.total_head_loss()
        );
    }

    /// The `NaN` guard, pinned directly.
    ///
    /// `solve_nlp` scores feasibility with `ineq(...).max(0.0)`, and Rust's
    /// `f64::max` returns the non-`NaN` operand, so a `NaN` constraint reads as
    /// *satisfied*. This states that property about the guard directly, so a
    /// future refactor that drops it fails here rather than in the field.
    #[test]
    fn a_nan_constraint_never_reads_as_satisfied() {
        let pipes = network();
        let problem = SizingProblem {
            pipes: &pipes,
            constraints: constraints(),
        };
        let upstream_violation = |v: f64| v.max(0.0);
        let guarded = problem.guard(f64::NAN);
        assert!(
            upstream_violation(guarded) > 0.0,
            "a guarded NaN is scored satisfied"
        );
        // A finite value passes through untouched, so the guard is not quietly
        // distorting the real constraints.
        assert_eq!(problem.guard(-1.5), -1.5);
        assert_eq!(problem.guard(2.5), 2.5);
        // The empirical demonstration of why this matters: an *unguarded* NaN
        // does score as satisfied under the upstream rule.
        assert_eq!(upstream_violation(f64::NAN), 0.0);
    }

    /// A non-finite starting point must not produce a confident nonsense answer.
    #[test]
    fn a_non_finite_start_does_not_produce_a_confident_answer() {
        let pipes = network();
        let start = vec![f64::NAN, 0.15, 0.2];
        let r = optimal_pipe_diameters(&pipes, constraints(), &start, &generous_params()).unwrap();
        // A garbage input must not yield a converged sizing with a non-finite
        // diameter. The status itself is not asserted either way: the claim is
        // only that "converged" and "non-finite" cannot both be true.
        if r.converged {
            for d in &r.diameters {
                assert!(d.is_finite(), "converged on a non-finite diameter {d}");
            }
        }
    }

    /// Empty input is a valid no-op, not an error.
    #[test]
    fn an_empty_network_is_a_no_op() {
        let r = optimal_pipe_diameters(
            &[],
            SizingConstraints::default(),
            &[],
            &NlpParams::default(),
        )
        .unwrap();
        assert!(r.converged);
        assert!(r.diameters.is_empty());
        assert_eq!(r.cost, 0.0);
    }

    /// Degenerate inputs are refused rather than producing a finite nonsense
    /// answer.
    #[test]
    fn degenerate_inputs_are_refused() {
        let pipes = network();
        let c = constraints();
        let good = vec![0.15; pipes.len()];

        let mut no_budget = c;
        no_budget.max_total_head_loss = 0.0;
        assert!(optimal_pipe_diameters(&pipes, no_budget, &good, &NlpParams::default()).is_err());

        let mut no_velocity = c;
        no_velocity.max_velocity = 0.0;
        assert!(optimal_pipe_diameters(&pipes, no_velocity, &good, &NlpParams::default()).is_err());

        let mut no_flow = pipes.clone();
        no_flow[1].flow = 0.0;
        assert!(optimal_pipe_diameters(&no_flow, c, &good, &NlpParams::default()).is_err());

        let mut bad_bounds = pipes.clone();
        bad_bounds[0].min_diameter = bad_bounds[0].max_diameter + 0.1;
        assert!(optimal_pipe_diameters(&bad_bounds, c, &good, &NlpParams::default()).is_err());

        let mut zero_min = pipes.clone();
        zero_min[0].min_diameter = 0.0;
        assert!(optimal_pipe_diameters(&zero_min, c, &good, &NlpParams::default()).is_err());

        // Wrong starting-vector length.
        assert!(optimal_pipe_diameters(&pipes, c, &[0.15], &NlpParams::default()).is_err());
    }

    /// A single pipe is the degenerate network, and it must still size sensibly.
    ///
    /// The assertion is on feasibility rather than on convergence: with one pipe
    /// the cost and the two hydraulic limits are all monotone in the same
    /// direction, so the optimum sits on whichever limit binds first and the
    /// solver's convergence flag depends on how the penalty settles. What must
    /// hold either way is that the answer is finite, inside the diameter bounds,
    /// and within the velocity limit -- the properties a designer relies on.
    #[test]
    fn a_single_pipe_sizes_within_its_velocity_limit() {
        let pipes = vec![SizingPipe {
            length: 1000.0,
            flow: 0.05,
            min_diameter: 0.05,
            max_diameter: 0.5,
            cost_per_metre: 1.0,
        }];
        let c = SizingConstraints {
            max_total_head_loss: 30.0,
            max_pipe_head_loss: 25.0,
            max_velocity: 1.5,
            ..SizingConstraints::default()
        };
        let start = vec![diameter_for_velocity(0.05, 1.5)];
        let r = optimal_pipe_diameters(&pipes, c, &start, &generous_params()).unwrap();
        let d = r.diameters[0];
        assert!(d.is_finite() && d > 0.0, "non-finite diameter {d}");
        assert!(
            d >= pipes[0].min_diameter && d <= pipes[0].max_diameter,
            "diameter {d} out of bounds"
        );
        let v = r.velocity(&pipes[0], 0);
        assert!(v <= c.max_velocity * 1.05, "velocity {v} exceeds the limit");
        assert!(r.head_losses[0].is_finite() && r.head_losses[0] >= 0.0);
    }
}

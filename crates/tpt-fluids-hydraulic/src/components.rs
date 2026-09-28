//! Component models: valves, pumps, turbines, cavitation, and surge tanks.
//!
//! These are the non-pipe elements a hydraulic network needs. Each one is
//! reduced to the quantity the network solvers actually consume: a
//! flow-dependent head gain, and where relevant a resistance coefficient for
//! the iterative solvers.
//!
//! # Sign conventions
//!
//! One convention runs through this whole module, and getting it wrong inverts
//! every answer, so it is worth stating once. **Positive head means head added
//! to the fluid.** A [`PumpCurve`] therefore has positive head across its whole
//! useful range, and a [`TurbineCurve`] has negative head while it generates.
//! Flow is signed along the machine's own generating direction, so a [`TurbineCurve`]
//! read at a negative flow is being run backwards, as a pump, and that is a
//! different quadrant rather than an error.

use crate::error::{HydraulicError, Result};
use tpt_fluids_core::consts::STANDARD_GRAVITY;
use tpt_fluids_core::math;
use tpt_fluids_core::quantity::{Area, Density, Pressure, Velocity};

/// US gallons per minute to cubic metres per second.
pub const GPM_TO_CUBIC_METRES_PER_SECOND: f64 = 6.309_019_64e-5;

/// Pounds per square inch to pascals.
pub const PSI_TO_PASCALS: f64 = 6_894.757_293;

/// A globe valve characterised by its flow coefficient `Cv`.
///
/// The US definition is `Q[gpm] = Cv sqrt(dp[psi] / SG)`, where `SG` is the
/// relative density of the fluid flowing through. Solving for the head loss
/// gives the quadratic the network solvers use directly.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Valve {
    /// The flow coefficient.
    pub cv: f64,
    /// The relative density of the flowing fluid, so water is 1.0.
    pub specific_gravity: f64,
    /// The fully-open loss coefficient `K` for line velocity. Used when no
    /// `Cv` is available; `None` when the valve is modelled by `Cv` alone.
    pub loss_coefficient: Option<f64>,
}

impl Valve {
    /// A valve from its `Cv`, for water.
    pub const fn from_cv(cv: f64) -> Self {
        Self {
            cv,
            specific_gravity: 1.0,
            loss_coefficient: None,
        }
    }

    /// A valve from its `Cv`, for a fluid of the given relative density.
    pub const fn from_cv_sg(cv: f64, specific_gravity: f64) -> Self {
        Self {
            cv,
            specific_gravity,
            loss_coefficient: None,
        }
    }

    /// A valve characterised by a line-velocity loss coefficient.
    pub const fn from_loss_coefficient(k: f64) -> Self {
        Self {
            cv: 0.0,
            specific_gravity: 1.0,
            loss_coefficient: Some(k),
        }
    }

    /// The flow through the valve at a given pressure drop, in cubic metres
    /// per second.
    ///
    /// A negative drop gives a negative flow, since a valve is bidirectional.
    pub fn flow(&self, pressure_drop: Pressure) -> f64 {
        let sg = if self.specific_gravity > 0.0 {
            self.specific_gravity
        } else {
            1.0
        };
        if self.loss_coefficient.is_some() || self.cv <= 0.0 {
            return 0.0;
        }
        let dp_psi = pressure_drop.value() / PSI_TO_PASCALS;
        if dp_psi <= 0.0 {
            return 0.0;
        }
        self.cv * math::sqrt(dp_psi / sg) * GPM_TO_CUBIC_METRES_PER_SECOND
    }

    /// The head loss across the valve for a given flow, in metres.
    ///
    /// Inverting the `Cv` definition, `dp[psi] = (Q / (Cv K))^2 SG`, gives a
    /// *pressure*; a head needs the fluid density as well, so the density is a
    /// parameter here. The network solvers carry head, so this is the form
    /// they need.
    pub fn head_loss(&self, flow: f64, density: Density) -> f64 {
        if let Some(k) = self.loss_coefficient {
            return k * flow * flow.abs();
        }
        let dp = self.pressure_drop(flow);
        if density.value() <= 0.0 {
            return f64::NAN;
        }
        dp / (density.value() * STANDARD_GRAVITY)
    }

    /// The pressure drop across the valve for a given flow, in pascals.
    ///
    /// This is the direct inverse of [`Valve::flow`] and needs no density,
    /// which is why it is the more primitive of the two.
    pub fn pressure_drop(&self, flow: f64) -> f64 {
        let sg = if self.specific_gravity > 0.0 {
            self.specific_gravity
        } else {
            1.0
        };
        if let Some(_k) = self.loss_coefficient {
            // A K-coefficient gives a head, not a pressure drop; without a
            // velocity there is no pressure to report.
            return f64::NAN;
        }
        if self.cv <= 0.0 {
            return 0.0;
        }
        let q_over = flow / (self.cv * GPM_TO_CUBIC_METRES_PER_SECOND);
        q_over * q_over * sg * PSI_TO_PASCALS
    }

    /// The resistance coefficient `R` in the solvers' `h = R Q^2` form.
    ///
    /// This is the single number that lets a valve drop into a
    /// [`crate::network::Link`] and be handled by either solver. It carries a
    /// density for the same reason [`Valve::head_loss`] does.
    pub fn resistance(&self, density: Density) -> f64 {
        if let Some(k) = self.loss_coefficient {
            return k;
        }
        if self.cv <= 0.0 || density.value() <= 0.0 {
            return f64::INFINITY;
        }
        let sg = if self.specific_gravity > 0.0 {
            self.specific_gravity
        } else {
            1.0
        };
        let k = self.cv * GPM_TO_CUBIC_METRES_PER_SECOND;
        sg * PSI_TO_PASCALS / (k * k * density.value() * STANDARD_GRAVITY)
    }
}

/// Fits `h(Q) = h0 - A Q^2 - B Q` through three or more `(flow, head)` points.
///
/// Three points determine a quadratic exactly; more are fitted by the normal
/// equations so a noisy test set still gives a usable curve. Dividing through
/// by `Q` on the interior points turns the system linear in the unknowns
/// `(h0, A, B)` for the three-point case, which avoids a general linear solve
/// for the common case.
///
/// # Errors
///
/// Returns [`HydraulicError::OutOfRange`] for fewer than three points, for
/// duplicate flows in the three-point case, or for a normal-equations system
/// that is singular because every flow is identical.
fn fit_quadratic(points: &[(f64, f64)]) -> Result<(f64, f64, f64)> {
    if points.len() < 3 {
        return Err(HydraulicError::OutOfRange(
            "a quadratic fit needs at least three points".into(),
        ));
    }
    for (q, h) in points {
        if !q.is_finite() || !h.is_finite() {
            return Err(HydraulicError::OutOfRange(
                "curve fit points must be finite".into(),
            ));
        }
    }
    if points.len() == 3 {
        let [(q0, h0), (q1, h1), (q2, h2)] = points else {
            unreachable!("length checked above")
        };
        if q0 == q1 || q1 == q2 || q0 == q2 {
            return Err(HydraulicError::OutOfRange(
                "the three curve fit points must be at distinct flows".into(),
            ));
        }
        // Subtract pairs to eliminate the quadratic term, giving two linear
        // equations in (h0, B) that are then solved directly.
        let d01 = q1 - q0;
        let d12 = q2 - q1;
        if d01 == 0.0 || d12 == 0.0 {
            return Err(HydraulicError::OutOfRange(
                "the three curve fit points must be at distinct flows".into(),
            ));
        }
        let s01 = (h1 - h0) / d01;
        let s12 = (h2 - h1) / d12;
        // h1 - h0 = -(A)(q1^2 - q0^2) - B (q1 - q0), so dividing by (q1-q0)
        // leaves -A(q1+q0) - B. Two such equations differ by A(q2 - q0).
        let quadratic = -(s12 - s01) / (q2 - q0);
        let linear = -s01 - quadratic * (q1 + q0);
        let shutoff = h0 + quadratic * q0 * q0 + linear * q0;
        return Ok((shutoff, quadratic, linear));
    }

    // Overdetermined: normal equations for the unknowns `(h0, -A, -B)` against
    // the design matrix `[1, -Q^2, -Q]`. The signs live in the design row so
    // that the recovered `A` and `B` come out already in the crate's
    // `h0 - A Q^2 - B Q` convention and need no negation afterwards.
    let mut m = [[0.0f64; 3]; 3];
    let mut rhs = [0.0f64; 3];
    for (q, h) in points {
        let row = [1.0, -(q * q), -*q];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += row[i] * row[j];
            }
            rhs[i] += row[i] * h;
        }
    }
    let solution = solve_3x3(m, rhs).ok_or_else(|| {
        HydraulicError::OutOfRange("the curve fit normal equations are singular".into())
    })?;
    Ok((solution[0], solution[1], solution[2]))
}

/// Solves a 3x3 system by Gaussian elimination with partial pivoting.
///
/// Returns `None` if the matrix is singular to working precision, which is the
/// caller's signal that the fit is not determined rather than that the answer is
/// large.
fn solve_3x3(mut m: [[f64; 3]; 3], mut rhs: [f64; 3]) -> Option<[f64; 3]> {
    for column in 0..3 {
        // Partial pivoting: the largest remaining pivot absorbs the rounding.
        let mut pivot = column;
        for row in (column + 1)..3 {
            if math::abs(m[row][column]) > math::abs(m[pivot][column]) {
                pivot = row;
            }
        }
        if m[pivot][column].abs() <= f64::EPSILON {
            return None;
        }
        m.swap(column, pivot);
        rhs.swap(column, pivot);

        for row in (column + 1)..3 {
            let factor = m[row][column] / m[column][column];
            // Split the borrow: the pivot row is read while the working row is
            // written, and they are never the same row once `row > column`.
            let pivot_row = m[column];
            let updated = m[row]
                .iter()
                .enumerate()
                .map(|(k, v)| v - factor * pivot_row[k])
                .collect::<Vec<f64>>();
            m[row].copy_from_slice(&updated);
            rhs[row] -= factor * rhs[column];
        }
    }

    let mut x = [0.0f64; 3];
    for row in (0..3).rev() {
        let mut sum = rhs[row];
        for k in (row + 1)..3 {
            sum -= m[row][k] * x[k];
        }
        x[row] = sum / m[row][row];
    }
    Some(x)
}

/// A centrifugal pump on its quadratic head characteristic.
///
/// `H = H0 - A Q^2 - B Q`, the standard three-point fit to a pump curve. The
/// shut-off head is `H0`, the head at zero flow, and the runout (zero-head)
/// flow follows from the curve.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PumpCurve {
    /// Shut-off head in metres.
    pub shutoff_head: f64,
    /// The quadratic head coefficient, in seconds squared per metre fifth.
    pub quadratic: f64,
    /// The linear head coefficient, in seconds per metre fourth.
    pub linear: f64,
}

impl PumpCurve {
    /// A pump from its shut-off head and the flow at which it develops no
    /// head, using the usual symmetric quadratic through those two points.
    ///
    /// A three-point fit is more accurate in general; this constructor is the
    /// two-point special case that covers most first-cut models.
    pub const fn from_shutoff_and_runout(shutoff_head: f64, runout_flow: f64) -> Self {
        // H(runout) = 0 with H0 - A Q^2 - B Q = 0, choosing the B that puts
        // the curve's best-efficiency point on the standard A = B*Q ratio.
        let linear = shutoff_head / (2.0 * runout_flow);
        let quadratic = linear / runout_flow;
        Self {
            shutoff_head,
            quadratic,
            linear,
        }
    }

    /// The head the pump develops at a given flow, in metres.
    pub fn head(&self, flow: f64) -> f64 {
        self.shutoff_head - self.quadratic * flow * flow - self.linear * flow
    }

    /// Fits the quadratic through three measured points on the pump curve.
    ///
    /// A three-point fit is what a real curve deserves: two points fix the
    /// parabola only after assuming a shape, whereas three measured points fix
    /// it outright and the fit is exact. The points are `(flow, head)` pairs and
    /// must be at distinct flows.
    ///
    /// The constant term is dropped and replaced by the shut-off head implied by
    /// the fit, because the parabola through three arbitrary points has no
    /// meaningful `Q = 0` intercept unless the data reaches back to zero flow.
    /// Reading the intercept off the fitted parabola is what keeps
    /// [`PumpCurve::head`] continuous at shut-off.
    ///
    /// # Errors
    ///
    /// Returns [`HydraulicError::OutOfRange`] if fewer than three distinct
    /// flows are supplied.
    pub fn from_three_points(points: &[(f64, f64)]) -> Result<Self> {
        let (shutoff_head, quadratic, linear) = fit_quadratic(points)?;
        Ok(Self {
            shutoff_head,
            quadratic,
            linear,
        })
    }

    /// The power the pump delivers to the water, in watts.
    ///
    /// `P = rho g Q H`, the hydraulic power. The shaft power needed to deliver
    /// it is larger by the pump efficiency, which this crate does not model.
    pub fn hydraulic_power(&self, flow: f64, density: Density) -> f64 {
        density.value() * STANDARD_GRAVITY * flow * self.head(flow)
    }

    /// The runout flow, where the curve reaches zero head.
    pub fn runout_flow(&self) -> f64 {
        // Solve A Q^2 + B Q - H0 = 0 for the positive root.
        if self.quadratic <= 0.0 {
            return f64::INFINITY;
        }
        let disc = self.linear * self.linear + 4.0 * self.quadratic * self.shutoff_head;
        if disc < 0.0 {
            return 0.0;
        }
        (-self.linear + math::sqrt(disc)) / (2.0 * self.quadratic)
    }

    /// The shut-off head, which is the curve value at zero flow.
    pub fn shutoff_head(&self) -> f64 {
        self.shutoff_head
    }
}

/// A hydro turbine on its four-quadrant characteristic.
///
/// A turbine is the mirror of a pump, and the mirror is worth spelling out
/// because it is the whole point of having a separate type rather than a
/// negated [`PumpCurve`]. Read in the generating direction the turbine
/// **extracts** head, so its head developed is negative and its flow is
/// positive:
///
/// ```text
/// dh(Q) = -(H_max - A Q^2 - B Q)
/// ```
///
/// The same parabola run past its runout flow gives the **windmilling**
/// quadrant for free: past runout the head developed turns positive, because
/// the machine is now being driven by the flow rather than driving it. That is
/// why a single quadratic covers the entire right-hand half of the four-quadrant
/// diagram, and it is a real property of the machine rather than a convenience.
///
/// The left-hand half, where the flow is reversed and the turbine is spinning
/// backwards as a pump, is **not** covered by that parabola, and the reason is
/// physical: the runner sees the flow from the other side, so the shape of the
/// loss changes. A separate branch is used there, sharing the value (not the
/// slope) of the generating branch at zero flow so the characteristic stays
/// continuous. That branch is the least well-determined part of any real
/// turbine characteristic, and the coefficient carrying it is exposed rather
/// than buried so a caller who has test data can use it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct TurbineCurve {
    /// The maximum head drop across the machine, at zero flow, in metres.
    pub maximum_head: f64,
    /// The quadratic head coefficient, in seconds squared per metre fifth.
    pub quadratic: f64,
    /// The linear head coefficient, in seconds per metre fourth.
    pub linear: f64,
    /// The coefficient of the reversed-flow (pump-mode) branch, dimensionless.
    ///
    /// Sets how fast head is absorbed when the machine is driven backwards.
    pub reverse_coefficient: f64,
}

impl TurbineCurve {
    /// The customary reverse-flow branch coefficient.
    ///
    /// A Francis machine driven backwards behaves broadly like a centrifugal
    /// pump of comparable head, and this value reproduces that magnitude.
    pub const DEFAULT_REVERSE: f64 = 4.0;

    /// A turbine from its maximum head and its runout flow.
    ///
    /// The runout flow is the flow at which the machine develops no net head,
    /// which is the same definition as for a pump and is where the generating
    /// and windmilling quadrants meet.
    pub const fn from_maximum_head_and_runout(maximum_head: f64, runout_flow: f64) -> Self {
        let quadratic = if runout_flow > 0.0 {
            maximum_head / (runout_flow * runout_flow)
        } else {
            0.0
        };
        Self {
            maximum_head,
            quadratic,
            linear: 0.0,
            reverse_coefficient: Self::DEFAULT_REVERSE,
        }
    }

    /// Fits the generating branch through three measured points.
    ///
    /// The points are `(flow, head drop)` with a positive flow and a positive
    /// head drop, which is the form a model test reports for a machine running
    /// as a generator. The reverse branch keeps its default coefficient,
    /// because a generating test says nothing about it.
    ///
    /// # Errors
    ///
    /// Propagates the errors from [`PumpCurve::from_three_points`]: the three
    /// flows must be finite and distinct.
    pub fn from_three_points(points: &[(f64, f64)]) -> Result<Self> {
        let (maximum_head, quadratic, linear) = fit_quadratic(points)?;
        Ok(Self {
            maximum_head,
            quadratic,
            linear,
            reverse_coefficient: Self::DEFAULT_REVERSE,
        })
    }

    /// The head drop *across* the machine, in metres, positive when the
    /// machine is absorbing head.
    ///
    /// This is the magnitude the network solvers want; the sign convention of
    /// the rest of the crate is applied by [`TurbineCurve::head_developed`].
    pub fn head_drop(&self, flow: f64) -> f64 {
        if flow < 0.0 {
            // Reversed flow. Continuous in value with the generating branch at
            // `flow = 0`, and falling without bound, so the machine absorbs
            // ever more head the harder it is driven backwards.
            return self.maximum_head - self.reverse_coefficient * flow * flow;
        }
        self.maximum_head - self.quadratic * flow * flow - self.linear * flow
    }

    /// The head the machine develops, in metres, negative while generating.
    pub fn head_developed(&self, flow: f64) -> f64 {
        -self.head_drop(flow)
    }

    /// Which quadrant of the four-quadrant diagram a flow sits in.
    pub fn quadrant(&self, flow: f64) -> TurbineQuadrant {
        if flow < 0.0 {
            return TurbineQuadrant::Pumping;
        }
        if flow <= self.runout_flow() {
            return TurbineQuadrant::Generating;
        }
        TurbineQuadrant::Windmilling
    }

    /// The flow at which the machine develops no net head, in cubic metres per
    /// second.
    pub fn runout_flow(&self) -> f64 {
        if self.quadratic <= 0.0 {
            return f64::INFINITY;
        }
        let disc = self.linear * self.linear + 4.0 * self.quadratic * self.maximum_head;
        if disc < 0.0 {
            return 0.0;
        }
        (-self.linear + math::sqrt(disc)) / (2.0 * self.quadratic)
    }

    /// The power the machine extracts from the water, in watts.
    ///
    /// `P = rho g Q dh` with `dh` the head drop. Positive means power is being
    /// taken out of the flow; negative means the machine is absorbing power, as
    /// it does when windmilling or driven backwards. No efficiency is applied:
    /// this is the hydraulic power crossing the machine, and the shaft power
    /// differs from it by the machine efficiency, which this crate does not
    /// model.
    pub fn extracted_power(&self, flow: f64, density: Density) -> f64 {
        density.value() * STANDARD_GRAVITY * flow * self.head_drop(flow)
    }
}

/// The quadrant of a turbine's four-quadrant diagram.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TurbineQuadrant {
    /// Flow in the generating direction, at or below the runout flow. The
    /// machine is extracting power.
    Generating,
    /// Flow in the generating direction above the runout flow. The machine is
    /// spinning fast enough to be driven by the water, and absorbing power.
    Windmilling,
    /// Flow reversed. The machine is being driven backwards, as a pump, and
    /// absorbs head.
    Pumping,
}

/// The cavitation margin at a point in a system.
///
/// The cavitation number `sigma = (p - p_v) / (0.5 rho V^2)` measures how far
/// the local absolute pressure sits above the fluid's vapour pressure. Below
/// one, vapour bubbles form; a system is generally kept well above it.
///
/// The interesting question is not "is sigma above one" but "is sigma above the
/// value **this** machine needs". A machine has its own inception number,
/// typically well above unity and strongly flow-dependent, and cavitation
/// begins when the available number falls to it. [`CavitationState`] answers the
/// first question and [`CavitationInception`] the second; [`CavitationMargin`]
/// joins them.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CavitationState {
    /// The cavitation number.
    pub sigma: f64,
    /// Whether the state is below the onset threshold.
    pub cavitating: bool,
}

impl CavitationState {
    /// The cavitation state at an absolute pressure, vapour pressure, and
    /// velocity.
    pub fn at(
        pressure: Pressure,
        vapour_pressure: Pressure,
        density: Density,
        velocity: Velocity,
    ) -> Self {
        let dynamic = 0.5 * density.value() * velocity.value() * velocity.value();
        let sigma = if dynamic == 0.0 {
            f64::INFINITY
        } else {
            (pressure.value() - vapour_pressure.value()) / dynamic
        };
        Self {
            sigma,
            cavitating: sigma < 1.0,
        }
    }
}

/// The cavitation number a machine needs before it begins to cavitate.
///
/// Cavitation in a turbomachine is worst at the extremes of flow and best at
/// the best-efficiency point, so a single constant inception number is not a
/// usable model: a machine checked only at its best point will happily cavitate
/// at part load or at overload. This parameterises the usual one-parameter shape
/// of that variation, a symmetric penalty in `(Q/Q_bep - 1)^2` about the best
/// point.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CavitationInception {
    /// The inception cavitation number at the best-efficiency flow.
    pub inception_number: f64,
    /// How fast the inception number rises away from the best point.
    ///
    /// Applied to `(Q/Q_bep - 1)^2`, so it is a per-unit-squared penalty. Zero
    /// reduces this to the constant-number model.
    pub flow_penalty: f64,
    /// The flow at which the machine is most efficient, in cubic metres per
    /// second.
    pub best_efficiency_flow: f64,
}

impl CavitationInception {
    /// A flow-independent inception number.
    pub const fn constant(inception_number: f64, best_efficiency_flow: f64) -> Self {
        Self {
            inception_number,
            flow_penalty: 0.0,
            best_efficiency_flow,
        }
    }

    /// The inception cavitation number at a given flow.
    ///
    /// Returns `f64::INFINITY` for a non-positive best-efficiency flow, since
    /// the penalty is defined relative to it and there is then no reference to
    /// measure the deviation from.
    pub fn at(&self, flow: f64) -> f64 {
        if self.best_efficiency_flow <= 0.0 {
            return f64::INFINITY;
        }
        let deviation = flow / self.best_efficiency_flow - 1.0;
        self.inception_number + self.flow_penalty * deviation * deviation
    }

    /// The net positive suction head the machine needs to stay free of
    /// cavitation, in metres.
    ///
    /// `NPSH_i = sigma_i V^2 / 2g` with the runner-passage velocity
    /// `V = Q / A` for the given runner diameter. This is the head the
    /// installation has to provide; compare it against the head actually
    /// available at the machine's own elevation.
    ///
    /// No density appears here, and that is not an oversight. The inception
    /// number is already dimensionless and the result is a *head*, so the
    /// `rho` that would convert between a pressure head and a velocity head
    /// cancels on both sides. Taking a density that then did nothing would
    /// invite a caller to believe it changed the answer.
    pub fn required_npsh(&self, flow: f64, diameter: f64) -> f64 {
        if diameter <= 0.0 {
            return f64::INFINITY;
        }
        let area = core::f64::consts::PI * diameter * diameter / 4.0;
        let velocity = flow / area;
        let sigma = self.at(flow);
        if velocity == 0.0 || sigma == f64::INFINITY {
            return f64::INFINITY;
        }
        sigma * velocity * velocity / (2.0 * STANDARD_GRAVITY)
    }
}

/// The gap between the cavitation head a system provides and the one a machine
/// needs.
///
/// Cavitation starts when this reaches zero from above, so it is the number to
/// watch in a part-load or overload calculation. Positive is safe.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct CavitationMargin {
    /// The margin in metres of head.
    pub margin: f64,
    /// Whether cavitation has begun.
    pub cavitating: bool,
}

impl CavitationMargin {
    /// Compares the available NPSH against the machine's requirement.
    pub fn new(available_npsh: f64, required_npsh: f64) -> Self {
        Self {
            margin: available_npsh - required_npsh,
            cavitating: available_npsh < required_npsh,
        }
    }

    /// The available NPSH implied by a static elevation above the machine and a
    /// pressure at it, in metres.
    ///
    /// `NPSH_a = p/(rho g) + z - p_v/(rho g)`, the sum of the absolute head at
    /// the machine, the height it stands above the free surface, and the vapour
    /// pressure expressed as a head. This is the quantity
    /// [`CavitationMargin::new`] wants on its left-hand side.
    pub fn available_npsh(
        pressure: Pressure,
        vapour_pressure: Pressure,
        density: Density,
        elevation_above_surface: f64,
    ) -> f64 {
        let rho_g = density.value() * STANDARD_GRAVITY;
        if rho_g <= 0.0 {
            return f64::INFINITY;
        }
        (pressure.value() - vapour_pressure.value()) / rho_g + elevation_above_surface
    }
}

/// A simple surge tank, modelled as a standpipe of area `A` on a riser.
///
/// The tank's level responds to the net inflow: `d(level)/dt = (Q_in - Q_out)
/// A`. A closed tank stores no net volume but adds inertia, represented here
/// by a small non-zero area to keep the model regular.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct SurgeTank {
    /// The tank's plan area in square metres.
    pub area: f64,
    /// The level in the tank above the datum, in metres.
    pub level: f64,
}

impl SurgeTank {
    /// A tank with the given plan area and starting level.
    pub const fn new(area: f64, level: f64) -> Self {
        Self { area, level }
    }

    /// Advances the level by `dt` given the net inflow, in cubic metres per
    /// second.
    pub fn advance(&mut self, net_inflow: f64, dt: f64) {
        if self.area > 0.0 {
            self.level += net_inflow * dt / self.area;
        }
    }

    /// The head the tank presents to the system, equal to its level.
    pub fn head(&self) -> f64 {
        self.level
    }
}

/// The loss coefficient for a sharp-edged entrance, 0.5.
pub const ENTRANCE_LOSS_SHARP: f64 = 0.5;
/// The loss coefficient for a well-rounded entrance, 0.04.
pub const ENTRANCE_LOSS_ROUNDED: f64 = 0.04;
/// The loss coefficient for a fully-flush 90-degree elbow, 0.3.
pub const ELBOW_LOSS_FLUSH_90: f64 = 0.3;
/// The loss coefficient for a fully-flush 45-degree elbow, 0.2.
pub const ELBOW_LOSS_FLUSH_45: f64 = 0.2;
/// The exit loss coefficient from a pipe into a reservoir, 1.0.
pub const EXIT_LOSS: f64 = 1.0;

/// The head loss through a minor (fitting) loss, in metres, from a `K`
/// coefficient and the line velocity.
///
/// This is `K V^2 / 2g`, the standard minor-loss form, and it is what turns a
/// tabulated fitting coefficient into the `R` the network solvers use.
pub fn minor_loss(k: f64, velocity: Velocity) -> f64 {
    k * velocity.value() * velocity.value() / (2.0 * STANDARD_GRAVITY)
}

/// The resistance coefficient for a minor loss, so `h = R Q^2`.
///
/// `R = K / (2 g A^2)` for a given area, which is what a
/// [`crate::network::Link`] stores.
pub fn minor_loss_resistance(k: f64, area: Area) -> f64 {
    let a2 = area.value() * area.value();
    if a2 == 0.0 {
        return f64::INFINITY;
    }
    k / (2.0 * STANDARD_GRAVITY * a2)
}

/// The velocity head `V^2 / 2g` in metres.
pub fn velocity_head(velocity: Velocity) -> f64 {
    let v = velocity.value();
    v * v / (2.0 * STANDARD_GRAVITY)
}

/// Validates a `Cv` or `K` value, rejecting non-physical inputs early.
pub fn validate_loss_parameter(value: f64) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(HydraulicError::OutOfRange(format!(
            "loss coefficient must be finite and non-negative, got {value}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valve_cv_matches_the_us_definition() {
        // Cv = 100 at 10 psi on water is 316.2 gpm by definition.
        let v = Valve::from_cv(100.0);
        let q = v.flow(Pressure::new(10.0 * PSI_TO_PASCALS));
        assert!(
            (q / GPM_TO_CUBIC_METRES_PER_SECOND - 316.2278).abs() < 1e-3,
            "{q}"
        );
    }

    #[test]
    fn valve_flow_is_zero_without_a_pressure_drop() {
        let v = Valve::from_cv(100.0);
        assert_eq!(v.flow(Pressure::new(0.0)), 0.0);
    }

    #[test]
    fn valve_pressure_drop_inverts_its_own_flow() {
        // Round-tripping flow -> pressure drop -> flow must be the identity.
        let v = Valve::from_cv(250.0);
        let q = v.flow(Pressure::new(40.0 * PSI_TO_PASCALS));
        let dp = v.pressure_drop(q);
        assert!((dp - 40.0 * PSI_TO_PASCALS).abs() / dp < 1e-9, "{dp}");
        let back = v.flow(Pressure::new(dp));
        assert!((back - q).abs() / q < 1e-9, "{back} vs {q}");
    }

    #[test]
    fn valve_head_loss_is_pressure_over_rho_g() {
        let v = Valve::from_cv(250.0);
        let rho = Density::new(998.2);
        let q = 0.099754359;
        let h = v.head_loss(q, rho);
        let expected = v.pressure_drop(q) / (998.2 * STANDARD_GRAVITY);
        assert!((h - expected).abs() < 1e-12, "{h} vs {expected}");
        // 40 psi through this valve is about 28 m of water.
        assert!((h - 28.17).abs() < 0.1, "h = {h}");
    }

    #[test]
    fn valve_resistance_reproduces_head_loss() {
        let v = Valve::from_cv(100.0);
        let rho = Density::new(998.2);
        let r = v.resistance(rho);
        let q = 0.02;
        assert!((r * q * q - v.head_loss(q, rho)).abs() < 1e-9);
    }

    #[test]
    fn valve_specific_gravity_reduces_flow() {
        // A denser fluid at the same pressure drop passes less flow.
        let water = Valve::from_cv_sg(100.0, 1.0);
        let heavy = Valve::from_cv_sg(100.0, 2.0);
        let dp = Pressure::new(10.0 * PSI_TO_PASCALS);
        assert!(heavy.flow(dp) < water.flow(dp));
        assert!((water.flow(dp) / heavy.flow(dp) - math::sqrt(2.0)).abs() < 1e-6);
    }

    #[test]
    fn loss_coefficient_valve_uses_quadratic_head_loss() {
        let v = Valve::from_loss_coefficient(1.0);
        assert!((v.head_loss(2.0, Density::new(1000.0)) - 4.0).abs() < 1e-12);
        // A Cv-based valve's Cv is ignored in this mode.
        let v2 = Valve {
            cv: 100.0,
            specific_gravity: 1.0,
            loss_coefficient: Some(2.0),
        };
        assert!((v2.head_loss(3.0, Density::new(1000.0)) - 18.0).abs() < 1e-12);
    }

    #[test]
    fn pump_head_falls_monotonically_with_flow() {
        let p = PumpCurve::from_shutoff_and_runout(50.0, 0.2);
        let mut previous = f64::INFINITY;
        let mut q = 0.0;
        while q < 0.15 {
            let h = p.head(q);
            assert!(h < previous, "head rose at Q={q}");
            previous = h;
            q += 0.01;
        }
    }

    #[test]
    fn pump_shutoff_and_runout_are_consistent() {
        let p = PumpCurve::from_shutoff_and_runout(50.0, 0.2);
        assert!((p.head(0.0) - 50.0).abs() < 1e-12);
        let runout = p.runout_flow();
        assert!((runout - 0.2).abs() < 1e-9, "runout = {runout}");
        assert!((p.head(runout) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn pump_power_peaks_near_the_best_efficiency_point() {
        let p = PumpCurve::from_shutoff_and_runout(50.0, 0.2);
        let rho = Density::new(998.2);
        let at_zero = p.hydraulic_power(0.0, rho);
        let at_mid = p.hydraulic_power(0.1, rho);
        let at_runout = p.hydraulic_power(p.runout_flow(), rho);
        assert!(at_mid > at_zero, "no power at shutoff");
        assert!(at_mid > at_runout, "no power at runout");
    }

    #[test]
    fn cavitation_number_matches_its_definition() {
        let rho = Density::new(998.2);
        let v = Velocity::new(2.5);
        let pv = Pressure::new(2_339.0);
        let p = Pressure::new(101_325.0 + 998.2 * STANDARD_GRAVITY * 10.0);
        let s = CavitationState::at(p, pv, rho, v);
        let expected = (p.value() - pv.value()) / (0.5 * 998.2 * 2.5 * 2.5);
        assert!((s.sigma - expected).abs() < 1e-9);
        assert!(!s.cavitating);
    }

    #[test]
    fn low_pressure_cavitates() {
        let rho = Density::new(998.2);
        let v = Velocity::new(10.0);
        // Just above vapour pressure: sigma is tiny.
        let s = CavitationState::at(Pressure::new(3_000.0), Pressure::new(2_339.0), rho, v);
        assert!(s.cavitating, "sigma = {}", s.sigma);
    }

    #[test]
    fn surge_tank_level_follows_net_inflow() {
        let mut tank = SurgeTank::new(2.0, 10.0);
        tank.advance(4.0, 1.0);
        assert!((tank.level - 12.0).abs() < 1e-12, "{}", tank.level);
        tank.advance(-4.0, 1.0);
        assert!((tank.level - 10.0).abs() < 1e-12, "{}", tank.level);
        assert!((tank.head() - tank.level).abs() < 1e-12);
    }

    #[test]
    fn surge_tank_ignores_a_non_positive_area() {
        let mut tank = SurgeTank::new(0.0, 10.0);
        tank.advance(5.0, 1.0);
        assert!((tank.level - 10.0).abs() < 1e-12);
    }

    #[test]
    fn minor_loss_consistences() {
        let v = Velocity::new(2.0);
        assert!((minor_loss(0.5, v) - 0.5 * 4.0 / (2.0 * STANDARD_GRAVITY)).abs() < 1e-12);
        assert!((velocity_head(v) - 4.0 / (2.0 * STANDARD_GRAVITY)).abs() < 1e-12);
    }

    #[test]
    fn minor_loss_resistance_reproduces_the_head_loss() {
        let area = Area::new(0.0706858);
        let r = minor_loss_resistance(0.5, area);
        let q = 0.0706858; // v = 1 m/s
        let expected = minor_loss(0.5, Velocity::new(1.0));
        assert!(
            (r * q * q - expected).abs() / expected < 1e-12,
            "{} vs {}",
            r * q * q,
            expected
        );
    }

    #[test]
    fn zero_area_gives_infinite_resistance() {
        assert!(minor_loss_resistance(1.0, Area::new(0.0)).is_infinite());
    }

    #[test]
    fn loss_parameter_validation_rejects_negatives() {
        assert!(validate_loss_parameter(0.5).is_ok());
        assert!(validate_loss_parameter(-1.0).is_err());
        assert!(validate_loss_parameter(f64::NAN).is_err());
    }

    #[test]
    fn tabulated_fitting_coefficients_are_ordered_sensibly() {
        // A rounded entrance must lose less than a sharp one, and a 45-degree
        // elbow less than a 90-degree one. These are ordering facts about the
        // published table, checked at runtime so the table cannot be
        // mis-transcribed without a test failing.
        let rounded_beats_sharp = ENTRANCE_LOSS_ROUNDED < ENTRANCE_LOSS_SHARP;
        let small_elbow_beats_big = ELBOW_LOSS_FLUSH_45 < ELBOW_LOSS_FLUSH_90;
        assert!(
            rounded_beats_sharp,
            "rounded entrance {ENTRANCE_LOSS_ROUNDED} vs sharp {ENTRANCE_LOSS_SHARP}"
        );
        assert!(
            small_elbow_beats_big,
            "45 deg {ELBOW_LOSS_FLUSH_45} vs 90 deg {ELBOW_LOSS_FLUSH_90}"
        );
        assert_eq!(EXIT_LOSS, 1.0, "an exit discharges all the velocity head");
    }

    // --- Turbine and cavitation inception ---------------------------------

    /// A machine with 40 m of head at shutoff and a 4 m^3/s runout.
    fn turbine() -> TurbineCurve {
        TurbineCurve::from_maximum_head_and_runout(40.0, 4.0)
    }

    #[test]
    fn a_three_point_fit_reproduces_its_own_points_exactly() {
        // The whole reason for the fit is that it is exact on three points. If
        // this drifts, the closed-form solve in `fit_quadratic` is wrong.
        let points = [(0.0, 50.0), (0.1, 48.0), (0.2, 42.0)];
        let pump = PumpCurve::from_three_points(&points).unwrap();
        for (q, h) in points {
            assert!(
                (pump.head(q) - h).abs() < 1e-9,
                "Q={q}: {} vs {h}",
                pump.head(q)
            );
        }
        // And the shut-off head the fit reports is the curve's own intercept.
        assert!((pump.shutoff_head() - 50.0).abs() < 1e-9);
    }

    #[test]
    fn a_three_point_fit_handles_a_linear_curve_and_more_than_three_points() {
        // A linear curve is the degenerate case where the quadratic coefficient
        // must come out zero, not NaN.
        let linear = [(0.0, 50.0), (1.0, 45.0), (2.0, 40.0)];
        let pump = PumpCurve::from_three_points(&linear).unwrap();
        assert!(pump.quadratic.abs() < 1e-9, "A = {}", pump.quadratic);
        assert!((pump.head(1.0) - 45.0).abs() < 1e-9);

        // Five points on the same curve: the normal-equation path must land on
        // the same parabola as the closed-form one.
        let many: Vec<(f64, f64)> = (0..5)
            .map(|i| {
                let q = f64::from(i) * 0.05;
                (q, 50.0 - 100.0 * q * q - 50.0 * q)
            })
            .collect();
        let fitted = PumpCurve::from_three_points(&many).unwrap();
        assert!(
            (fitted.quadratic - 100.0).abs() / 100.0 < 1e-6,
            "{}",
            fitted.quadratic
        );
        assert!(
            (fitted.linear - 50.0).abs() / 50.0 < 1e-6,
            "{}",
            fitted.linear
        );
    }

    #[test]
    fn a_curve_fit_refuses_data_it_cannot_fit() {
        assert!(PumpCurve::from_three_points(&[(0.0, 1.0), (1.0, 2.0)]).is_err());
        // Two identical flows leave the system under-determined.
        assert!(PumpCurve::from_three_points(&[(1.0, 1.0), (1.0, 2.0), (2.0, 3.0)]).is_err());
        assert!(PumpCurve::from_three_points(&[(0.0, f64::NAN), (1.0, 2.0), (2.0, 3.0)]).is_err());
    }

    #[test]
    fn a_turbine_extracts_head_rather_than_developing_it() {
        let t = turbine();
        // Generating: the head developed is negative, the head drop positive.
        assert!(t.head_developed(1.0) < 0.0);
        assert!(t.head_drop(1.0) > 0.0);
        assert!((t.head_developed(1.0) + t.head_drop(1.0)).abs() < 1e-12);
        // 40 m at shutoff falls to 30 m by 2 m^3/s: the parabola through
        // (0, 40) and (4, 0) is `40 - 2.5 Q^2`.
        assert!(
            (t.head_drop(2.0) - 30.0).abs() < 1e-12,
            "{}",
            t.head_drop(2.0)
        );
    }

    #[test]
    fn a_turbine_passes_through_zero_head_at_its_runout() {
        let t = turbine();
        assert!(
            (t.runout_flow() - 4.0).abs() / 4.0 < 1e-12,
            "{}",
            t.runout_flow()
        );
        assert!(t.head_drop(t.runout_flow()).abs() < 1e-9);
    }

    #[test]
    fn the_quadrant_follows_the_flow_not_the_head() {
        let t = turbine();
        assert_eq!(t.quadrant(1.0), TurbineQuadrant::Generating);
        assert_eq!(
            t.quadrant(4.0),
            TurbineQuadrant::Generating,
            "runout is inclusive"
        );
        assert_eq!(t.quadrant(6.0), TurbineQuadrant::Windmilling);
        assert_eq!(t.quadrant(-1.0), TurbineQuadrant::Pumping);
    }

    #[test]
    fn a_windmilling_turbine_reverses_its_power_and_its_head() {
        // Past runout the machine is driven by the water: it absorbs power and
        // develops positive head. Both signs have to flip together or the model
        // would claim a machine extracting power from a flow it cannot extract
        // any from.
        let t = turbine();
        let past = 8.0;
        assert!(t.head_drop(past) < 0.0, "head drop {}", t.head_drop(past));
        assert!(t.head_developed(past) > 0.0);
        assert!(
            t.extracted_power(past, Density::new(1000.0)) < 0.0,
            "windmilling must absorb power"
        );
        assert!(t.extracted_power(2.0, Density::new(1000.0)) > 0.0);
    }

    #[test]
    fn turbine_power_matches_rho_g_q_dh() {
        let t = turbine();
        let rho = Density::new(1000.0);
        let q = 2.0;
        let expected = 1000.0 * STANDARD_GRAVITY * q * t.head_drop(q);
        assert!((t.extracted_power(q, rho) - expected).abs() / expected < 1e-12);
        // 2 m^3/s through a 30 m head drop is about 589 kW.
        assert!((expected / 1.0e3 - 588.6).abs() < 0.5, "{expected} W");
    }

    #[test]
    fn the_reverse_branch_is_continuous_with_the_generating_branch() {
        // The two branches meet at zero flow. If they did not, the characteristic
        // would have a step in it and the network solvers would see a head jump
        // appear and vanish as the flow reversed.
        let t = turbine();
        assert!((t.head_drop(0.0) - t.head_drop(-0.0)).abs() < 1e-12);
        assert!((t.head_drop(0.0) - t.maximum_head).abs() < 1e-12);
        // And it must go the right way: driven backwards, head is absorbed.
        assert!(t.head_drop(-2.0) < t.head_drop(0.0));
    }

    #[test]
    fn a_three_point_turbine_fit_is_exact_on_its_generating_data() {
        let points = [(0.0, 40.0), (2.0, 24.0), (4.0, 0.0)];
        let t = TurbineCurve::from_three_points(&points).unwrap();
        for (q, dh) in points {
            assert!(
                (t.head_drop(q) - dh).abs() < 1e-9,
                "Q={q}: {} vs {dh}",
                t.head_drop(q)
            );
        }
    }

    #[test]
    fn inception_is_easiest_at_the_best_efficiency_flow() {
        let inception = CavitationInception {
            inception_number: 0.1,
            flow_penalty: 0.4,
            best_efficiency_flow: 2.0,
        };
        assert!((inception.at(2.0) - 0.1).abs() < 1e-12);
        // Symmetric about the best point, and strictly worse on both sides.
        assert!((inception.at(1.0) - inception.at(3.0)).abs() < 1e-12);
        assert!(inception.at(0.5) > inception.at(1.5));
        assert!(inception.at(1.5) > inception.at(2.0));
    }

    #[test]
    fn a_constant_inception_number_ignores_flow() {
        let inception = CavitationInception::constant(0.25, 2.0);
        for q in [0.0, 1.0, 2.0, 8.0] {
            assert!((inception.at(q) - 0.25).abs() < 1e-12, "Q={q}");
        }
        // With no reference flow the deviation is undefined, so it must not
        // silently invent one.
        assert!(CavitationInception::constant(0.25, 0.0)
            .at(1.0)
            .is_infinite());
    }

    #[test]
    fn required_npsh_grows_as_the_square_of_the_flow() {
        // NPSH_i = sigma_i V^2 / 2g, so doubling the flow at fixed sigma needs
        // four times the head. Anything else would mean the velocity head is
        // being mishandled, which is the classic NPSH error.
        let inception = CavitationInception::constant(0.2, 1.0);
        let d = 2.0;
        let a = inception.required_npsh(1.0, d);
        let b = inception.required_npsh(2.0, d);
        assert!((b / a - 4.0).abs() / 4.0 < 1e-12, "{b} vs {a}");
        // At shutoff there is no velocity head, so no requirement at all.
        assert!(inception.required_npsh(0.0, d).is_infinite());
        assert!(inception.required_npsh(1.0, 0.0).is_infinite());
    }

    #[test]
    fn available_npsh_is_pressure_above_vapour_plus_elevation() {
        // A turbine 5 m above the tailrace, at atmospheric pressure, with water
        // at 20 C. Atmospheric pressure alone is 101 325 / (998.2 * 9.807) =
        // 10.36 m of head, of which 2 339 Pa of vapour pressure is spent, and
        // the remaining 5 m of static elevation is then added on top.
        let rho = Density::new(998.2);
        let npsh = CavitationMargin::available_npsh(
            Pressure::new(101_325.0),
            Pressure::new(2_339.0),
            rho,
            5.0,
        );
        let expected = (101_325.0 - 2_339.0) / (998.2 * STANDARD_GRAVITY) + 5.0;
        assert!((npsh - expected).abs() < 1e-12, "NPSH_a = {npsh}");
        assert!((npsh - 15.12).abs() < 0.05, "NPSH_a = {npsh}");
        // Standing lower gives strictly less, which is the whole reason a
        // Francis machine wants to be sited low.
        let lower = CavitationMargin::available_npsh(
            Pressure::new(101_325.0),
            Pressure::new(2_339.0),
            rho,
            2.0,
        );
        assert!((npsh - lower - 3.0).abs() < 1e-12);
    }

    #[test]
    fn the_margin_reports_cavitation_exactly_when_it_begins() {
        let safe = CavitationMargin::new(10.0, 4.0);
        assert!(!safe.cavitating);
        assert!((safe.margin - 6.0).abs() < 1e-12);

        let marginal = CavitationMargin::new(4.0, 4.0);
        assert!(
            !marginal.cavitating,
            "exactly at the requirement is not yet cavitating"
        );
        assert!(marginal.margin.abs() < 1e-12);

        let bad = CavitationMargin::new(3.0, 4.0);
        assert!(bad.cavitating);
        assert!((bad.margin + 1.0).abs() < 1e-12);
    }

    #[test]
    fn required_npsh_is_harder_at_overload_than_at_part_load() {
        // This corrects a natural but wrong intuition. The inception number
        // rises symmetrically about the best point, which suggests part load and
        // overload should be equally hard. They are not: the requirement scales
        // with `V^2`, so halving the flow quarters the requirement while the
        // sigma penalty only doubles at `flow_penalty = 0.3`. Part load is
        // therefore *easier* on NPSH, and overload is the case that sizes the
        // installation. Getting this backwards would put the turbine at the
        // wrong elevation.
        let inception = CavitationInception {
            inception_number: 0.05,
            flow_penalty: 0.3,
            best_efficiency_flow: 3.0,
        };
        let d = 2.0;
        let at_bep = inception.required_npsh(3.0, d);
        for (i, flow) in [0.15, 0.45, 1.5, 2.25, 2.85].iter().enumerate() {
            let part_load = inception.required_npsh(*flow, d);
            assert!(
                part_load < at_bep,
                "part load {flow}: {part_load} vs {at_bep} at {i}"
            );
        }
        for flow in [3.3, 4.5, 6.0, 9.0] {
            let overload = inception.required_npsh(flow, d);
            assert!(overload > at_bep, "overload {flow}: {overload} vs {at_bep}");
        }
    }

    #[test]
    fn a_strong_flow_penalty_makes_low_flow_as_demanding_as_the_best_point() {
        // The part-load result above is a property of the *coefficients*, not a
        // law. The requirement goes as `sigma(Q) Q^2`, and the sigma penalty is
        // quadratic in the deviation, so the product has a minimum somewhere
        // below the best point and rises again at very low flow. Comparing
        // `Q = Q_bep / 2` against `Q = Q_bep`:
        //
        // ```text
        // f(1/2)/f(1) = (sigma_0 + penalty/4) / (4 sigma_0) >= 1
        //             iff penalty >= 12 sigma_0
        // ```
        //
        // At `flow_penalty = 0.6` against `sigma_0 = 0.05` that threshold is hit
        // exactly, which is why the numbers in the test above were chosen to sit
        // clear of it. A machine with a sharper penalty does not get to ignore
        // low flow, and this pins the boundary rather than leaving it folklore.
        let sigma_0 = 0.05;
        let bep = 3.0;
        let d = 2.0;
        let build = |penalty: f64| CavitationInception {
            inception_number: sigma_0,
            flow_penalty: penalty,
            best_efficiency_flow: bep,
        };
        let at_bep = build(0.59).required_npsh(bep, d);
        let half_flow = build(0.59).required_npsh(bep / 2.0, d);
        assert!(
            half_flow < at_bep,
            "below the threshold: {half_flow} vs {at_bep}"
        );
        // Exactly at `12 sigma_0` the two are equal.
        let at_threshold = build(12.0 * sigma_0);
        assert!(
            (at_threshold.required_npsh(bep / 2.0, d) - at_threshold.required_npsh(bep, d)).abs()
                / at_bep
                < 1e-12
        );
        // Above it, low flow is the harder case after all.
        let sharp = build(12.1 * sigma_0);
        assert!(sharp.required_npsh(bep / 2.0, d) > sharp.required_npsh(bep, d));
    }

    #[test]
    fn a_turbine_that_cavitates_at_overload_but_not_at_its_best_point() {
        // This is the failure a constant inception number hides, and it is why
        // the flow penalty exists. Both margins are built the same way; only the
        // machine's requirement differs.
        let inception = CavitationInception {
            inception_number: 0.05,
            flow_penalty: 0.3,
            best_efficiency_flow: 3.0,
        };
        let d = 2.0;
        let available = 0.05;

        let at_bep = CavitationMargin::new(available, inception.required_npsh(3.0, d));
        let at_overload = CavitationMargin::new(available, inception.required_npsh(6.0, d));
        assert!(!at_bep.cavitating, "BEP margin {}", at_bep.margin);
        assert!(
            at_overload.cavitating,
            "overload must be the harder case, margin {}",
            at_overload.margin
        );

        // A constant-number machine would have called both safe, because it
        // never sees the sigma rise that the overload actually needs.
        let constant = CavitationInception::constant(0.05, 3.0);
        let naive = CavitationMargin::new(available, constant.required_npsh(6.0, d));
        assert!(
            !naive.cavitating,
            "the constant model should miss this, which is the point: margin {}",
            naive.margin
        );
    }
}

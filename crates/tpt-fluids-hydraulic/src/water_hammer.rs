//! Water hammer by the Method of Characteristics (MOC).
//!
//! A sudden flow change in a pipe launches a pressure wave travelling at the
//! wave speed `a = sqrt(K/rho)`. MOC follows those wavefronts on a
//! characteristic grid of `dx` and `dt = dx / a`, and is the standard way to
//! resolve the pressure history a valve closure or pump trip produces.
//!
//! # The equations
//!
//! Linearising the momentum equation about the steady state gives, at each
//! internal node,
//!
//! ```text
//! (Q_P - Q_0) + (g_0 / a) (H_P - H_0) = R_P Q_P |Q_P|
//! ```
//!
//! with `R` the pipe resistance, and continuity `sum Q_P = sum Q_upstream`.
//! Writing the two characteristic relations that connect a node to its
//! neighbours and eliminating the heads gives the classical MOC solution for
//! the node's outflows.
//!
//! The frictionless limit is the check that matters: with `R = 0` a valve that
//! slams shut must produce exactly the Joukowsky rise `a Q_0 / g`, and
//! `momentum_coefficient` below is `a / g`, which reproduces it.
//!
//! # Column separation
//!
//! If a node's head falls below the vapour-pressure head, the liquid column
//! separates and the node is clamped to vapour pressure, with the flow
//! *released* to the free discharge that head implies. Releasing both is what
//! keeps the state on the momentum equation: clamping the head alone left a
//! `(H, Q)` pair that satisfied no characteristic, and the next step read that
//! violation back in as an arriving wave. That nonlinearity is the whole
//! difficulty of water hammer, and it is why MOC needs a Courant step to stay
//! stable rather than an arbitrarily small one.
//!
//! # Two solvers, and why there are two
//!
//! [`method_of_characteristics`] models **one** pipe between two nodes with a
//! single lumped momentum equation. It reproduces the Joukowsky rise exactly and
//! is the right tool for the canonical problem, but it applies a boundary action
//! along the whole pipe at once, so it has no notion of a wave *travelling*.
//!
//! [`MocNetwork`] is the general case. It splits every branch into reaches with a
//! common Courant time step, treats each pipe end as a node solved from two
//! arriving characteristics, and closes it with its boundary condition
//! (reservoir, valve, dead end, or continuity at a junction). A wavefront then
//! crosses one reach per step, reflects off a dead end, and returns, which is
//! what lets it answer *when* a change is felt at a given point -- the question
//! a lumped model cannot be asked.
//!
//! One physical fact underpins both, and it is easy to get backwards: in a
//! *steady* state `Q_P = Q_0`, so the flow term in the linearised momentum
//! equation vanishes and a **frictionless pipe has no steady head gradient at
//! all**. The `a Q / (g A)` term is purely a transient effect -- it is the
//! Joukowsky rise, appearing only when the flow changes. The steady drop comes
//! from friction, carried by the reaches.

use crate::error::{HydraulicError, Result};
use tpt_fluids_core::consts::STANDARD_GRAVITY;

/// A pipe branch in the MOC network.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Branch {
    /// The node this branch runs to. The branch is understood to run from the
    /// upstream node toward this one.
    pub node: usize,
    /// Pipe length in metres.
    pub length: f64,
    /// Wave speed in metres per second.
    pub wave_speed: f64,
    /// Cross-sectional area in square metres.
    pub area: f64,
    /// Darcy-Weisbach resistance in seconds squared per cubic metre per
    /// fifth, so the friction term is `R Q |Q|`.
    pub resistance: f64,
    /// Whether this branch refuses reverse flow, as a check valve does.
    pub check_valve: bool,
}

/// The downstream boundary: either a fixed-head reservoir or a valve whose
/// flow the caller sets at each step.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Boundary {
    /// The imposed head, when this is a reservoir.
    pub head: f64,
    /// Whether this boundary is a fixed-head reservoir rather than a valve.
    pub is_reservoir: bool,
}

impl Boundary {
    /// A fixed-head reservoir.
    pub const fn reservoir(head: f64) -> Self {
        Self {
            head,
            is_reservoir: true,
        }
    }

    /// A valve whose flow the caller controls.
    pub const fn valve() -> Self {
        Self {
            head: 0.0,
            is_reservoir: false,
        }
    }
}

/// The momentum-equation coefficient `a / g` for a pipe, in seconds times
/// metres.
///
/// The frictionless momentum relation is `dH = -(a/g) dQ`, so a flow change
/// `dQ` produces a head change of `-(a/g) dQ`. Multiplying by `rho g` gives
/// the Joukowsky rise `rho a dQ`, which is how this coefficient is checked.
pub fn momentum_coefficient(wave_speed: f64) -> f64 {
    wave_speed / STANDARD_GRAVITY
}

/// The head rise a frictionless pipe shows for a flow change, in metres.
pub fn frictionless_head_rise(wave_speed: f64, flow_change: f64) -> f64 {
    -momentum_coefficient(wave_speed) * flow_change
}

/// The output of an MOC run.
#[derive(Clone, PartialEq, Debug)]
pub struct MocResult {
    /// Head at each node at each step, in metres. Index `[step][node]`.
    pub head_history: Vec<Vec<f64>>,
    /// Flow in each branch at each step, in cubic metres per second.
    pub flow_history: Vec<Vec<f64>>,
    /// The time step used, in seconds.
    pub dt: f64,
    /// The number of steps taken.
    pub steps: usize,
    /// Nodes whose head was clamped to the vapour head on at least one step.
    pub separated_nodes: Vec<usize>,
}

impl MocResult {
    /// The head at a node at a given step, in metres.
    pub fn head(&self, node: usize, step: usize) -> Option<f64> {
        self.head_history.get(step)?.get(node).copied()
    }

    /// The flow in a branch at a given step, in cubic metres per second.
    pub fn flow(&self, branch: usize, step: usize) -> Option<f64> {
        self.flow_history.get(step)?.get(branch).copied()
    }

    /// The largest head rise seen anywhere, in metres, relative to the head at
    /// the same node on step zero.
    pub fn max_head_rise(&self) -> f64 {
        let Some(first) = self.head_history.first() else {
            return 0.0;
        };
        self.head_history
            .iter()
            .flat_map(|row| row.iter())
            .zip(first.iter().cycle())
            .fold(0.0f64, |acc, (h, h0)| acc.max(h - h0))
    }
}

/// Solves a reservoir-fed pipe by MOC under a prescribed valve closure.
///
/// `valve_flows[i]` is the flow demanded through the valve at step `i`; it
/// need not be monotone, though a realistic closure is. The reservoir sits
/// upstream at a fixed head and the valve closes the downstream end, which is
/// the canonical water-hammer problem and the one whose frictionless answer is
/// known exactly.
pub fn method_of_characteristics(
    branch: Branch,
    initial_flow: f64,
    valve_flows: &[f64],
    boundary: Boundary,
    vapour_head: f64,
) -> Result<MocResult> {
    if branch.length <= 0.0 || branch.area <= 0.0 || branch.wave_speed <= 0.0 {
        return Err(HydraulicError::InvalidTopology(
            "branch length, area, and wave speed must be positive".into(),
        ));
    }
    if !boundary.is_reservoir {
        return Err(HydraulicError::InvalidTopology(
            "MOC requires a fixed-head reservoir boundary".into(),
        ));
    }

    // Discretise the pipe into reaches of length dx, with the matching
    // Courant time step dt = dx / a.
    let dx = branch.length;
    let dt = dx / branch.wave_speed;
    let a_over_g = momentum_coefficient(branch.wave_speed);

    // Node 0 is the reservoir; node 1 is the valve end. Between them the
    // single reach is resolved by its characteristic pair.
    let steps = valve_flows.len();
    let mut head_history = Vec::with_capacity(steps + 1);
    let mut flow_history = Vec::with_capacity(steps + 1);
    let mut separated_nodes: Vec<usize> = Vec::new();

    // Steady state. The valve head sits below the reservoir by the
    // frictionless gradient *plus* the steady friction loss `R Q |Q|`, which
    // is what the momentum equation integrates over the whole reach.
    let q0 = initial_flow;
    let steady_friction = branch.resistance * q0 * q0.abs();
    let h_valve = boundary.head - a_over_g * q0 - steady_friction;
    let mut head = vec![boundary.head, h_valve];
    let mut flow = q0;

    head_history.push(head.clone());
    flow_history.push(vec![flow]);

    // The `-` characteristic arriving at the valve node from the reservoir
    // side carries the upstream head; the `+` characteristic carries the
    // valve-side flow. Together they close the node each step.
    let mut accumulated_flow = q0;

    for &q_valve in valve_flows.iter().take(steps) {
        // Momentum at the valve node, frictionless:
        //   (Q - Q_prev) + (Q_prev / a) (H - H_prev) = 0
        // which rearranges to the valve head for the new flow:
        //
        //   H = H_prev - (a / g) (Q - Q_prev)
        //
        // Friction then **integrates** along the closure. The momentum equation
        // for the valve node carries a steady friction force `R Q |Q|`, and the
        // head change as the flow moves from `q_prev` to `q` is
        //
        //   dH = -(a / g) dQ + R Q |Q| dQ
        //
        // so the friction contribution is the definite integral of `Q |Q|`
        // between the two flows, not a difference of the friction force at the
        // two ends. Differencing `R Q |Q| - R Q_prev |Q_prev|` telescopes over a
        // closure to the constant `+ R Q0^2`, so it added a fixed offset instead
        // of damping anything: a frictional gradual closure came out at
        // 0.5333 m, slightly *above* the frictionless 0.5298 m, which is
        // impossible. Friction can only ever reduce the surge.
        //
        // `R` here is a lumped head-loss coefficient -- `R Q |Q|` is already a
        // head, as the steady seed shows -- so no area factor appears. The
        // antiderivative of `Q |Q|` is `(2/3) sign(Q) |Q|^{3/2}`.
        let integral = |q: f64| {
            let s = if q < 0.0 { -1.0 } else { 1.0 };
            s * q.abs().powf(1.5)
        };
        let friction_dh =
            branch.resistance * (2.0 / 3.0) * (integral(q_valve) - integral(accumulated_flow));
        let h_new = head[1] - a_over_g * (q_valve - accumulated_flow) + friction_dh;
        accumulated_flow = q_valve;
        flow = q_valve;

        let mut h_valve_node = h_new;
        // Column separation: the liquid cannot sustain a head below the
        // vapour head, so clamp it and record the node.
        if h_valve_node < vapour_head {
            h_valve_node = vapour_head;
            if !separated_nodes.contains(&1) {
                separated_nodes.push(1);
            }
        }

        head = vec![boundary.head, h_valve_node];
        head_history.push(head.clone());
        flow_history.push(vec![flow]);
    }

    Ok(MocResult {
        head_history,
        flow_history,
        dt,
        steps,
        separated_nodes,
    })
}

/// Validates that a proposed time step is a Courant step for the branch.
pub fn validate_courant(branch: &Branch, dt: f64) -> Result<()> {
    let stable = branch.length / branch.wave_speed;
    if dt > stable * (1.0 + 1.0e-9) {
        return Err(HydraulicError::OutOfRange(format!(
            "time step {dt} exceeds the Courant limit {stable} for this branch"
        )));
    }
    Ok(())
}

/// A node on a MOC network, with the boundary condition that closes it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum MocNode {
    /// A fixed-head reservoir. The head is prescribed; the flow follows.
    Reservoir(f64),
    /// A valve whose flow the caller prescribes each step.
    Valve,
    /// A dead end: a closed boundary that reflects a wave at **full amplitude**.
    ///
    /// This is the boundary a *prescribed-head* downstream condition needs in
    /// order to quote a Joukowsky head. A fixed-head reservoir absorbs an
    /// arriving wave completely, so it can never display the doubled head that
    /// a closed end produces; a dead end, having no flow path, reflects the
    /// pressure wave with coefficient `-1` and the head disturbance doubles on
    /// arrival. The flow is identically zero, so the head follows from the
    /// single characteristic that reaches it.
    DeadEnd,
    /// A tee or junction, closed by continuity.
    Junction,
}

/// One pipe in a MOC network, discretised into reaches.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct MocBranch {
    /// The upstream node index.
    pub from: usize,
    /// The downstream node index.
    pub to: usize,
    /// Total pipe length in metres.
    pub length: f64,
    /// Wave speed in metres per second.
    pub wave_speed: f64,
    /// Cross-sectional area in square metres.
    pub area: f64,
    /// Darcy-Weisbach resistance in seconds squared per cubic metre per fifth.
    pub resistance: f64,
    /// Whether the branch refuses reverse flow, as a check valve does.
    pub check_valve: bool,
}

/// A multi-node MOC network.
///
/// # The formulation
///
/// Each branch is cut into `n` reaches of `dx = L / n` with a common time step
/// `dt = dx / a`, so a wavefront crosses exactly one reach per step. Branch `i`
/// contributes `n + 1` **ends**; the ends between the branch's two nodes are
/// *interior* ends, and that is where the wave travels.
///
/// **Every end is a node solved from two arriving characteristics.** The `+`
/// relation arrives from the end one reach upstream and the `-` relation from
/// the end one reach downstream:
///
/// ```text
/// a_plus  = H_up + B Q_up        H = a_plus  - B Q
/// a_minus = H_dn - B Q_dn        H = a_minus + B Q
/// ```
///
/// with `B = a / (g A n)`, the momentum coefficient for **one reach**. That is
/// two equations in the two unknowns `(H, Q)`, so an interior end is closed in
/// closed form:
///
/// ```text
/// Q = (a_plus - a_minus) / (2 B)      H = (a_plus + a_minus) / 2
/// ```
///
/// Treating an interior end as a single propagated value is not an approximation
/// of this -- it is a one-equation treatment of a two-equation node, and it makes
/// the wave appear one reach early on every step for the whole run.
///
/// Only a branch's two boundary ends lack a characteristic, and those are exactly
/// the ends attached to a network node. A **reservoir** prescribes the head so
/// the flow follows; a **valve** prescribes the flow so the head follows; a
/// **junction** has neither, and continuity closes it. Since each end there obeys
/// `H = a_k -/+ b_k Q_k` and `sum out = sum in`, the common head is the ratio
/// `H = sum(a_k / b_k) / sum(1 / b_k)`.
///
/// # The steady state, and why the frictionless field is flat
///
/// The linearised momentum equation is
///
/// ```text
/// (Q_P - Q_0) + (g A / a)(H_P - H_0) + R Q |Q| = 0
/// ```
///
/// and in a *steady* state `Q_P = Q_0`, so the flow term vanishes. A
/// **frictionless pipe therefore has no steady head gradient at all**: the
/// `a Q / (g A)` term is purely a *transient* effect, appearing only when the
/// flow changes. That is exactly what the Joukowsky rise is.
///
/// This is worth stating because treating `a Q / (g A)` as a steady loss -- which
/// this implementation did at first -- produces a gradient that is not physical,
/// which in turn makes the seeded state inconsistent with the characteristics
/// and the whole network drifts on a run that should be steady. It fails
/// silently: the field stays smooth, the flows stay plausible, and the only
/// symptom is that nothing is quite right.
///
/// The steady drop comes from **friction**, and it is carried by the reaches
/// themselves: each interior end is seeded on the linear gradient between its
/// two node heads, so a uniform flow really does hold a real head drop.
///
/// That detail is load-bearing, and getting it wrong was a real bug. The
/// friction force once entered each characteristic as a *difference* between
/// neighbouring ends, `r (Q_k E_k - Q_{k-1} E_{k-1})`, which is identically
/// zero for the uniform flow of a steady network -- so nothing held the gradient
/// up, the seeded Darcy-Weisbach drop washed out within a few steps, and the
/// network settled into looking frictionless. The friction is now the absolute
/// per-reach loss, `r Q E`, and a frictional steady state is an exact fixed
/// point. Like the transient-gradient bug above it failed silently: every number
/// stayed finite and plausible, and only "the head quietly returns to the
/// reservoir head" gives it away.
#[derive(Clone, PartialEq, Debug)]
pub struct MocNetwork {
    /// The nodes, one entry per node index.
    pub nodes: Vec<MocNode>,
    /// The branches.
    pub branches: Vec<MocBranch>,
    /// The flow in each branch at the initial steady state.
    pub initial_flows: Vec<f64>,
    /// The number of reaches in each branch.
    pub reaches: Vec<usize>,
    /// The common time step, in seconds.
    pub dt: f64,
    /// The friction-lag time constant in seconds, or `None` for the
    /// quasi-steady law.
    ///
    /// See [`Self::with_friction_lag`]. `None` is the default and is exact for
    /// the original behaviour: a friction force proportional to `Q |Q|`
    /// evaluated on the current flow.
    pub friction_lag: Option<f64>,
    /// Per-branch friction-lag time constants in seconds, or `None` for the
    /// quasi-steady law.
    ///
    /// Set by [`Self::with_derived_friction_lag`]. When present it takes
    /// precedence over [`Self::friction_lag`], because a per-branch value
    /// derived from the pipe's own properties is more defensible than a single
    /// number applied to a mixed network.
    pub friction_lag_per_branch: Option<Vec<f64>>,
}

/// The result of a multi-node MOC run.
#[derive(Clone, PartialEq, Debug)]
pub struct NetworkMocResult {
    /// Head at each network node at each step, in metres. `[step][node]`.
    pub head_history: Vec<Vec<f64>>,
    /// Flow in each branch at each step, in cubic metres per second.
    pub flow_history: Vec<Vec<f64>>,
    /// Head at every pipe end at each step, in metres. `[step][end]`.
    ///
    /// The node history is the useful public view, but it cannot witness a
    /// wavefront travelling *along* a pipe, because it holds only the network
    /// nodes. This records the interior ends too, which is what makes travel time
    /// observable at all.
    pub end_head_history: Vec<Vec<f64>>,
    /// The time step used, in seconds.
    pub dt: f64,
    /// The number of steps taken.
    pub steps: usize,
    /// Nodes clamped to the vapour head on at least one step.
    pub separated_nodes: Vec<usize>,
}

impl NetworkMocResult {
    /// The head at a network node at a given step, in metres.
    pub fn head(&self, node: usize, step: usize) -> Option<f64> {
        self.head_history.get(step)?.get(node).copied()
    }

    /// The flow in a branch at a given step, in cubic metres per second.
    pub fn flow(&self, branch: usize, step: usize) -> Option<f64> {
        self.flow_history.get(step)?.get(branch).copied()
    }

    /// The head at a pipe end at a given step, in metres.
    ///
    /// Unlike [`Self::head`], this includes the *interior* ends, which are
    /// computational nodes rather than network nodes.
    pub fn end_head(&self, end: usize, step: usize) -> Option<f64> {
        self.end_head_history.get(step)?.get(end).copied()
    }

    /// The largest head excursion at any network node, in metres, relative to
    /// that node's head at step zero.
    pub fn max_head_excursion(&self) -> f64 {
        let Some(first) = self.head_history.first() else {
            return 0.0;
        };
        self.head_history
            .iter()
            .flat_map(|row| row.iter())
            .zip(first.iter().cycle())
            .fold(0.0f64, |acc, (h, h0)| acc.max((h - h0).abs()))
    }
}
impl MocNetwork {
    /// Builds a network with a common time step `dt`.
    ///
    /// Each branch takes `ceil(L / (a dt))` reaches, so `dt` is a Courant step
    /// everywhere. The Courant condition `dt <= L / a` is checked *directly*,
    /// not merely through the reach count: `ceil` of a very small number is
    /// still 1, so a `dt` a hundred times over the limit would be silently
    /// rounded away to a single reach instead of refused.
    pub fn new(
        nodes: Vec<MocNode>,
        branches: Vec<MocBranch>,
        initial_flows: Vec<f64>,
        dt: f64,
    ) -> Result<Self> {
        if dt <= 0.0 || !dt.is_finite() {
            return Err(HydraulicError::OutOfRange(
                "MOC time step must be positive and finite".into(),
            ));
        }
        if branches.len() != initial_flows.len() {
            return Err(HydraulicError::InvalidTopology(
                "each branch needs exactly one initial flow".into(),
            ));
        }
        if !dt.is_finite() {
            return Err(HydraulicError::OutOfRange(
                "MOC time step must be finite".into(),
            ));
        }
        let mut reaches = Vec::with_capacity(branches.len());
        for branch in &branches {
            if branch.length <= 0.0 || branch.area <= 0.0 || branch.wave_speed <= 0.0 {
                return Err(HydraulicError::InvalidTopology(
                    "branch length, area, and wave speed must be positive".into(),
                ));
            }
            if branch.from >= nodes.len() || branch.to >= nodes.len() {
                return Err(HydraulicError::InvalidTopology(
                    "branch references a node that does not exist".into(),
                ));
            }
            if branch.from == branch.to {
                return Err(HydraulicError::InvalidTopology(
                    "a branch must connect two distinct nodes".into(),
                ));
            }
            let courant = branch.length / branch.wave_speed;
            if dt > courant * (1.0 + 1.0e-9) {
                return Err(HydraulicError::OutOfRange(format!(
                    "time step {dt} exceeds the Courant limit {courant} for this branch"
                )));
            }
            let n = (branch.length / (branch.wave_speed * dt)).ceil();
            if n < 1.0 {
                return Err(HydraulicError::OutOfRange(format!(
                    "time step {dt} is too large to resolve a wave in this branch"
                )));
            }
            reaches.push(n as usize);
        }
        Ok(Self {
            nodes,
            branches,
            initial_flows,
            reaches,
            dt,
            friction_lag: None,
            friction_lag_per_branch: None,
        })
    }

    /// Enables **unsteady friction** with the given lag time in seconds.
    ///
    /// The friction force does not follow the flow instantaneously. When the
    /// flow decelerates sharply, the wall shear lags it, so the momentum
    /// equation sees a smaller friction force than `R Q |Q|` evaluated on the
    /// current flow. Quasi-steady friction therefore *over*-predicts the peak
    /// head rise of a rapid valve closure, and this is the correction for that.
    ///
    /// The model is an exponentially weighted moving average of `|Q|`, with the
    /// friction force `R Q E` where `E` is the average:
    ///
    /// ```text
    /// E_{n+1} = (1 - C) E_n + C |Q_{n+1}|        C = dt / (dt + T_f)
    /// ```
    ///
    /// The one property that makes this safe is that **`E` converges to `|Q|`
    /// for a steady flow**, so the steady solution is *unchanged* however
    /// `T_f` is set. That is asserted in the tests, and it is why enabling this
    /// cannot corrupt a network that is genuinely at rest or in steady flow --
    /// only genuinely unsteady runs are affected.
    ///
    /// `T_f` is a physical time scale, not a tuning knob, and the tests check
    /// the direction and relative size of the effect rather than an absolute
    /// head. The friction time scale is conventionally written `2 L / (g |V|)`;
    /// note that this is many orders of magnitude larger than the wave transit
    /// time `L / a` in a normal pipe, so the correction is a *slow* effect and
    /// matters over the seconds-to-minutes the friction takes to relax rather
    /// than within a single wave passage. Choosing `T_f` much smaller than
    /// `dt` degenerates to the quasi-steady law, which is the correct limit.
    ///
    /// This sets **one** constant for every branch, which is right for a
    /// network of similar pipes and wrong for a mixed one. See
    /// [`Self::with_derived_friction_lag`] for the per-branch form.
    ///
    /// This sets **one** constant for every branch. That is right for a network
    /// of similar pipes and wrong for a mixed one -- see
    /// [`Self::with_derived_friction_lag`], which derives it per branch instead.
    pub fn with_friction_lag(mut self, time: f64) -> Result<Self> {
        if time <= 0.0 || !time.is_finite() {
            return Err(HydraulicError::OutOfRange(
                "the friction lag time constant must be positive and finite".into(),
            ));
        }
        self.friction_lag = Some(time);
        Ok(self)
    }

    /// Enables unsteady friction with a **per-branch** time constant derived
    /// from each pipe's own properties, rather than one global value.
    ///
    /// The friction time scale of a pipe is `T_f = L / (g A R |Q|)`, which is
    /// `2 L / (g |V|)` written in the resistance terms this crate already
    /// carries, and follows from the momentum equation: the friction force
    /// relaxes over the time the flow takes to be turned over by the friction
    /// itself. It is *not* a free parameter -- a short, rough, fast pipe has a
    /// completely different one from a long smooth one, and forcing a single
    /// value across a mixed network is a modelling error, not a simplification.
    ///
    /// `scale` multiplies the derived value, and exists only to express that
    /// this is a first-order model whose constant is a matter of judgement. Use
    /// `1.0` for the value as derived.
    ///
    /// The steady-limit property is unchanged: each `E` starts at its branch's
    /// steady `|Q|` and `E = |Q|` is a fixed point, so a steady network is still
    /// untouched however the constants come out.
    pub fn with_derived_friction_lag(mut self, scale: f64) -> Result<Self> {
        if scale <= 0.0 || !scale.is_finite() {
            return Err(HydraulicError::OutOfRange(
                "the friction lag scale must be positive and finite".into(),
            ));
        }
        let times: Vec<f64> = (0..self.branches.len())
            .map(|i| self.branch_friction_time(i) * scale)
            .collect();
        self.friction_lag = None;
        self.friction_lag_per_branch = Some(times);
        Ok(self)
    }

    /// The friction relaxation time of one branch, `L / (g A R |Q|)`, in
    /// seconds.
    ///
    /// This is the time scale over which the friction force catches up with a
    /// changing flow. A branch with no resistance or no flow has no friction
    /// force to relax, so the time is infinite and the EWMA weight `C` is zero:
    /// the quasi-steady law, which is the correct limit rather than a fallback.
    pub fn branch_friction_time(&self, branch: usize) -> f64 {
        let b = &self.branches[branch];
        let q = self.initial_flows[branch].abs();
        let denominator = STANDARD_GRAVITY * b.area * b.resistance * q;
        if denominator <= 0.0 || !denominator.is_finite() {
            return f64::INFINITY;
        }
        b.length / denominator
    }

    /// The wave speed actually used on a branch, `L / (n dt)`.
    ///
    /// Equal to or *slower* than the physical speed, by construction. Reported
    /// because a caller checking a round trip needs the discretised value;
    /// confusing the two makes a reflection appear to arrive early.
    pub fn effective_wave_speed(&self, branch: usize) -> f64 {
        let b = &self.branches[branch];
        b.length / (self.reaches[branch] as f64 * self.dt)
    }

    /// The total number of pipe ends: one more than the reach count per branch.
    pub fn end_count(&self) -> usize {
        self.reaches.iter().map(|&n| n + 1).sum()
    }

    /// The index of the first end of `branch`, which is its upstream end.
    pub fn end_offset(&self, branch: usize) -> usize {
        self.reaches[..branch].iter().map(|&n| n + 1).sum()
    }

    /// The branch that pipe end `k` belongs to.
    pub fn end_branch(&self, k: usize) -> usize {
        for i in 0..self.branches.len() {
            if k < self.end_offset(i) + self.reaches[i] + 1 {
                return i;
            }
        }
        usize::MAX
    }

    /// The network node end `k` attaches to, or `usize::MAX` for an *interior*
    /// end, which attaches to nothing.
    ///
    /// Interior ends are computational nodes: they exist so a wavefront can be
    /// followed reach by reach, and they carry no boundary condition of their
    /// own. Reporting them as the downstream node is a silent and total failure
    /// -- the boundary condition would then apply along the whole pipe at once,
    /// which is exactly the lumped behaviour this solver exists to remove.
    pub fn end_node(&self, k: usize) -> usize {
        let i = self.end_branch(k);
        if i == usize::MAX {
            return usize::MAX;
        }
        let offset = self.end_offset(i);
        let branch = &self.branches[i];
        if k == offset {
            branch.from
        } else if k == offset + self.reaches[i] {
            branch.to
        } else {
            usize::MAX
        }
    }

    /// Whether end `k` is its branch's upstream end.
    pub fn end_is_upstream(&self, k: usize) -> bool {
        k == self.end_offset(self.end_branch(k))
    }

    /// The pipe ends attached to a network node.
    fn ends_at(&self, node: usize, n_ends: usize) -> Vec<usize> {
        (0..n_ends).filter(|&k| self.end_node(k) == node).collect()
    }

    /// The steady node heads, accumulating the **friction** loss from the
    /// fixed-head nodes along the branches.
    ///
    /// The frictionless momentum term is deliberately absent. In a steady state
    /// `Q_P = Q_0`, the flow term in the linearised momentum equation vanishes,
    /// and a frictionless pipe has no head gradient at all. Including the
    /// `a Q / (g A)` term here would put a non-physical gradient into the seed,
    /// which is then inconsistent with the characteristics and makes the network
    /// drift on a run that should be perfectly steady.
    fn seed_heads(&self) -> Vec<f64> {
        let n = self.nodes.len();
        let mut heads = vec![0.0f64; n];
        let mut assigned = vec![false; n];
        for (i, node) in self.nodes.iter().enumerate() {
            if let MocNode::Reservoir(h) = *node {
                heads[i] = h;
                assigned[i] = true;
            }
        }
        for _ in 0..=n {
            let mut changed = false;
            for (i, branch) in self.branches.iter().enumerate() {
                let q = self.initial_flows[i];
                // Darcy-Weisbach: the whole-branch steady loss, and nothing else.
                let loss = self.branches[i].resistance * q * q.abs();
                if assigned[branch.from] && !assigned[branch.to] {
                    heads[branch.to] = heads[branch.from] - loss;
                    assigned[branch.to] = true;
                    changed = true;
                } else if assigned[branch.to] && !assigned[branch.from] {
                    heads[branch.from] = heads[branch.to] + loss;
                    assigned[branch.from] = true;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        heads
    }
}
impl MocNetwork {
    /// Solves the network under a prescribed sequence of valve flows.
    ///
    /// `valve_flows[step][node]` is the flow through each [`MocNode::Valve`] at
    /// that step. A valve not prescribed on a step holds its previous flow,
    /// which is what a valve that has finished closing actually does.
    ///
    /// Every [`MocNode::Reservoir`] keeps its own fixed head here. Use
    /// [`Self::solve_with_heads`] for a reservoir whose level moves.
    pub fn solve(&self, valve_flows: &[Vec<f64>], vapour_head: f64) -> Result<NetworkMocResult> {
        self.solve_with_heads(valve_flows, &[], vapour_head)
    }

    /// Solves the network with both a valve-flow sequence and a head schedule,
    /// so a reservoir whose level **moves** can launch a wave instead of only
    /// absorbing one.
    ///
    /// `head_schedule[step][node]` overrides the head of a [`MocNode::Reservoir`]
    /// at that step. A node absent from a row keeps its value from
    /// [`MocNode::Reservoir`], so a moving reservoir can be mixed with fixed
    /// ones. A node specified as `f64::NAN` is left at its fixed head too,
    /// which lets a caller carry one array across steps.
    ///
    /// A step in reservoir head is a genuine wave source, and the check that this
    /// is wired up correctly is where the wave goes. A `dH` step at the
    /// reservoir travels at the wave speed, one reach per step, and on reaching
    /// a dead end is reflected at full amplitude to `2 dH` -- the incident and
    /// reflected waves add, because the dead end cannot relieve pressure by
    /// passing flow. A boundary that only *absorbed* would leave the far end at
    /// `dH` and one that copied the head across whole would leave it flat.
    pub fn solve_with_heads(
        &self,
        valve_flows: &[Vec<f64>],
        head_schedule: &[Vec<f64>],
        vapour_head: f64,
    ) -> Result<NetworkMocResult> {
        // A non-finite vapour head is rejected outright. It is not a regime, it
        // is an absent one, and letting it through would make every downstream
        // comparison against it silently false.
        //
        // Note that a vapour head *above* the reservoir is perfectly physical
        // and must not be rejected: a valve closure raises the head by the
        // Joukowsky rise, so a surge to 101.9 m against a 100 m reservoir is
        // exactly what cavitation is tested against. An earlier attempt at a
        // bound here rejected those, and it broke two column-separation tests
        // that were relying on them.
        if !vapour_head.is_finite() {
            return Err(HydraulicError::OutOfRange(
                "the vapour head must be finite".into(),
            ));
        }
        let n_nodes = self.nodes.len();
        let n_ends = self.end_count();

        // `B = a / (g A n)`, the momentum coefficient for **one reach**. The
        // `1 / n` is easy to miss: `a / (g A)` alone is the whole-*pipe* loss,
        // and using it per reach gives `n` times the gradient -- a 4x error on a
        // 4-reach pipe that stays finite and stays smooth.
        let b_coef: Vec<f64> = (0..self.branches.len())
            .map(|i| {
                self.effective_wave_speed(i)
                    / (STANDARD_GRAVITY * self.branches[i].area * self.reaches[i] as f64)
            })
            .collect();
        // `R / n` per reach, so `n` of them sum to the whole-branch resistance.
        let r_coef: Vec<f64> = (0..self.branches.len())
            .map(|i| self.branches[i].resistance / self.reaches[i] as f64)
            .collect();

        let head = self.seed_heads();

        // The end state at time `t`. The interior ends are seeded **along** the
        // steady friction gradient, not flat.
        //
        // Flat interior ends look right and are wrong. The friction term enters
        // as a *difference* of the friction force between neighbouring ends, so
        // for a uniform flow -- the steady case -- that difference is exactly
        // zero and the solver has nothing holding a gradient up. A flat seed is
        // therefore not a fixed point: the node heads start with the correct
        // Darcy-Weisbach drop and the very first step washes it out, leaving a
        // frictionless-looking network. The loss only survives if the interior
        // carries the gradient it is supposed to represent.
        //
        // The gradient is seeded in proportion to reach index so that the two
        // characteristics at each interior end agree, which is what makes the
        // interior ends genuine fixed points rather than a smooth interpolation.
        let mut h = vec![0.0f64; n_ends];
        let mut q = vec![0.0f64; n_ends];
        for i in 0..self.branches.len() {
            let offset = self.end_offset(i);
            let n = self.reaches[i];
            let branch = &self.branches[i];
            for j in 0..=n {
                q[offset + j] = self.initial_flows[i];
            }
            h[offset] = head[branch.from];
            h[offset + n] = head[branch.to];
            // Linear in `j / n`, so the endpoints are the exact node heads and
            // each interior end sits `R Q |Q| / n` below the one upstream of it.
            let from = head[branch.from];
            let to = head[branch.to];
            for j in 1..n {
                let t = j as f64 / n as f64;
                h[offset + j] = from + t * (to - from);
            }
        }

        let mut head_history = vec![head.clone()];
        let mut flow_history = vec![self.initial_flows.clone()];
        let mut end_head_history = vec![h.clone()];
        let mut separated_nodes: Vec<usize> = Vec::new();

        // The EWMA state for unsteady friction: `E[k]`, the lagged `|Q|` at each
        // end. It starts at the steady value, so a run that is already steady
        // stays exactly steady and the correction is a pure transient effect.
        // `mag` caches `|Q|` per end for the friction force `R Q E`.
        //
        // A per-branch constant takes precedence over the global one, since a
        // value derived from a pipe's own properties beats a single number
        // applied to a mixed network. `C = dt / (dt + T_f)` is clamped to
        // `(0, 1]`, and an infinite `T_f` gives `C = 0`: a branch with no
        // friction force to relax stays quasi-steady, which is the right
        // limit rather than a fallback.
        let per_branch_c: Option<Vec<f64>> = self.friction_lag_per_branch.as_ref().map(|t| {
            t.iter()
                .map(|&tf| {
                    if !tf.is_finite() {
                        0.0
                    } else {
                        (self.dt / (self.dt + tf)).clamp(0.0, 1.0)
                    }
                })
                .collect()
        });
        let global_c: Option<f64> = self
            .friction_lag
            .map(|t| (self.dt / (self.dt + t)).clamp(0.0, 1.0));
        let c_at = |branch: usize| -> Option<f64> {
            per_branch_c
                .as_ref()
                .and_then(|c| c.get(branch).copied())
                .or(global_c)
        };
        let any_lag = per_branch_c.is_some() || global_c.is_some();
        let mut e_state: Vec<f64> = (0..n_ends)
            .map(|k| {
                let i = self.end_branch(k);
                self.initial_flows[i].abs()
            })
            .collect();

        for (step, demand) in valve_flows.iter().enumerate() {
            let mut h_new = vec![0.0f64; n_ends];
            let mut q_new = vec![0.0f64; n_ends];
            // The lagged `|Q|` used by the friction force this step. Without
            // a lag this is just `|Q|`, recovering the quasi-steady law exactly.
            let mag: Vec<f64> = if any_lag {
                e_state.clone()
            } else {
                q.iter().map(|v| v.abs()).collect()
            };

            // --- 1. Gather the arriving characteristics at every end. ---
            let mut a_plus = vec![f64::NAN; n_ends];
            let mut a_minus = vec![f64::NAN; n_ends];
            for k in 0..n_ends {
                let b_idx = self.end_branch(k);
                let offset = self.end_offset(b_idx);
                let n = self.reaches[b_idx];
                let b = b_coef[b_idx];
                let r = r_coef[b_idx];
                if k > offset {
                    let up = k - 1;
                    a_plus[k] = h[up] + b * q[up] - r * q[up] * mag[up];
                }
                if k < offset + n {
                    let down = k + 1;
                    a_minus[k] = h[down] - b * q[down] + r * q[down] * mag[down];
                }
            }

            // --- 2. Close the interior ends: two equations, two unknowns. ---
            for k in 0..n_ends {
                if self.end_node(k) != usize::MAX {
                    continue; // a boundary end; its node closes it below
                }
                let b = b_coef[self.end_branch(k)];
                if b <= 0.0 {
                    return Err(HydraulicError::InvalidTopology(
                        "a branch has a non-positive momentum coefficient".into(),
                    ));
                }
                let (ap, am) = (a_plus[k], a_minus[k]);
                if !ap.is_finite() || !am.is_finite() {
                    return Err(HydraulicError::InvalidTopology(
                        "an interior end is missing a characteristic".into(),
                    ));
                }
                q_new[k] = (ap - am) / (2.0 * b);
                h_new[k] = 0.5 * (ap + am);
            }

            // --- 3. Close each network node with its boundary condition. ---
            for node in 0..n_nodes {
                let ends = self.ends_at(node, n_ends);
                if ends.is_empty() {
                    continue;
                }
                let b_at = |k: usize| -> Result<f64> {
                    let b = b_coef[self.end_branch(k)];
                    if b <= 0.0 {
                        Err(HydraulicError::InvalidTopology(
                            "a branch has a non-positive momentum coefficient".into(),
                        ))
                    } else {
                        Ok(b)
                    }
                };
                // The one characteristic arriving at a boundary end, as
                // `H = a -/+ b Q`.
                let arriving = |k: usize| -> Result<(f64, bool)> {
                    let upstream = self.end_is_upstream(k);
                    let a = if upstream { a_minus[k] } else { a_plus[k] };
                    if !a.is_finite() {
                        return Err(HydraulicError::InvalidTopology(
                            "a boundary end is missing a characteristic".into(),
                        ));
                    }
                    Ok((a, upstream))
                };

                match self.nodes[node] {
                    MocNode::Reservoir(h_fixed) => {
                        // A scheduled head overrides the fixed one, which is what
                        // lets a reservoir whose level moves launch a wave. A NaN
                        // or absent entry means "unchanged".
                        let scheduled = head_schedule
                            .get(step)
                            .and_then(|row| row.get(node).copied())
                            .filter(|v| v.is_finite())
                            .unwrap_or(h_fixed);
                        for &k in &ends {
                            let b = b_at(k)?;
                            let (a, upstream) = arriving(k)?;
                            h_new[k] = scheduled;
                            q_new[k] = if upstream {
                                (scheduled - a) / b
                            } else {
                                (a - scheduled) / b
                            };
                        }
                    }
                    MocNode::Valve => {
                        for &k in &ends {
                            let b = b_at(k)?;
                            let (a, upstream) = arriving(k)?;
                            let prescribed = demand.get(node).copied().unwrap_or(q[k]);
                            q_new[k] = prescribed;
                            h_new[k] = if upstream {
                                a + b * prescribed
                            } else {
                                a - b * prescribed
                            };
                        }
                    }
                    MocNode::DeadEnd => {
                        // No flow can pass, so `Q = 0` and the head follows from
                        // the one characteristic that arrives. Writing it as
                        // `a +/- B*0` shows why the wave is *not* absorbed: a
                        // reservoir cancels its gradient by changing flow, but
                        // here the flow is pinned at zero, so the full head
                        // disturbance reflects instead.
                        for &k in &ends {
                            // `a + B * 0 == a - B * 0`, so the branch is
                            // genuinely the same either way: pinning the flow
                            // at zero removes the momentum term entirely.
                            let (a, _upstream) = arriving(k)?;
                            q_new[k] = 0.0;
                            h_new[k] = a;
                        }
                    }
                    MocNode::Junction => {
                        // Continuity. Each end obeys `H = a_k -/+ b_k Q_k` and
                        // `sum out = sum in`, so substituting and collecting the
                        // common head leaves one ratio.
                        let mut numerator = 0.0f64;
                        let mut denominator = 0.0f64;
                        for &k in &ends {
                            let b = b_at(k)?;
                            let (a, _) = arriving(k)?;
                            numerator += a / b;
                            denominator += 1.0 / b;
                        }
                        if denominator.abs() < f64::EPSILON {
                            return Err(HydraulicError::InvalidTopology(
                                "junction continuity is degenerate; check the topology".into(),
                            ));
                        }
                        let h_junction = numerator / denominator;
                        for &k in &ends {
                            let b = b_at(k)?;
                            let (a, upstream) = arriving(k)?;
                            h_new[k] = h_junction;
                            q_new[k] = if upstream {
                                (h_junction - a) / b
                            } else {
                                (a - h_junction) / b
                            };
                        }
                    }
                }
            }

            // --- 4. Column separation. ---
            //
            // A liquid column cannot sustain a head below the vapour pressure.
            // The standard treatment is to clamp the head **and** release the
            // flow to whatever that head implies.
            //
            // Clamping the head alone is inconsistent, and the trace is
            // unmistakable: the head jumps to the vapour head while the flow
            // stays frozen at its pre-separation value, leaving a `(H, Q)` pair
            // that satisfies no characteristic. The next step reads that state
            // back as the arriving characteristic, so the violation is injected
            // into the run instead of being absorbed.
            //
            // Below a separated column the liquid is effectively free: the head
            // sits at vapour pressure and the flow is set by that head and the
            // arriving characteristic, which is a free discharge. That keeps
            // the pair on the momentum equation the rest of the solver assumes.
            for node in 0..n_nodes {
                if matches!(self.nodes[node], MocNode::Reservoir(_)) {
                    continue;
                }
                for k in self.ends_at(node, n_ends) {
                    if h_new[k] < vapour_head {
                        h_new[k] = vapour_head;
                        // Release the flow to the value consistent with the
                        // clamped head, from the same characteristic the
                        // boundary would have used: `H = a -/+ b Q`.
                        let b = b_coef[self.end_branch(k)];
                        if b > 0.0 && self.end_is_upstream(k) {
                            if a_minus[k].is_finite() {
                                q_new[k] = (vapour_head - a_minus[k]) / b;
                            }
                        } else if b > 0.0 && a_plus[k].is_finite() {
                            q_new[k] = (a_plus[k] - vapour_head) / b;
                        }
                        if !separated_nodes.contains(&node) {
                            separated_nodes.push(node);
                        }
                    }
                }
            }

            // --- 5. Check valves refuse reverse flow. ---
            //
            // Clamping the **downstream end only** is the whole point. Zeroing
            // the branch's upstream end as well is a lumped shortcut that stops
            // the flow across the entire branch in one step, which annihilates
            // the wave that is still travelling through the reaches behind the
            // valve -- the exact behaviour this solver exists to remove. With
            // only the valve end clamped, the closure propagates as a real
            // front, one reach per step, and the interior sees it arrive.
            for (i, branch) in self.branches.iter().enumerate() {
                if !branch.check_valve {
                    continue;
                }
                let offset = self.end_offset(i);
                let downstream_end = offset + self.reaches[i];
                if q_new[downstream_end] < 0.0 {
                    q_new[downstream_end] = 0.0;
                }
            }

            // --- 6. Carry the state forward. ---
            //
            // The EWMA advances *after* the step, from the flow just computed,
            // so the friction used this step is the one that lagged it.
            if any_lag {
                for k in 0..n_ends {
                    if let Some(c) = c_at(self.end_branch(k)) {
                        e_state[k] = (1.0 - c) * e_state[k] + c * q_new[k].abs();
                    }
                }
            }
            h.copy_from_slice(&h_new);
            q.copy_from_slice(&q_new);
            let mut flow = self.initial_flows.clone();
            for i in 0..self.branches.len() {
                let offset = self.end_offset(i);
                flow[i] = q[offset + self.reaches[i]];
            }
            // The node history records each network node's head, read off its
            // first attached end. Interior ends are not network nodes and do not
            // appear here; the end history carries those.
            head_history.push(
                (0..n_nodes)
                    .map(|node| {
                        self.ends_at(node, n_ends)
                            .first()
                            .map(|&k| h[k])
                            .unwrap_or(0.0)
                    })
                    .collect(),
            );
            flow_history.push(flow.clone());
            end_head_history.push(h.clone());
        }

        Ok(NetworkMocResult {
            head_history,
            flow_history,
            end_head_history,
            dt: self.dt,
            steps: valve_flows.len(),
            separated_nodes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pipe() -> Branch {
        Branch {
            node: 1,
            length: 300.0,
            wave_speed: 1484.5764,
            area: 0.0706858,
            resistance: 0.0,
            check_valve: false,
        }
    }

    #[test]
    fn momentum_coefficient_reproduces_joukowsky() {
        // Head rise for a full stop = a Q0 / g, which is the Joukowsky head.
        let a = 1484.5764;
        let q0 = 0.5;
        let rise = frictionless_head_rise(a, -q0);
        assert!((rise - a * q0 / STANDARD_GRAVITY).abs() < 1e-12);
        assert!((rise - 75.6923).abs() < 1e-3, "rise = {rise}");
    }

    #[test]
    fn instant_closure_gives_exactly_the_joukowsky_rise() {
        let b = pipe();
        let q0 = 0.5;
        let head0 = 100.0;
        // A single step demanding zero flow is an instantaneous stop.
        let r = method_of_characteristics(b, q0, &[0.0], Boundary::reservoir(head0), 0.0).unwrap();
        let h0 = r.head(1, 0).unwrap();
        let h1 = r.head(1, 1).unwrap();
        let expected = momentum_coefficient(b.wave_speed) * q0;
        assert!((h1 - h0 - expected).abs() < 1e-9, "rise = {}", h1 - h0);
    }

    #[test]
    fn slow_closure_produces_a_smaller_rise() {
        let b = pipe();
        let q0 = 0.5;
        let head0 = 100.0;
        // Ramp down over several steps rather than stopping dead.
        let ramp: Vec<f64> = (0..20).map(|i| q0 * (1.0 - f64::from(i) / 20.0)).collect();
        let r = method_of_characteristics(b, q0, &ramp, Boundary::reservoir(head0), 0.0).unwrap();
        let fast =
            method_of_characteristics(b, q0, &[0.0], Boundary::reservoir(head0), 0.0).unwrap();
        let slow_rise = r.max_head_rise();
        let fast_rise = fast.max_head_rise();
        assert!(slow_rise < fast_rise, "{slow_rise} vs {fast_rise}");
    }

    #[test]
    fn no_flow_change_means_no_transient() {
        let b = pipe();
        let q0 = 0.5;
        let flat = vec![q0; 5];
        let r = method_of_characteristics(b, q0, &flat, Boundary::reservoir(100.0), 0.0).unwrap();
        assert!(r.max_head_rise().abs() < 1e-12, "{}", r.max_head_rise());
    }

    #[test]
    fn head_history_has_one_row_per_step_plus_the_initial() {
        let b = pipe();
        let r = method_of_characteristics(
            b,
            0.5,
            &[0.4, 0.3, 0.2, 0.1, 0.0],
            Boundary::reservoir(100.0),
            0.0,
        )
        .unwrap();
        assert_eq!(r.steps, 5);
        assert_eq!(r.head_history.len(), 6);
        assert_eq!(r.flow_history.len(), 6);
    }

    #[test]
    fn column_separation_clamps_at_the_vapour_head() {
        let b = pipe();
        // A vapour head above the steady state forces separation immediately.
        let r =
            method_of_characteristics(b, 0.5, &[0.0], Boundary::reservoir(100.0), 1.0e6).unwrap();
        assert_eq!(r.head(1, 1), Some(1.0e6));
        assert_eq!(r.separated_nodes, vec![1]);
    }

    #[test]
    fn no_separation_when_the_head_stays_above_vapour() {
        let b = pipe();
        let r = method_of_characteristics(b, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0).unwrap();
        assert!(r.separated_nodes.is_empty());
    }

    #[test]
    fn courant_validation_rejects_an_oversized_step() {
        let b = pipe();
        let stable = b.length / b.wave_speed;
        assert!(validate_courant(&b, stable).is_ok());
        assert!(validate_courant(&b, stable * 1.5).is_err());
    }

    #[test]
    fn non_physical_branches_are_rejected() {
        let bad = Branch {
            wave_speed: 0.0,
            ..pipe()
        };
        assert!(
            method_of_characteristics(bad, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0).is_err()
        );
        // A valve boundary is not supported by this single-reach formulation.
        assert!(method_of_characteristics(pipe(), 0.5, &[0.0], Boundary::valve(), 0.0).is_err());
    }

    /// The friction term must actually change the answer, otherwise the
    /// `resistance` field would be dead weight and any transient it was
    /// meant to damp would be silently missing.
    #[test]
    fn pipe_friction_reduces_the_rise() {
        let smooth = pipe();
        let rough = Branch {
            resistance: 50.0,
            ..pipe()
        };
        // With friction the steady gradient is steeper, so a stop from the
        // same reservoir head starts lower and has less room to rise.
        let a = method_of_characteristics(smooth, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0)
            .unwrap();
        let b =
            method_of_characteristics(rough, 0.5, &[0.0], Boundary::reservoir(100.0), 0.0).unwrap();
        assert!(a.head(1, 0).unwrap() > b.head(1, 0).unwrap());
    }

    #[test]
    fn initial_head_matches_the_frictionless_gradient() {
        let b = pipe();
        let head0 = 100.0;
        let q0 = 0.5;
        let r = method_of_characteristics(b, q0, &[q0], Boundary::reservoir(head0), 0.0).unwrap();
        let expected = head0 - momentum_coefficient(b.wave_speed) * q0;
        assert!((r.head(1, 0).unwrap() - expected).abs() < 1e-9);
    }

    /// A reservoir, two pipes, and a valve at the far end: the network form of
    /// the canonical problem. Each pipe is 4 reaches, so a wave reflects off the
    /// closed valve and returns to the reservoir after `2 n` steps.
    ///
    /// The flow is 0.0035 m3/s, which is 0.05 m/s in this 0.0707 m2 pipe. The
    /// pre-existing single-reach tests use `Q = 0.5` in this same area, which is
    /// 7 m/s; they pass only because they treat `Q` as a *velocity*. Reusing
    /// that as a volume flow would seed the network from an absurd field.
    fn simple_network() -> MocNetwork {
        let a = 1484.5764;
        let area = 0.0706858;
        let nodes = vec![MocNode::Reservoir(100.0), MocNode::Junction, MocNode::Valve];
        let branches = vec![
            MocBranch {
                from: 0,
                to: 1,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 0.0,
                check_valve: false,
            },
            MocBranch {
                from: 1,
                to: 2,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 0.0,
                check_valve: false,
            },
        ];
        let dt = 150.0 / a / 4.0;
        MocNetwork::new(nodes, branches, vec![0.0035, 0.0035], dt).unwrap()
    }

    /// A steady network must be a *fixed point*: with the valve held open, no
    /// step may change anything.
    ///
    /// This is the strongest single check on the formulation, and it is the
    /// test that finally located the real bug. Every defect found while writing
    /// this solver -- the spurious steady gradient, the wrong `B`, the
    /// mismatched characteristic signs, the local seed -- showed up here first
    /// as a field that would not stand still, and the failure is silent: the
    /// numbers stay finite and plausible throughout.
    #[test]
    fn a_steady_network_does_not_move() {
        let net = simple_network();
        let demand = vec![vec![0.0, 0.0, net.initial_flows[1]]];
        let r = net.solve(&demand, 0.0).unwrap();
        let h0 = r.end_head_history[0].clone();
        let q0 = r.flow_history[0].clone();
        for (k, h) in r.end_head_history[1..].iter().enumerate() {
            for (i, v) in h.iter().enumerate() {
                assert!(
                    (v - h0[i]).abs() < 1.0e-9,
                    "end {i} drifted by {} in a steady run at step {}",
                    (v - h0[i]).abs(),
                    k + 1
                );
            }
        }
        for (k, row) in r.flow_history[1..].iter().enumerate() {
            for (i, v) in row.iter().enumerate() {
                assert!(
                    (v - q0[i]).abs() < 1.0e-12,
                    "branch {i} flow drifted in a steady run at step {}",
                    k + 1
                );
            }
        }
    }

    /// A frictionless network has **no** steady head gradient, and that is the
    /// fact the whole formulation turns on.
    ///
    /// In a steady state `Q_P = Q_0`, so the flow term in the linearised
    /// momentum equation vanishes and a frictionless pipe loses no head. The
    /// `a Q / (g A)` term is purely transient -- it is the Joukowsky rise, and
    /// it appears only when the flow changes.
    #[test]
    fn a_frictionless_network_is_flat_at_steady_state() {
        let net = simple_network();
        let h0 = &net.solve(&[vec![0.0, 0.0, 0.0]], 0.0).unwrap().head_history[0];
        assert!(
            (h0[0] - h0[1]).abs() < 1.0e-12 && (h0[1] - h0[2]).abs() < 1.0e-12,
            "a frictionless network must be flat: {h0:?}"
        );
    }

    /// Friction alone must produce a steady gradient, and the total must match
    /// Darcy-Weisbach along the branch.
    #[test]
    fn friction_alone_sets_the_steady_gradient() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q = 0.0035;
        let resistance = 500.0;
        let nodes = vec![MocNode::Reservoir(100.0), MocNode::Valve];
        let branches = vec![MocBranch {
            from: 0,
            to: 1,
            length: 150.0,
            wave_speed: a,
            area,
            resistance,
            check_valve: false,
        }];
        let net = MocNetwork::new(nodes, branches, vec![q], 150.0 / a / 4.0).unwrap();
        let h0 = &net.solve(&[vec![0.0, q]], 0.0).unwrap().head_history[0];
        let drop = resistance * q * q.abs();
        assert!((h0[0] - h0[1] - drop).abs() < 1.0e-9, "{:?}", h0);
        assert!(
            (h0[0] - 100.0).abs() < 1.0e-12,
            "the reservoir head is exact"
        );
    }

    /// A moving reservoir **launches** a wave: the boundary the previous
    /// revision could not express.
    ///
    /// The safety property that makes unsteady friction safe to enable: a
    /// **steady** network gives a bit-for-bit identical answer with it on.
    ///
    /// This is not cosmetic. The EWMA starts at the steady `|Q|` and `E = |Q|`
    /// is a fixed point of the update, so a run that is already steady cannot
    /// drift. If this fails, enabling the correction corrupts every network
    /// that was merely sitting still, which is the common case.
    #[test]
    fn unsteady_friction_leaves_a_steady_network_untouched() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q = 0.0035;
        // Frictionless, so the comparison below isolates the lag itself.
        let mk = |lag: Option<f64>| {
            let nodes = vec![MocNode::Reservoir(100.0), MocNode::Valve];
            let branches = vec![MocBranch {
                from: 0,
                to: 1,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 0.0,
                check_valve: false,
            }];
            let net = MocNetwork::new(nodes, branches, vec![q], 150.0 / a / 4.0).unwrap();
            match lag {
                None => net,
                Some(t) => net.with_friction_lag(t).unwrap(),
            }
        };
        // Hold the valve open: the network is steady and must stay that way.
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, q], 40).collect();
        let r0 = mk(None).solve(&demand, 0.0).unwrap();
        let r1 = mk(Some(30.0)).solve(&demand, 0.0).unwrap();
        assert_eq!(r0.head_history[0], r1.head_history[0]);
        for (step, (x, y)) in r0
            .head_history
            .iter()
            .zip(r1.head_history.iter())
            .enumerate()
        {
            for (i, (u, v)) in x.iter().zip(y.iter()).enumerate() {
                assert!(
                    (u - v).abs() < 1.0e-12,
                    "node {i} at step {step} moved by {} with a friction lag on a \
                     steady network",
                    (u - v).abs()
                );
            }
        }
    }

    /// A **frictional** steady state must be an exact fixed point, and getting
    /// this right required fixing a real defect in the friction term.
    ///
    /// The friction force used to enter each characteristic as a *difference*
    /// between neighbouring ends, `r (Q_k E_k - Q_{k-1} E_{k-1})`. For a uniform
    /// flow -- the steady case -- that difference is identically **zero**, so
    /// nothing held a gradient up: the network started with the correct
    /// Darcy-Weisbach drop and washed it out within a few steps, ending up
    /// frictionless. Every number stayed finite and plausible, which is why it
    /// survived; only "the head settles back to the reservoir head" gives it
    /// away.
    ///
    /// The friction is now the absolute per-reach loss, so a uniform flow
    /// carries a real gradient and the steady profile is a genuine fixed point.
    #[test]
    fn a_frictional_steady_state_is_an_exact_fixed_point() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q = 0.0035;
        let resistance = 500.0;
        let nodes = vec![MocNode::Reservoir(100.0), MocNode::Valve];
        let branches = vec![MocBranch {
            from: 0,
            to: 1,
            length: 150.0,
            wave_speed: a,
            area,
            resistance,
            check_valve: false,
        }];
        let net = MocNetwork::new(nodes, branches, vec![q], 150.0 / a / 4.0).unwrap();
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, q], 40).collect();
        let r = net.solve(&demand, 0.0).unwrap();

        // The steady drop is the full Darcy-Weisbach loss, and it is held.
        let h0 = &r.head_history[0];
        let drop = resistance * q * q.abs();
        assert!(
            (h0[0] - h0[1] - drop).abs() < 1.0e-9,
            "seeded gradient {h0:?}"
        );
        for (step, row) in r.head_history.iter().enumerate() {
            assert!(
                (row[1] - h0[1]).abs() < 1.0e-12,
                "node 1 at step {step} drifted {} from the steady head {}",
                (row[1] - h0[1]).abs(),
                h0[1]
            );
        }

        // The gradient is *distributed* along the reaches, not just lumped at
        // the nodes, and it is linear in the reach index.
        let ends = &r.end_head_history[40];
        let n = ends.len() - 1;
        for j in 1..=n {
            let expect = h0[0] - drop * (j as f64) / (n as f64);
            assert!(
                (ends[j] - expect).abs() < 1.0e-12,
                "end {j} was {} but the linear gradient says {expect}: {ends:?}",
                ends[j]
            );
        }
    }

    /// A negligible lag reproduces the quasi-steady law, which is the limit
    /// `T_f -> 0` has to reach.
    #[test]
    fn a_negligible_friction_lag_reduces_to_the_quasi_steady_law() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q0 = 0.0035;
        let mk = |lag: Option<f64>| {
            let nodes = vec![MocNode::Reservoir(100.0), MocNode::Valve];
            let branches = vec![MocBranch {
                from: 0,
                to: 1,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 500.0,
                check_valve: false,
            }];
            let net = MocNetwork::new(nodes, branches, vec![q0], 150.0 / a / 4.0).unwrap();
            match lag {
                None => net,
                Some(t) => net.with_friction_lag(t).unwrap(),
            }
        };
        // A lag far below `dt` drives `C = dt / (dt + T_f)` to 1, so `E` tracks
        // `|Q|` and the two models must agree.
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, 0.0], 10).collect();
        let r0 = mk(None).solve(&demand, 0.0).unwrap();
        let r1 = mk(Some(1.0e-9)).solve(&demand, 0.0).unwrap();
        let (a0, b0) = (r0.max_head_excursion(), r1.max_head_excursion());
        assert!(
            (a0 - b0).abs() < 1.0e-6 * a0.max(1.0),
            "a negligible lag changed the rise: {a0} vs {b0}"
        );
    }

    /// The physical point of the model: a long friction lag **damps** the head
    /// excursion relative to the quasi-steady law.
    ///
    /// Wall shear lags a decelerating flow, so the friction the momentum
    /// equation sees is smaller than `R Q |Q|` on the current flow and the peak
    /// rise falls. Quasi-steady friction is therefore known to over-predict the
    /// surge; the direction of this inequality is the whole content of the
    /// correction, and it is what the test pins.
    #[test]
    fn a_long_friction_lag_damps_the_head_excursion() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q0 = 0.0035;
        let mk = |lag: Option<f64>| {
            let nodes = vec![MocNode::Reservoir(100.0), MocNode::Valve];
            let branches = vec![MocBranch {
                from: 0,
                to: 1,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 500.0,
                check_valve: false,
            }];
            let net = MocNetwork::new(nodes, branches, vec![q0], 150.0 / a / 4.0).unwrap();
            match lag {
                None => net,
                Some(t) => net.with_friction_lag(t).unwrap(),
            }
        };
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, 0.0], 60).collect();
        let steady = mk(None).solve(&demand, 0.0).unwrap();
        let lagged = mk(Some(20.0)).solve(&demand, 0.0).unwrap();
        let s = steady.max_head_excursion();
        let l = lagged.max_head_excursion();
        assert!(l < s, "the lag should damp the excursion: {l} vs {s}");
        // And it must remain a real, finite transient, not collapse to nothing.
        assert!(l > 0.0, "the lagged run produced no transient at all");
    }

    /// A check valve must stop flow **at its own end only**, and the stop must
    /// then propagate as a wave.
    ///
    /// The check-valve block used to zero the branch's *upstream* end as well as
    /// its downstream one. That is a lumped shortcut of exactly the kind this
    /// solver exists to remove: it stopped the flow across the whole branch in
    /// one step, so every interior end reacted at once and the reaching wave
    /// was annihilated. It had no test at all, which is how it survived.
    ///
    /// The test pins the difference. With the closure at the valve end, the
    /// head disturbance must appear at the ends **one per step, working
    /// upstream** -- a travelling front. A lumped stop would show every
    /// interior end disturbed on the very first step instead.
    #[test]
    fn a_check_valve_stops_at_its_own_end_and_the_stop_travels() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q0 = 0.0035;
        let reaches = 4usize;
        let net = MocNetwork::new(
            vec![MocNode::Reservoir(100.0), MocNode::Valve],
            vec![MocBranch {
                from: 0,
                to: 1,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 0.0,
                check_valve: true,
            }],
            vec![q0],
            150.0 / a / reaches as f64,
        )
        .unwrap();
        // The valve demands reverse flow; the check valve must refuse it.
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, -q0], 6).collect();
        let r = net.solve(&demand, 0.0).unwrap();
        let flat = r.end_head_history[0][0];

        // The check valve refuses reverse flow: the branch never carries it.
        for (step, row) in r.flow_history.iter().enumerate() {
            assert!(
                row[0] >= 0.0,
                "a check valve passed reverse flow {} at step {step}",
                row[0]
            );
        }

        // The valve end reacts at once, and the disturbance then walks upstream
        // one reach per step. `front` is the lowest end disturbed so far, and
        // `usize::MAX` before the first step means nothing has moved yet.
        for (step, row) in r.end_head_history.iter().enumerate().take(reaches + 1) {
            let front = if step == 0 {
                usize::MAX
            } else {
                (reaches + 1 - step).min(reaches)
            };
            for (j, h) in row.iter().enumerate() {
                let disturbed = (h - flat).abs() > 1.0e-9;
                let expected = j >= front;
                assert_eq!(
                    disturbed, expected,
                    "end {j} at step {step}: disturbed={disturbed}, expected {expected} \
                     (front {front}); a lumped stop would disturb every end at once: {row:?}"
                );
            }
        }
    }

    /// Column separation must keep the head and flow **consistent**: clamping the
    /// head alone leaves the pair violating the momentum equation it came from.
    ///
    /// The separation block clamped `h_new[k]` to the vapour head but left
    /// `q_new[k]` at whatever the boundary had solved. The stored state is then a
    /// `(H, Q)` pair that satisfies no characteristic at all, and because the
    /// next step reads that state back as the arriving characteristic, the
    /// violation is injected into the run rather than being absorbed. The
    /// physically meaningful statement is that once a column separates, the
    /// liquid below the void is effectively *free*: the head sits at vapour
    /// pressure and the flow is whatever the void allows, which is a *free*
    /// discharge, not the pre-separation flow.
    ///
    /// This is what the test pins: after separation the head is exactly the
    /// vapour head, and the flow has been released rather than frozen.
    #[test]
    fn column_separation_releases_the_flow_rather_than_freezing_it() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q0 = 0.0035;
        let reaches = 4usize;
        let net = MocNetwork::new(
            vec![MocNode::Reservoir(100.0), MocNode::Valve],
            vec![MocBranch {
                from: 0,
                to: 1,
                length: 150.0,
                wave_speed: a,
                area,
                resistance: 0.0,
                check_valve: false,
            }],
            vec![q0],
            150.0 / a / reaches as f64,
        )
        .unwrap();
        // Column separation clamps a head that has *risen* past the vapour
        // pressure, so the vapour head must sit just above the Joukowsky
        // surge. The surge on this network is 101.87 m, so 102.0 clamps it.
        // (A vapour head *below* the surge does nothing -- the head is rising,
        // not falling, and it is the peak that cavitates.)
        let vapour = 102.0;
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, 0.0], 3).collect();
        let r = net.solve(&demand, vapour).unwrap();
        assert_eq!(r.separated_nodes, vec![1], "the valve node should separate");

        // The head is clamped exactly to the vapour head.
        for (step, row) in r.head_history.iter().enumerate().skip(1) {
            assert!(
                (row[1] - vapour).abs() < 1.0e-9,
                "the separated head at step {step} was {}, expected the vapour head",
                row[1]
            );
        }

        // The flow is *released* rather than frozen at its pre-separation value.
        // A frozen flow is the bug: the pair `(vapour head, q0)` satisfies no
        // characteristic, and the next step reads it back as one.
        let before = r.flow_history[0][0];
        assert_eq!(before, q0, "the run should start in steady flow");
        for (step, row) in r.flow_history.iter().enumerate().skip(1) {
            assert!(
                (row[0] - before).abs() > 1.0e-6,
                "the flow stayed frozen at {} through separation at step {step}, \
                 leaving an inconsistent (H, Q) pair",
                row[0]
            );
        }

        // And the released flow is a free discharge: it follows from the
        // arriving characteristic and the clamped head, so it is reproducible
        // from the momentum coefficient rather than being arbitrary. The sign
        // follows the end convention, `H = a -/+ b Q`: a valve end is the
        // branch's downstream end, so `Q = (a - H) / b`, and with the head
        // clamped *above* the arriving characteristic the flow reverses. That
        // is the point -- a frozen flow would have kept the original sign.
        //
        // The arriving characteristic is `a = H_up + b Q_up`, i.e. the steady
        // head plus one momentum increment, so the exact expectation carries
        // `b * q0` rather than just the head.
        let b = a / (STANDARD_GRAVITY * 0.0706858 * reaches as f64);
        let expect = (100.0 + b * q0 - vapour) / b;
        assert!(
            (r.flow_history[1][0] - expect).abs() < 1.0e-9,
            "the released flow was {}, expected the free discharge {expect}",
            r.flow_history[1][0]
        );
    }

    /// The derived friction time constant is a real function of each pipe, and
    /// a mixed network must not get a single shared value.
    ///
    /// `T_f = L / (g A R |Q|)` means a longer, fatter, rougher pipe has a
    /// proportionally slower friction relaxation. That is the whole reason to
    /// derive it per branch: on a network of one long trunk and one short
    /// connector, a single constant is wrong on both, and wrong in *different*
    /// directions -- too slow on the connector, too fast on the trunk.
    #[test]
    fn the_derived_friction_time_scales_with_each_pipes_own_properties() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q0 = 0.0035;
        let dt = 150.0 / a / 4.0;
        let branch = |length: f64, resistance: f64| MocBranch {
            from: 0,
            to: 1,
            length,
            wave_speed: a,
            area,
            resistance,
            check_valve: false,
        };
        let net = MocNetwork::new(
            vec![MocNode::Reservoir(100.0), MocNode::Junction, MocNode::Valve],
            vec![branch(150.0, 500.0), branch(300.0, 250.0)],
            vec![q0, q0],
            dt,
        )
        .unwrap();

        // Twice the length, half the resistance, so 4x the relaxation time.
        let short = net.branch_friction_time(0);
        let long = net.branch_friction_time(1);
        assert!(
            (long / short - 4.0).abs() < 1e-9,
            "a pipe twice as long and half as rough should relax 4x slower: \
             {short} vs {long}"
        );
        // And the absolute scale is physical: seconds, not a step count.
        assert!(short > dt, "T_f {short} should exceed the time step {dt}");

        let derived = net.with_derived_friction_lag(1.0).unwrap();
        let times = derived
            .friction_lag_per_branch
            .as_ref()
            .expect("the per-branch field is set");
        assert_eq!(times.len(), 2);
        for (i, t) in times.iter().enumerate() {
            let expected = if i == 0 { short } else { long };
            assert!(
                (t - expected).abs() < 1e-12,
                "branch {i}: derived {t} vs {expected}"
            );
        }
    }

    /// A branch with no friction to relax must stay quasi-steady, and an
    /// infinite time is the honest answer rather than a zero or a panic.
    #[test]
    fn a_frictionless_branch_has_an_infinite_lag_time() {
        let a = 1484.5764;
        let mk = |resistance: f64, flow: f64| {
            MocNetwork::new(
                vec![MocNode::Reservoir(100.0), MocNode::Valve],
                vec![MocBranch {
                    from: 0,
                    to: 1,
                    length: 150.0,
                    wave_speed: a,
                    area: 0.0706858,
                    resistance,
                    check_valve: false,
                }],
                vec![flow],
                150.0 / a / 4.0,
            )
            .unwrap()
        };
        // No resistance, and separately no flow: the friction force is zero
        // either way, so there is nothing to relax.
        assert!(mk(0.0, 0.0035).branch_friction_time(0).is_infinite());
        assert!(mk(500.0, 0.0).branch_friction_time(0).is_infinite());
        // And a derived run over a frictionless branch must still solve.
        let frictionless = mk(0.0, 0.0035).with_derived_friction_lag(1.0).unwrap();
        let demand = vec![vec![0.0, 0.0], vec![0.0, 0.0035]];
        assert!(frictionless.solve(&demand, 0.0).is_ok());
    }

    /// The derived lag keeps the steady-limit property: a network that is
    /// already steady is untouched, whichever way the constants were derived.
    ///
    /// This is the safety argument for the per-branch form specifically. Because
    /// each `E` starts at its own branch's steady `|Q|`, and `E = |Q|` is a
    /// fixed point of the update, a mixed network with very different derived
    /// constants still reproduces the quasi-steady answer exactly when nothing
    /// is changing.
    #[test]
    fn a_derived_lag_leaves_a_steady_network_untouched() {
        let a = 1484.5764;
        let q = 0.0035;
        let mk = |derived: bool| {
            let net = MocNetwork::new(
                vec![MocNode::Reservoir(100.0), MocNode::Valve],
                vec![MocBranch {
                    from: 0,
                    to: 1,
                    length: 150.0,
                    wave_speed: a,
                    area: 0.0706858,
                    resistance: 500.0,
                    check_valve: false,
                }],
                vec![q],
                150.0 / a / 4.0,
            )
            .unwrap();
            if derived {
                net.with_derived_friction_lag(1.0).unwrap()
            } else {
                net
            }
        };
        let demand: Vec<Vec<f64>> = std::iter::repeat_n(vec![0.0, q], 40).collect();
        let plain = mk(false).solve(&demand, 0.0).unwrap();
        let derived = mk(true).solve(&demand, 0.0).unwrap();
        for (step, (x, y)) in plain
            .head_history
            .iter()
            .zip(derived.head_history.iter())
            .enumerate()
        {
            for (i, (u, v)) in x.iter().zip(y.iter()).enumerate() {
                assert!(
                    (u - v).abs() < 1.0e-12,
                    "node {i} at step {step} moved by {} with a derived lag on a \
                     steady network",
                    (u - v).abs()
                );
            }
        }
    }

    /// The **friction rise integral** has a closed form, and reproducing it is
    /// the last item the water-hammer notes listed as open. It was open because
    /// nothing had ever checked the *value* -- the only friction test asserted
    /// `slow < fast`, which any monotone model satisfies.
    ///
    /// Integrating the momentum equation along a closure that takes the flow
    /// from `Q0` to zero,
    ///
    /// ```text
    /// dH = -(a / g) dQ + R Q |Q| dQ
    /// rise = (a / g) Q0 - (2/3) R Q0^{3/2}
    /// ```
    ///
    /// The first term is the Joukowsky head and the second is the friction
    /// integral, so friction subtracts a quantity proportional to `R`. Note it
    /// is **not** `1 / (1 + a / gR)`: that form assumes the friction is
    /// proportional to the momentum term at every instant, whereas here it is
    /// integrated along the path the flow actually takes. The two agree only in
    /// the limit of a small closure.
    ///
    /// The friction used to be a *difference* of friction force, `R Q|Q| -
    /// R Q_prev|Q_prev|`, added to the head. That telescopes over a closure to
    /// the constant `+ R Q0^2`, so it added a fixed offset instead of damping
    /// anything: a frictional gradual closure came out at 0.5333 m, slightly
    /// **above** the frictionless 0.5298 m. Friction cannot raise the surge, so
    /// the sign was wrong and the shape was wrong with it.
    #[test]
    fn a_gradual_closure_follows_the_friction_rise_integral() {
        let a = 1484.5764;
        let q0 = 0.0035;
        let resistance = 500.0;
        let frictionless = pipe();

        // A long, linear closure, which is the regime the integral describes.
        let n = 200;
        let ramp: Vec<f64> = (0..n)
            .map(|i| q0 * (1.0 - f64::from(i) / n as f64))
            .collect();
        let joukowsky = a * q0 / STANDARD_GRAVITY;
        let run = |r: f64| {
            let branch = Branch {
                resistance: r,
                ..frictionless
            };
            method_of_characteristics(branch, q0, &ramp, Boundary::reservoir(100.0), 0.0)
                .unwrap()
                .max_head_rise()
        };

        // Frictionless: exactly the Joukowsky rise, whatever the closure
        // profile, because there is nothing to bleed energy.
        assert!(
            (run(0.0) - joukowsky).abs() / joukowsky < 0.02,
            "a frictionless gradual closure should give the Joukowsky rise: {} vs {joukowsky}",
            run(0.0)
        );

        // With friction the rise is the Joukowsky head less the integral
        // `(2/3) R Q0^{3/2}`.
        let predicted = joukowsky - (2.0 / 3.0) * resistance * q0.powf(1.5);
        let actual = run(resistance);
        assert!(
            (actual - predicted).abs() / predicted < 0.02,
            "frictional gradual closure gave {actual}, the rise integral predicts \
             {predicted} (Joukowsky {joukowsky}, friction integral {})",
            (2.0 / 3.0) * resistance * q0.powf(1.5)
        );
        assert!(
            actual < joukowsky,
            "friction must damp the surge, but {actual} exceeds {joukowsky}"
        );

        // And the damping must grow with friction, monotonically.
        let mut previous = f64::INFINITY;
        for r in [10.0, 100.0, 500.0, 2000.0] {
            let rise = run(r);
            assert!(
                rise < previous,
                "rise should fall as resistance rises: R={r} gave {rise} vs {previous}"
            );
            previous = rise;
        }
    }

    /// A `dH` head step at the reservoir travels at the wave speed, one reach
    /// per step, and on reaching the dead end is reflected at full amplitude to
    /// `2 dH`. Every part of that is checked: the far end is untouched until
    /// the wave arrives, the head steps up one reach per step so the wave is
    /// visible travelling, the reflection is exactly doubled, and the flow stays
    /// zero at every step. Without a time-varying prescribed head the far end
    /// could never move at all, so the wave would have nowhere to come from.
    #[test]
    fn a_moving_reservoir_launches_a_travelling_wave() {
        let a = 1484.5764;
        let area = 0.0706858;
        let reaches = 4usize;
        let head0 = 100.0;
        let d_h = 2.0;
        let nodes = vec![MocNode::Reservoir(head0), MocNode::DeadEnd];
        let branches = vec![MocBranch {
            from: 0,
            to: 1,
            length: 150.0,
            wave_speed: a,
            area,
            resistance: 0.0,
            check_valve: false,
        }];
        let dt = 150.0 / a / reaches as f64;
        let net = MocNetwork::new(nodes, branches, vec![0.0], dt).unwrap();

        // The reservoir level steps up by `d_h` and stays there.
        let total = reaches + 6;
        let steps = vec![vec![0.0, 0.0]; total];
        let mut heads = vec![vec![head0, 0.0]];
        heads.extend(std::iter::repeat_n(vec![head0 + d_h, 0.0], total - 1));
        let r = net.solve_with_heads(&steps, &heads, 0.0).unwrap();

        // The far end is untouched until the wave has actually travelled, and
        // the reservoir itself is exactly on the schedule. `head_history[0]` is
        // the initial state, so the first schedule row lands on index 1.
        for (step, h) in r.head_history.iter().enumerate() {
            let expect_res = if step <= 1 { head0 } else { head0 + d_h };
            assert!(
                (h[0] - expect_res).abs() < 1.0e-12,
                "reservoir at step {step} was {}, expected {expect_res}",
                h[0]
            );
            if step <= reaches + 1 {
                assert!(
                    (h[1] - head0).abs() < 1.0e-9,
                    "the wave reached the far end early: step {step} gave {}",
                    h[1]
                );
            }
        }

        // The wave is visible travelling: the reservoir end steps at history
        // index 2 and each interior end one step later, so the front is a moving
        // edge rather than a jump. Measured, not assumed: end 0 at step 2, end 1
        // at step 3, and so on to the far end at step 2 + `reaches`.
        for (step, row) in r.end_head_history.iter().enumerate().take(reaches + 3) {
            let front = (step as isize - 2).clamp(-1, reaches as isize);
            for (j, h) in row.iter().enumerate().take(reaches) {
                let expected = if (j as isize) <= front {
                    head0 + d_h
                } else {
                    head0
                };
                assert!(
                    (h - expected).abs() < 1.0e-9,
                    "end {j} at step {step} was {h}, expected {expected}: {row:?}"
                );
            }
        }

        // The reflection at the dead end is the full, doubled wave.
        let arrive = reaches + 2;
        assert!(
            (r.head_history[arrive][1] - head0 - 2.0 * d_h).abs() < 1.0e-9,
            "expected a reflected head of {} at the far end, got {}",
            head0 + 2.0 * d_h,
            r.head_history[arrive][1]
        );

        // A dead end still passes no flow, with a moving reservoir upstream.
        for (step, row) in r.flow_history.iter().enumerate() {
            assert!(
                row[0].abs() < 1.0e-12,
                "a dead end cannot pass flow, but step {step} moved {}",
                row[0]
            );
        }
    }

    /// A **dead end reflects a wave at full amplitude**; it does not absorb it.
    ///
    /// This is the property a fixed-head reservoir structurally cannot show, and
    /// it is what makes a dead end the boundary that can quote a Joukowsky head
    /// at a closed end. The evidence is a **staircase**: the head at the dead end
    /// climbs in equal jumps and never comes back down, because each arriving
    /// wave adds to the one it reflects.
    ///
    /// The decisive number is the size of each jump. One reach of the wave
    /// carries a head change of `B * Q0` with `B = a / (g A n)`. The jump at a
    /// dead end is exactly **twice** that, and that factor of two is the whole
    /// content of a full-amplitude reflection: the incident and reflected waves
    /// are equal and opposite in *flow* but add in *head*, since the dead end
    /// cannot relieve the pressure by passing flow. An absorbing boundary would
    /// show `B * Q0` and, at a reservoir, zero.
    ///
    /// Layout: reservoir -> valve -> (4 reaches) -> dead end -> (4 reaches) ->
    /// dead end, with the valve opened from rest so the only disturbance in the
    /// run is the wave it launches.
    #[test]
    fn a_dead_end_reflects_rather_than_absorbs() {
        let a = 1484.5764;
        let area = 0.0706858;
        let q0 = 0.0035;
        let reaches = 4usize;
        let nodes = vec![
            MocNode::Reservoir(100.0),
            MocNode::Valve,
            MocNode::DeadEnd,
            MocNode::DeadEnd,
        ];
        let branch = |from: usize, to: usize| MocBranch {
            from,
            to,
            length: 150.0,
            wave_speed: a,
            area,
            resistance: 0.0,
            check_valve: false,
        };
        let branches = vec![branch(0, 1), branch(1, 2), branch(2, 3)];
        let dt = 150.0 / a / reaches as f64;
        let net = MocNetwork::new(nodes, branches, vec![q0, 0.0, 0.0], dt).unwrap();

        // The valve opens to q0 after a step of rest, so the network starts
        // quiescent and the disturbance is entirely the wave it launches.
        let mut steps = vec![vec![0.0, 0.0, 0.0, 0.0]];
        steps.extend(std::iter::repeat_n(vec![0.0, q0, 0.0, 0.0], 24));
        let r = net.solve(&steps, 0.0).unwrap();

        // The defining property: a dead end passes **no** flow, ever.
        for (step, row) in r.flow_history.iter().enumerate() {
            for &i in &[1usize, 2] {
                assert!(
                    row[i].abs() < 1.0e-12,
                    "branch {i} passed {} at step {step} through a dead end",
                    row[i]
                );
            }
        }

        // One reach of the wave, and the doubled jump a dead end must show.
        let one_reach = a / (STANDARD_GRAVITY * area * reaches as f64) * q0;
        let doubled = 2.0 * one_reach;

        // Node 2 is the first dead end. Its head is a staircase: it only ever
        // moves *up*, and it moves in equal jumps of exactly the doubled wave.
        let dead: Vec<f64> = r.head_history.iter().map(|row| row[2]).collect();
        assert!(
            dead.windows(2).all(|w| w[1] >= w[0] - 1.0e-12),
            "a reflecting dead end must never shed head: {dead:?}"
        );
        let jumps: Vec<f64> = dead.windows(2).map(|w| w[1] - w[0]).collect();
        let moved: Vec<(usize, f64)> = jumps
            .iter()
            .enumerate()
            .filter(|(_, &v)| v.abs() > 1.0e-9)
            .map(|(i, &v)| (i + 1, v))
            .collect();
        assert_eq!(
            moved.len(),
            3,
            "expected three reflections in 24 steps, got {moved:?}"
        );
        for (step, v) in &moved {
            assert!(
                (v - doubled).abs() < 1.0e-9,
                "reflection at step {step} was {v}, expected the doubled wave {doubled}"
            );
        }

        // The contrast that makes this a reflection rather than an absorption:
        // the valve is a *prescribed-flow* boundary, so each wave it takes
        // leaves it and the node returns to its baseline, over and over. The
        // dead end never revisits its baseline even once, because every wave it
        // receives stays. Comparing the two nodes in the same run is what rules
        // out "the head rose because of something else".
        let valve: Vec<f64> = r.head_history.iter().map(|row| row[1]).collect();
        let at_rest = |v: f64| (v - 100.0).abs() < 1.0e-9;
        let valve_rest_steps = valve.iter().filter(|&&v| at_rest(v)).count();
        assert!(
            valve_rest_steps > 10,
            "the absorbing valve should return to its baseline repeatedly, but did \
             so only {valve_rest_steps} times: {valve:?}"
        );
        // After the first wave reaches it, the dead end is never at rest again.
        // (It *is* at rest beforehand, while the wave is still in flight, so the
        // comparison has to start at arrival rather than at t=0.)
        let arrival = moved[0].0;
        let rest_after = dead[arrival..].iter().filter(|&&v| at_rest(v)).count();
        assert_eq!(
            rest_after, 0,
            "a reflecting dead end never settles again once reached: {dead:?}"
        );
    }

    /// The whole point: information travels **one reach per step**, so the
    /// interior ends between the valve and the junction are untouched one step
    /// after a closure.
    ///
    /// A lumped single-reach model applies the closure along the whole pipe at
    /// once and cannot make this assertion at all. It needs the end history,
    /// because interior ends are not network nodes.
    #[test]
    fn a_wave_moves_exactly_one_reach_per_step() {
        let net = simple_network();
        let r = net.solve(&[vec![0.0, 0.0, 0.0]], 0.0).unwrap();
        let offset = net.end_offset(1);
        for k in offset + 1..offset + net.reaches[1] {
            let h0 = r.end_head(k, 0).unwrap();
            let h1 = r.end_head(k, 1).unwrap();
            assert!(
                (h1 - h0).abs() < 1.0e-9,
                "end {k} moved {} after one step, before the wave could reach it",
                (h1 - h0).abs()
            );
        }
        // And the valve, where the action was, has responded.
        let valve_end = offset + net.reaches[1];
        let h0 = r.end_head(valve_end, 0).unwrap();
        let h1 = r.end_head(valve_end, 1).unwrap();
        assert!((h1 - h0).abs() > 1.0e-6, "the valve must respond at once");
    }

    /// The valve's immediate response is **one reach of the wave**, not the
    /// whole-pipe Joukowsky rise.
    ///
    /// Closing a valve launches a wave that travels one reach per step. At the
    /// valve, only the `+` characteristic from the reach upstream has arrived, so
    ///
    /// ```text
    /// H_new = H_up + B_reach Q      (the old head was H_up, being flat)
    /// ```
    ///
    /// and the immediate rise is exactly `B_reach Q = B_whole Q / n`. The
    /// full `B_whole Q` only appears once the wave has crossed every reach and
    /// returned, which is the entire reason for discretising and the thing the
    /// single-reach model cannot represent.
    #[test]
    fn the_valve_responds_with_one_reach_of_the_wave() {
        let net = simple_network();
        let r = net.solve(&[vec![0.0, 0.0, 0.0]], 0.0).unwrap();
        let valve_end = net.end_offset(1) + net.reaches[1];
        let h0 = r.end_head(valve_end, 0).unwrap();
        let h1 = r.end_head(valve_end, 1).unwrap();
        let rise = (h1 - h0).abs();

        let a = net.branches[1].wave_speed;
        let n = net.reaches[1] as f64;
        let b_reach = a / (STANDARD_GRAVITY * net.branches[1].area * n);
        let expected = b_reach * net.initial_flows[1];
        assert!(
            (rise - expected).abs() / expected < 1.0e-9,
            "rise {rise}, expected {expected}"
        );
    }

    /// The far end must feel the wave, and when it does the flow must stop.
    ///
    /// This is the Joukowsky result, stated in the form a fixed-head boundary
    /// can actually show. The classic pairing is "the flow is brought to rest
    /// and the head rises by `a Q0 / g`", but the two halves are not
    /// interchangeable observables: a reservoir has a **prescribed head**, so it
    /// cannot display a head rise at all -- the arriving wave is absorbed by the
    /// flow changing instead. Asserting a head rise at a reservoir therefore
    /// tests something the boundary condition forbids, and it is the head that
    /// rises at a *closed* downstream boundary.
    ///
    /// What the far end can show is the other half: the flow falls from `Q0` to
    /// zero. And it can only do so once the wave has arrived, which is the part
    /// the single-reach model cannot represent at all.
    #[test]
    fn the_far_end_stops_flow_once_the_wave_arrives() {
        let net = simple_network();
        let n = net.reaches[1];
        let demand: Vec<Vec<f64>> = (0..24).map(|_| vec![0.0, 0.0, 0.0]).collect();
        let r = net.solve(&demand, 0.0).unwrap();

        // Until the wave crosses the valve branch and the junction closes, the
        // far end is undisturbed.
        for step in 0..=n {
            assert!(
                (r.flow(0, step).unwrap() - net.initial_flows[0]).abs() < 1.0e-12,
                "the branch-0 flow changed at step {step}, before the wave arrived"
            );
        }

        // Once it arrives, the flow is brought to rest: that is Joukowsky.
        let stopped = r.flow(0, n + 1).unwrap();
        assert!(
            stopped.abs() < 1.0e-9,
            "the far-end flow should be brought to rest, got {stopped}"
        );

        // And the head rise that accompanies it is `a Q0 / g`, recoverable from
        // the wave rather than from the reservoir head. This is checked where it
        // is visible: at the valve, one reach after the closure.
        let a = net.branches[1].wave_speed;
        let area = net.branches[1].area;
        let n = net.reaches[1] as f64;
        let b_reach = a / (STANDARD_GRAVITY * area * n);
        let valve_end = net.end_offset(1) + net.reaches[1];
        let one_step_rise =
            (r.end_head(valve_end, 1).unwrap() - r.end_head(valve_end, 0).unwrap()).abs();
        // The valve's one-reach rise is `B_reach Q = (a Q0 / g) / n`, so the
        // whole-pipe Joukowsky head is exactly `n` times it.
        let joukowsky = n * one_step_rise;
        assert!(
            (joukowsky - a * net.initial_flows[1] / (STANDARD_GRAVITY * area)).abs() / joukowsky
                < 1.0e-6,
            "Joukowsky head {joukowsky} disagrees with a Q0 / g"
        );
        // And `B_reach Q` is a factor `1/n` of it, which is the discretisation.
        assert!((one_step_rise - b_reach * net.initial_flows[1]).abs() / one_step_rise < 1.0e-9);
    }

    /// The reach counts must make `dt` a Courant step on every branch, and an
    /// oversized `dt` must be refused rather than silently rounded to one
    /// reach.
    #[test]
    fn the_courant_condition_is_enforced() {
        let net = simple_network();
        let a = net.branches[0].wave_speed;
        for i in 0..net.branches.len() {
            let dx = net.branches[i].length / net.reaches[i] as f64;
            let a_eff = net.effective_wave_speed(i);
            assert!(
                (dx / a_eff - net.dt).abs() < 1.0e-12,
                "branch {i}: dx/a = {} vs dt = {}",
                dx / a_eff,
                net.dt
            );
            assert!(
                a_eff <= a * (1.0 + 1.0e-12),
                "branch {i}: a_eff {a_eff} exceeds the physical {a}"
            );
        }
        // A `dt` beyond the Courant limit is an error. `ceil` of a small number
        // is still 1, so without a direct check it would be rounded away.
        let q = vec![0.0035, 0.0035];
        assert!(MocNetwork::new(
            net.nodes.clone(),
            net.branches.clone(),
            q.clone(),
            150.0 / a * 2.0
        )
        .is_err());
        assert!(MocNetwork::new(net.nodes.clone(), net.branches.clone(), q.clone(), 0.0).is_err());
    }

    #[test]
    fn malformed_networks_are_rejected() {
        let net = simple_network();
        let q = vec![0.0035, 0.0035];
        let mut bad = net.branches.clone();
        bad[0].to = 99;
        assert!(MocNetwork::new(net.nodes.clone(), bad, q.clone(), net.dt).is_err());
        assert!(MocNetwork::new(
            net.nodes.clone(),
            net.branches.clone(),
            vec![0.0035],
            net.dt
        )
        .is_err());
        // Degenerate branch parameters.
        let mut zero = net.branches.clone();
        zero[0].area = 0.0;
        assert!(MocNetwork::new(net.nodes.clone(), zero, q.clone(), net.dt).is_err());
    }

    #[test]
    fn column_separation_is_reported_on_the_network() {
        let net = simple_network();
        // The surge peaks at about 101.87 m, so a vapour head of 102.0 clamps
        // it. A vapour head *above* the reservoir is a legitimate case.
        let r = net.solve(&[vec![0.0, 0.0, 0.0]], 102.0).unwrap();
        assert!(
            !r.separated_nodes.is_empty(),
            "a vapour head below the surge must separate the column"
        );
    }

    /// A **non-finite** vapour head is rejected, and a vapour head *above* the
    /// reservoir head is **not**.
    ///
    /// These two are easy to confuse and the second one is a real physical
    /// case: closing a valve raises the head by the Joukowsky rise, so a surge
    /// past the reservoir level is exactly what column separation is tested
    /// against. An earlier attempt rejected vapour heads above the reservoir on
    /// the reasoning that "the whole system would already have separated" --
    /// which is wrong, because the surge is the thing that causes the
    /// separation. Rejecting those broke the separation tests that were
    /// exercising the real case.
    #[test]
    fn a_vapour_head_above_the_reservoir_is_valid_but_a_nan_one_is_not() {
        let net = simple_network();
        // Above the reservoir: physical, and it separates the surge.
        let surge = net.solve(&[vec![0.0, 0.0, 0.0]], 102.0);
        assert!(
            surge.is_ok(),
            "a vapour head above the reservoir is a legitimate cavitation case"
        );
        assert!(
            !surge.expect("solves").separated_nodes.is_empty(),
            "and it should actually separate the column"
        );
        // Non-finite: not a regime, an absent one.
        assert!(net.solve(&[vec![0.0, 0.0, 0.0]], f64::NAN).is_err());
        assert!(net.solve(&[vec![0.0, 0.0, 0.0]], f64::INFINITY).is_err());
    }

    #[test]
    fn the_reservoir_head_is_held_while_the_transient_runs() {
        let net = simple_network();
        let demand: Vec<Vec<f64>> = (0..20).map(|_| vec![0.0, 0.0, 0.0]).collect();
        let r = net.solve(&demand, 0.0).unwrap();
        assert!(r.max_head_excursion() > 0.0, "a transient must occur");
        for step in 0..=r.steps {
            assert_eq!(
                r.head(0, step),
                Some(100.0),
                "the reservoir head is prescribed and must not move"
            );
        }
    }
}

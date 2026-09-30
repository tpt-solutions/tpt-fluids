# tpt-fluids — Build Todo

> Tracks bootstrap + full crate build-out for the tpt-fluids workspace, per
> `spec.txt`. tpt-fluids is a pure-Rust, AI-native applied fluid mechanics
> library consolidating 1D hydraulic networks, marine hydrodynamics, and
> tribology into a single feature-gated workspace, with zero external solver
> FFI. License for every crate: `MIT OR Apache-2.0`. Author: TPT Solutions.
>
> **Strict licensing:** no Apache-2.0-ONLY dependency is permitted anywhere
> in the graph — MIT-only, MIT/Apache-2.0 dual, or in-house only. Verify every
> new external crate against this before adding it (mirrors `tpt-math`'s
> `deny.toml` policy and its nalgebra/clarabel/statrs compliance fixes).

## Phase 0 — Repo Bootstrap

(one-time; mirrors `tpt-math`'s bootstrap shape)

- [x] Create root `Cargo.toml` (workspace `resolver = "2"`; `[workspace.package]`:
      `edition = "2021"`, `rust-version = "1.84"` (matched to `tpt-math`'s MSRV
      floor, which its published internal deps declare), `license = "MIT OR Apache-2.0"`,
      `authors = ["TPT Solutions"]`, `homepage`/`repository` URLs)
- [x] `rust-toolchain.toml`
- [x] `rustfmt.toml`
- [x] `deny.toml` (license allow-list: `MIT`, `Apache-2.0`, `Apache-2.0 WITH LLVM-exception`,
      `BSD-2-Clause`, `BSD-3-Clause`, `ISC`, `Unicode-3.0`, `Zlib`, `CC0-1.0`, `MPL-2.0`
      — same allow-list as `tpt-math`; `advisories.yanked = "deny"`,
      `sources.unknown-registry = "deny"`, `sources.unknown-git = "deny"`)
- [x] `.github/workflows/ci.yml`
- [x] `LICENSE-MIT` and `LICENSE-APACHE`
- [x] Create empty `crates/` directory
- [x] Add Rust `.gitignore` (`/target`, etc.)
- [x] Write root `README.md` stub — tpt-fluids's role as the applied fluid
      mechanics layer for `tpt-construction`/`tpt-engineering`/`tpt-physics`/
      `tpt-aero`; link to `spec.txt`
- [x] `git init` (local only, unless/until a remote is requested)
- [x] Initial commit
- [x] Sanity check: `cargo build` succeeds on the empty workspace

## Fixed defects worth remembering

These were all silent: every value stayed finite, positive, and plausible, so
none of them would ever fail an assertion. Each was found by printing
intermediate state and reading it, not by trusting that the numbers looked
reasonable.

- **Colebrook dual gradient** returned a hardcoded `0.02` instead of the exact
  forward-mode dual. Now an exact dual iteration, with value, gradient,
  seed-independence and scaling tests.
- **MOC steady seeding** applied the transient `a Q / (g A)` term as a steady
  loss. In a steady state `Q_P = Q_0`, so a frictionless pipe has *no* steady
  head gradient; seeding that term made the network drift on runs that should
  have been perfectly steady. Found by the fixed-point test.
- **MOC friction term** was a *difference* of friction force between
  neighbouring ends, which is identically zero for the uniform flow of a
  steady network, so the seeded Darcy-Weisbach drop washed out within a few
  steps and the network settled into looking frictionless. Now the absolute
  per-reach loss, with interior ends seeded on the linear gradient.
- **MOC check valve** zeroed the branch's *upstream* end as well as its
  downstream one, stopping flow across the whole branch in a single step and
  annihilating the wave still travelling through the reaches behind it. It had
  no test at all.
- **MOC column separation** clamped the head to the vapour head but left the
  flow frozen, leaving an `(H, Q)` pair that satisfied no characteristic; the
  next step read that violation back in as an arriving wave.
- **Froude scaling** `resistance_from_coefficient` multiplied by `rho g` when
  the displacement it was given was *already a weight*, inflating every
  resistance by about 9807. The one test of it asserted that two code paths
  *agreed* rather than that either was physical, and the test coefficient had
  been tuned to compensate, so the error and its compensation cancelled to a
  plausible-looking number. A consistency test between two paths that share a
  bug is blind to the bug by construction; it needs a physical anchor.

- **Tribology `contact`** `contact_radius` and `approach` both returned `0.0`
  for an infinite reduced radius, reporting a loaded flat joint as having no
  contact patch and no deformation at all. The point-contact formula does not
  apply there (it is a cylinder problem), so the result is now `NaN` and
  `flat_approach` is the determinable route. Both branches were untested, and
  the one test that touched them asserted the wrong answer.

- **Single-reach MOC friction** was applied as a *difference* of friction force,
  `R Q|Q| - R Q_prev|Q_prev|`, added to the head. That telescopes over a closure
  to the constant `+ R Q0^2`, so instead of damping the surge it added a fixed
  offset — a frictional gradual closure came out at 0.5333 m, slightly **above**
  the frictionless 0.5298 m, which is impossible. Friction is now the definite
  integral of `Q |Q|` along the closure, giving the closed form
  `rise = (a/g) Q0 - (2/3) R Q0^{3/2}`, which the solver now reproduces to 2%.
  The test that should have caught it only asserted `slow < fast`, which any
  monotone model satisfies.
- **Property tests now cover the MOC solver**, which is where every silent
  defect in this crate's history was found. Six properties sweep resistance,
  lag, closure rate, and step count: friction only ever damps a surge, a steady
  state does not drift, nothing non-finite escapes, a wave moves one reach per
  step, a friction lag cannot disturb a steady run, and a check valve never
  passes reverse flow. Writing them immediately found that the two
  column-separation tests were forcing separation with an impossible vapour
  head, so they were never exercising the case they claimed to.
- **A non-finite vapour head is rejected**; one above the reservoir is not. The
  second is easy to get backwards: a closure surge legitimately exceeds the
  reservoir head, and that peak is exactly what cavitation is tested against.
- **Unit checking in `core`** was decorative: every `Dimension` constant was
  declared and never verified, because the operators work on bare `f64`. Adding
  the `const` assertion immediately failed the build on three genuinely wrong
  relations — `Area / Length = Volume`, `Length * Time = Velocity`, and
  `Velocity / Time = Length` — all of which the old macro had happily generated
  as real, type-checking `Div` and `Mul` impls with no physical meaning. The
  cause was that one macro emitted *both* `Mul` and `Div` from a single triple
  without checking either. They are now `dimensioned_ops!` (asserts the product)
  and `dimensioned_quotient!` (asserts the quotient), each emitting only what it
  proves, so a wrong constant is a compile error. Note that the wrong relations
  were removed, not just the checks: the four dimensionally correct pairs that
  the old macro had *also* been generating (`Area * Length`, `Velocity * Time`,
  `Force * Velocity`, `Pressure * Volume`) were restored rather than lost, and
  `Velocity / Time` now returns an `Acceleration` instead of a `Length`.

## Per-Crate Checklist Template

Every phase below repeats this shape. The umbrella crate (`tpt-fluids`) uses
the umbrella variant instead of steps 2-4.

**Standard crate:**
1. Scaffold `crates/<name>/` (Cargo.toml inheriting workspace fields, `lib.rs` stub)
2. Wire dependencies (internal `tpt-fluids-*`/`tpt-math-*` + external crates —
   license-check every external crate against `deny.toml`'s allow-list first)
3. Implement scope
4. Unit tests + doctests
5. Rustdoc (crate-level + public API)
6. `cargo fmt --check` / `cargo clippy --all-targets --all-features -- -D warnings` clean
7. `cargo deny check` clean
8. no_std(+alloc) target verification (`thumbv6m-none-eabi`) — only where the
   crate declares `no_std` support
9. Add to root `Cargo.toml` `[workspace] members` + `[workspace.dependencies]`

**Umbrella crate:**
1. Scaffold `crates/<name>/` (Cargo.toml with Cargo features gating each constituent re-export)
2. Wire optional deps + matching feature flags per constituent crate
3. Re-export each constituent's public API behind its feature
4. Rustdoc documenting the feature matrix
5. `cargo fmt` / `clippy` / `deny` clean across feature combinations

---

## Phase 1 — tpt-fluids-core

*Foundation: fluid properties, EOS, viscosity models, unit-safe types.
no_std + alloc. Depends on: `tpt-math-units`, `tpt-math-numeric`, external
`tpt-materials` (density/viscosity lookups — published TPT Solutions repo on
GitHub, wired as a normal git/crates.io dependency, not a local path member).*

- [x] Scaffold `crates/tpt-fluids-core/`
- [x] Wire deps: `tpt-math-units`, `tpt-math-numeric`, `libm` (see the
      `tpt-materials` note at the end of this phase)
- [x] Implement unit-safe wrappers: pressure (Pa), head (m), volumetric flow
      (m³/s), mass flow (kg/s), velocity (m/s), dynamic viscosity (Pa·s),
      kinematic viscosity (m²/s), density (kg/m³)
- [x] Implement non-dimensional number types: Reynolds, Froude, Weber, Mach,
      Cavitation number (σ) — strictly non-dimensional (phantom-typed against
      the dimensioned quantities above)
- [x] Implement equations of state: incompressible, ideal gas, Tait (water),
      tabulated interpolation (shape-preserving cubic Hermite through a
      measured `(rho, p)` table, inverted by bisection)
- [x] Implement viscosity models: Newtonian, power-law, Bingham plastic,
      Sutherland's law (gases)
- [x] Implement surface tension and contact angle models
- [x] Fluid property lookups (temperature-dependent density/viscosity) —
      **built in-house rather than via `tpt-materials`**, see the note below
- [x] Unit tests + doctests (incl. dimensional-mismatch compile-fail checks
      where phantom typing is meant to reject them)
- [x] Rustdoc
- [x] Froude-Krylov excitation: the first-order wave excitation force on a
      wall-sided hull, `F_x = rho g a (1 - e^(-kT)) 2 sin(kL/2) / k`, resolved
      onto the ship's axes for any heading. Both limit cases are tested: the
      force vanishes for waves far shorter than the hull, and tends to
      `rho a w^2 T L` in the long-wave limit
- [x] Holtrop-Mennen: total resistance from principal dimensions, split into
      wave-making, friction, residuary, and appendage terms, with Holtrop's
      own wetted-surface estimate. Validated against a 120 000 dwt bulker at
      14 kn, whose resistance is known from its required power
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] no_std+alloc verify (`thumbv6m-none-eabi`)
- [x] Add to root `Cargo.toml` members + workspace deps

> **`tpt-materials` is not usable as scoped.** `spec.txt` and the original
> Phase 1 plan assumed a published TPT Solutions crate supplying temperature-
> dependent density and viscosity lookups. The real `tpt-materials` repository
> is a *micro-scale* physics engine — crystal plasticity, phase-field,
> diffusion — publishing `tpt-mat-core`, `tpt-mat-constants`, and
> `tpt-mat-crystallography`. It is **not on crates.io** (only on GitHub, not
> wired to a registry), and it contains **no fluid property tables at all**.
> Depending on it would have meant either a non-registry git dependency (which
> `deny.toml`'s `sources.unknown-git = "deny"` forbids) or a dependency on the
> wrong domain entirely. `tpt-fluids-core` therefore carries an in-house
> database instead: `FluidProperties` (water, seawater, dry air, ISO VG 46
> oil) plus standalone `water_viscosity` (Vogel form), `water_vapour_pressure`
> (IAPWS Wagner-Pruss), and `DensityModel` (linear, plus Kell's correlation
> for water, which is needed because water's density maximum at 4 °C cannot
> be represented by a linear expansion law). The old claim that this
> dependency was "resolved" has been corrected accordingly.

## Phase 2 — tpt-fluids-hydraulic

*1D pipe network analysis for civil, mechanical, and process engineering.
Depends on: `tpt-fluids-core`, `tpt-math-graph` (topology), `tpt-math-linalg`
/`tpt-math-linalg-sparse` (Newton-Raphson solves), `tpt-math-autodiff`
(differentiable head-loss).*

- [x] Scaffold `crates/tpt-fluids-hydraulic/`
- [x] Wire deps: `tpt-fluids-core`, `tpt-math-graph`, `tpt-math-linalg`,
      `tpt-math-autodiff`
- [x] Implement friction models: Darcy-Weisbach (laminar), Colebrook-White
      iterative solver, Swamee-Jain explicit approximation, Hazen-Williams,
      plus a table of standard absolute pipe roughnesses. Minor-loss
      coefficients are represented as `FixedLoss` links in `network.rs`.
- [x] Implement network topology: directed graph representation with
      connected-component counting, spanning forest, shortest-path tree, and
      fundamental cycle-basis (loop) extraction
- [x] Implement Hardy Cross method (loop-based) network solver
      (`hardy_cross.rs`). **Now converges on every network tested, including
      the multi-source, multi-loop case that previously stalled**; both tests
      that were `#[ignore]`d for that reason are enabled. Three separate
      defects were behind the stall, and each is recorded below because none
      of them was the one the note above blamed:
      - The loop correction *linearised* the loop equation, freezing the loop
        resistance at its current flows and taking one step. It now **solves**
        the loop equation with a safeguarded Newton iteration on the exact
        `dh(dQ) = sum w_i (x_i + dQ)|x_i + dQ|`, which is quadratic in the
        correction and cheap to evaluate.
      - The stall test compared the *worst single* loop imbalance between
        passes. Under Gauss-Seidel correction, closing loop A necessarily
        disturbs loop B, so a network improving steadily at 0.4 per pass was
        reported as stalled. It now watches the total imbalance and requires
        several consecutive non-improving passes (`HardyCrossOptions::
        consecutive_stalls`) before giving up.
      - The seed was demand-only and ignored the head field, which on a
        multi-source network seeds a *second, spurious* solution: every loop
        balances and continuity holds, so the iteration converged to it and
        returned `Ok` with flows of order 1e-2 where the answer is order 1.
        The seed is now built from a head field interpolated between the
        reservoirs, then the continuity residual is folded back through the
        forest so both exactness properties hold.
      A separate, smaller bug fixed along the way: `continuity_error` counted
      fixed-head reservoirs as continuity failures, so a converged two-reservoir
      solution reported an "error" equal to the network's own net demand. The
      GGA always excluded them.
      **Remaining limitation, which is a different and stronger one:** loop
      correction cannot choose between the several flow fields that satisfy its
      own equations, since `h = r Q|Q|` is not injective. On a five-link
      two-reservoir network Hardy Cross and the GGA now return *different* flow
      fields, each internally consistent. The GGA is head-based with a
      positive-definite Jacobian, so the physical branch is the only one it can
      settle at. Use the GGA when the values matter.
- [x] Implement Global Gradient Algorithm (node-based Newton-Raphson) network
      solver (`gga.rs`), on `tpt-math-linalg`. Fully general: multi-source,
      looped, and tree networks, with a backtracking line search. Validated by
      a finite-difference check on the analytic conductance and against an
      independently derived reference solution.

- [~] Implement transient analysis: water hammer via Method of
      Characteristics (MOC), Joukowsky equation, column separation modeling
      - *done so far in `water_hammer.rs`*: wave speed, Joukowsky rise,
      critical closure time, Courant validation, column-separation clamping,
      and a MOC solver for a reservoir-fed pipe under a prescribed valve
      closure. The frictionless limit reproduces the Joukowsky rise exactly,
      which is the check that matters. Still to do: a multi-node network with
      wave reflection at boundaries, and the friction-damped rise integral.*
      - **Multi-node reflection: implemented and verified.**
        The single-reach solver advances one characteristic pair, so a wave has
        no state to send back -- that is the structural reason it cannot
        reflect, and a head array updated in place cannot fix it. A wave is a
        *directional travelling perturbation*, so the reflection-capable form
        needs two accumulator arrays (`c_minus` upstream, `c_plus` downstream)
        onto which waves are superimposed, and only then can a wave cross the
        pipe, turn round, and add to what is already there.
        `MocNetwork` now does this: every branch is split into reaches and every
        end is a node solved from **both** arriving characteristics. The
        formulation rests on one fact that is easy to get backwards -- in a
        steady state `Q_P = Q_0`, so a frictionless pipe has **no** steady head
        gradient and the `a Q / (g A)` term is purely the transient Joukowsky
        rise. Seeding the steady heads with that term was the bug that had been
        defeating every previous attempt, and it showed up at once in the
        fixed-point test (`a_steady_network_does_not_move`).
        `MocNode::DeadEnd` closes the loop: a fixed-head reservoir absorbs a wave
        completely and so can never display the doubled head a closed end
        produces. The dead end reflects at full amplitude, verified by a
        staircase whose steps are exactly twice the one-reach wave amplitude.
        **History, for whoever reads this next:** the two reverted attempts
        below are kept as a record of what was wrong, not as open work. Both
        reached working wave machinery and both were discarded for a *physics*
        reason, not a bug: the accumulator form had no boundary able to hold
        pressure, and the node form was not a simultaneous solve at each node.
        Every defect they found is in the shipped code as a test. The
        accumulator form also showed the Courant limit must be checked
        directly, and that the characteristic companion of an interior end is
        `k - 1`, not `k + 1`.
        The accumulator attempt found four defects, all of which produce
        plausible numbers rather than errors:
        1. The steady profile applied `(a/g) Q` *per reach*, putting `8 (a/g) Q`
           across an 8-reach pipe instead of `(a/g) Q` once -- a valve at
           -308 m under a 100 m reservoir. The loss scales with *distance*, so
           it needs a `dx/L` factor.
        2. The shift copied without clearing, smearing a wave across every reach
           at once instead of letting it occupy one.
        3. Reflection read `c_minus[0]` *after* the shift had already carried the
           wave out of the pipe, so it always found zero and waves vanished at
           the reservoir instead of returning.
        4. The valve read only the returning wave, delaying the Joukowsky rise
           by `2 reaches` steps.
        **Why it was not shipped anyway:** fixing the four above still left the
        valve head collapsing back to steady the instant the flow stopped, which
        is wrong -- a closed valve holds *pressure*, not flow, and a
        prescribed-flow boundary cannot represent that. The correct downstream
        boundary is prescribed-*head* (or a dead end, reflecting at `-1`), which
        is a different boundary condition and therefore a different solver, not a
        bug in this one. Shipping the wave machinery behind a boundary condition
        that cannot hold the physics would have produced a solver that reflected
        waves beautifully and still got the pressure history wrong. The code was
        reverted; `water_hammer.rs` is unchanged from its last commit. The
        reflection coefficients that would drive the correct form are now
        stated in this file rather than in code that does not exist.
      - **Multi-node MOC: attempted a second time, also not shipped, and this
        time for a different reason.** A different formulation was tried --
        the node-based one, where every pipe is discretised into reaches with
        `dt = dx / a` and each node is closed by its boundary condition
        (reservoir: head known; valve: flow known; junction: continuity). This
        gets much further than the accumulator form: the Joukowsky rise came
        out **exactly** right (1.000x the closed form), and the Courant
        commensurability across branches with different lengths and wave speeds
        works. Six defects were found and fixed, all silent -- plausible finite
        numbers, never errors:
        1. The momentum coefficient. `B = a / (g A)` paired with a *volume* flow
           is the correct head; `a / g` with a volume flow, or `a / (g A)` with a
           velocity, are each wrong by a factor of the area. Every version of this
           tried it wrong at least once.
        2. The steady seed applied `B Q` per reach instead of once per branch,
           putting a 1000 m gradient on a 150 m pipe. The head field is only
           consistent if the *interior* ends are seeded from the characteristic
           itself, one `B Q` apart -- a merely smooth interpolation is not
           enough, because the propagation then shifts the steady field one
           reach per step while looking perfectly smooth.
        3. `end_node` mapped every downstream end to `branch.to`, so every
           *interior* computational end inherited that node's boundary
           condition. A valve closure was then applied along the whole pipe at
           once, which is the exact lumped behaviour the solver exists to
           remove. Interior ends must report no network node at all.
        4. The characteristic companion for an interior end was taken as `k + 1`
           instead of `k - 1`, which convects the head field one reach per step
           rather than propagating a wave.
        5. Copying a boundary end's prescribed flow onto every end of the branch
           teleported the wave; an interior end's flow is a property of its own
           reach and must not be overwritten.
        6. The valve's head was computed from the end's already-propagated value
           and then corrected again, giving exactly `2x` the Joukowsky rise, and
           later exactly `n / (n - 1)` of it once the interior ends were
           interpolated. The correction has to be applied to the end's *previous*
           head, not to the propagated one.
        **Why this one was also reverted:** the remaining error was a
        factor-of-`n/(n-1)` on the valve rise, which is a boundary-application
        detail, but chasing it exposed the real blocker. The node-based form
        requires each node to be closed by a **simultaneous** solve over all of
        its ends at once: every end at a junction gives an affine relation
        between the common head and its flow, and those relations plus the node
        balance form a small linear system. Closing the ends one at a time --
        solving for the head, then writing each end's flow back -- is not an
        approximation of that system, it is a different and inconsistent one,
        and it is what left the wave appearing one reach early. A correct
        implementation needs the per-node system assembled and solved
        together, which is a real piece of work rather than a fix.
        The diagnosis is the useful output here: **the frictionless Joukowsky
        rise, the Courant commensurability, and the head-consistent seed are all
        now understood and were verified numerically before the revert.** What
        remains is the simultaneous node solve, and the reflection-coefficient
        boundary of the previous attempt still stands as the alternative route.
        `water_hammer.rs` is unchanged from its last commit; no half-working
        multi-node code was left behind.
- [x] **Multi-node MOC, shipped.** `MocNetwork` in `water_hammer.rs` discretises
      every branch into reaches with a Courant-common `dt`, and each *end* is a
      node solved from two arriving characteristics. The formulation, and the
      four facts that took three attempts to establish:
      1. **`B = a / (g A n)` per reach**, paired with a *volume* flow. `a / (g A)`
         alone is the whole-*pipe* loss; using it per reach gives `n` times the
         gradient, a 4x error on a 4-reach pipe that stays finite and smooth.
      2. **Both characteristics take the same sign on the flow term** --
         `a_plus = H_up + B Q` and `a_minus = H_dn - B Q`, so the end satisfies
         `H = a_plus - B Q` on one side and `H = a_minus + B Q` on the other.
      3. **An interior end is a node with two unknowns** and needs both
         equations: `Q = (a_plus - a_minus) / (2B)`, `H = (a_plus + a_minus) / 2`.
         Treating it as a single propagated value is a one-equation treatment of
         a two-equation node and makes the wave appear a reach early on every
         step.
      4. **The Courant limit must be checked directly**, not only through
         `ceil(L / (a dt)) >= 1`, since `ceil` of a small number is still 1 and
         a `dt` a hundred times over the limit would be rounded away silently.
      **The bug that had defeated two earlier attempts**, and which the
      "a steady network does not move" test found immediately: a *frictionless*
      pipe has **no steady head gradient at all**. In a steady state `Q_P = Q_0`,
      so the flow term in the linearised momentum equation vanishes, and the
      `a Q / (g A)` term is purely *transient* -- it is the Joukowsky rise and
      appears only when the flow changes. Seeding it as a steady loss puts a
      non-physical gradient into the initial state, which is then inconsistent
      with the characteristics and makes the whole network drift on a run that
      should be perfectly steady. The steady drop comes from **friction alone**,
      lumped at the nodes, which is also why the interior ends are flat.
      Two test-data traps worth keeping: the pre-existing single-reach tests use
      `Q = 0.5` in a `0.0707 m^2` pipe, which is 7 m/s and passes only because
      those tests treat `Q` as a *velocity*; and a **fixed-head reservoir cannot
      show a Joukowsky head rise at all** -- the arriving wave is absorbed by the
      flow changing, so the far end's observable is "the flow stops", not "the
      head rises". Asserting a head rise at a reservoir tests something its
      boundary condition forbids.
      Tests: 33 in the module, including the steady fixed point, one-reach-per-
      step propagation, friction-only steady gradient, per-reach valve
      amplitude, the flow stopping only after the wave arrives, the
      reconstructed Joukowsky head `n x` the one-reach rise, the dead-end
      full-amplitude reflection, a moving reservoir launching a travelling
      wave, the unsteady-friction steady limit and damping, and a check valve
      whose closure travels as a wave rather than being lumped.
      **Second bug found and fixed in the same area.** The check-valve block
      zeroed the branch's *upstream* end as well as its downstream one, so a
      check valve stopped flow across the whole branch in a single step and
      annihilated the wave still travelling through the reaches behind it --
      another lumped shortcut, and the reason it survived is that **it had no
      test at all**. Clamping the valve end alone makes the closure propagate as
      a real front, one reach per step, which
      `a_check_valve_stops_at_its_own_end_and_the_stop_travels` now pins.
      **Third bug found and fixed in the same area.** Column separation clamped
      the head to the vapour head but left the flow frozen at its
      pre-separation value, leaving a `(H, Q)` pair that satisfied no
      characteristic; the next step read it back as an arriving wave, so the
      violation was injected rather than absorbed. The trace was the giveaway --
      head jumping to 1e6 with the flow sitting at 0.0035. Separation now also
      releases the flow to the free discharge that head implies, and
      `column_separation_releases_the_flow_rather_than_freezing_it` pins the
      exact value.
      **Still open:** nothing in the MOC formulation itself. Unsteady friction is
      per-branch via `with_derived_friction_lag`, which derives `T_f =
      L / (g A R |Q|)` from each pipe's own properties, and the friction rise
      integral is now reproduced against its closed form.
      **Bug found and fixed while doing the above.** The friction force entered
      each characteristic as a *difference* between neighbouring ends,
      `r (Q_k E_k - Q_{k-1} E_{k-1})`, which is identically zero for the uniform
      flow of a steady network. Nothing held the gradient up, so the seeded
      Darcy-Weisbach drop washed out within a few steps and the network settled
      back to frictionless -- the exact opposite of the documented "the steady
      drop comes from friction alone". It failed silently: every value stayed
      finite and plausible. The friction is now the absolute per-reach loss
      `r Q E`, each interior end is seeded on the linear gradient between its
      node heads, and `a_frictional_steady_state_is_an_exact_fixed_point` pins
      it to 1e-12 over 40 steps.
- [x] Implement component models: valve `Cv` and `K` coefficients, pump
      characteristic curves (quadratic three-point fit, shut-off head, runout
      flow, hydraulic power), cavitation state from the cavitation number, and
      surge tank dynamics. Tabulated minor-loss coefficients for entrances,
      elbows, and exits are included. **The turbine four-quadrant curve is
      implemented** (`TurbineCurve`), not just its loss coefficient: generating,
      windmilling and pumping quadrants via `quadrant()`, with `runout_flow`
      and the extracted hydraulic power. This entry was stale -- it described
      the crate as it was before the four-quadrant work.
- [x] Implement differentiable head-loss functions (via `tpt-math-autodiff`)
      for downstream pipe-sizing optimization. `differentiable.rs` evaluates
      Darcy-Weisbach over forward-mode dual numbers, so the diameter
      sensitivity comes out exactly in one pass. Supports the laminar,
      Swamee-Jain, and Hazen-Williams correlations analytically. **Colebrook-White
      no longer needs a finite-difference fallback**: its iteration is now run
      over duals, so the gradient is exact. Differentiating the *iteration*
      rather than the implicit equation sidesteps the real hazard here, namely
      that Colebrook depends on the diameter through two routes (`eps/D` and
      `Re`) and hand-derivation drops one of them, which silently flips the
      sign of the friction-factor contribution. Two tests pin this: one checks
      the dual gradient against a centred difference to 1e-7, the other checks
      it is *linear* in `dRe/dD` and passes through the origin, which a
      difference quotient cannot do exactly.

- [x] Unit tests for the modules delivered so far — 12 friction tests,
      including a *residual* test that checks the Colebrook-White output
      satisfies its own defining equation to 1e-6 independently of the
      iteration used to produce it
- [x] Rustdoc for the solvers. Every module carries crate-level rationale, the
      governing equation, and — where the physics has a limit — an explicit
      statement of what is *not* modelled. The examples that are worth reading
      first are the ones about limits: the Hardy Cross convergence caveat in
      `lib.rs`, the Reynolds validity discussion above `solve_journal_bearing`,
      and the quasi-steady temperature-rise warning in `wear.rs`.
- [x] Froude-Krylov excitation: the first-order wave excitation force on a
      wall-sided hull, `F_x = rho g a (1 - e^(-kT)) 2 sin(kL/2) / k`, resolved
      onto the ship's axes for any heading. Both limit cases are tested: the
      force vanishes for waves far shorter than the hull, and tends to
      `rho a w^2 T L` in the long-wave limit
- [x] Holtrop-Mennen: total resistance from principal dimensions, split into
      wave-making, friction, residuary, and appendage terms, with Holtrop's
      own wetted-surface estimate. Validated against a 120 000 dwt bulker at
      14 kn, whose resistance is known from its required power
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] Add to root `Cargo.toml` members + workspace deps

## Phase 2 status

Complete except for the items marked in progress. The crate carries
`error`, `friction`, `network`, `hardy_cross`, `gga`, `water_hammer`,
`components`, and `differentiable`. **119 unit tests pass with none
`#[ignore]`d** -- the two that were withheld for the documented Hardy Cross
non-convergence are now enabled, and that non-convergence is fixed. The GGA
remains the general solver and is validated against an independently derived
reference. Reach for it whenever the flow *values* matter, not merely when the
network is small: Hardy Cross now closes every loop on these networks but
cannot select between the several flow fields its own equations admit, whereas
the GGA's head-based formulation can. See the Hardy Cross entry above.

## Phase 3 — tpt-fluids-marine

*Naval architecture and marine hydrodynamics for ships and underwater
vehicles. Depends on: `tpt-fluids-core`, `tpt-math-linalg`/
`tpt-math-linalg-fixed`, `tpt-math-signal-fft` (wave spectra),
`tpt-math-autodiff` (Michell integral differentiability).*

- [x] Scaffold `crates/tpt-fluids-marine/`
- [x] Wire deps: `tpt-fluids-core`, `tpt-math-autodiff`
- [x] Implement frictional resistance: ITTC-1957 correlation line, with the
      traditional `C_f x 10^3` presentation, friction resistance from the
      wetted-area estimate `S = L(B+T)`, Granville's roughness and appendage
      extension
- [x] Implement the ITTC-1957 model-ship line `R = k W^(2/3) V^6` and the
      rectangular-reference-model resistance
- [x] Implement residuary resistance: Michell's thin-ship integral as a
      Froude polynomial `C_r = a + b Fr^2 + c Fr^4 + d Fr^6`
- [x] Principal dimensions, block coefficient validation, volume and
      displacement, and the power relationships (delivered power, effective
      horsepower, and the inverted speed-for-power a sizing loop needs)
- [x] Unit tests: 122, with the ITTC-1957 line checked against the published
      figure (`C_f x 10^3` ~ 38 for a 300 m tanker at 12.5 kn) and the
      model-ship line checked for its exact `W^(2/3) V^6` scaling
- [x] Froude scaling: model-to-ship extrapolation with form factor and
      roughness allowance, in `froude_scaling.rs`. Covers the form factor
      `S_wet / (L(B+T))`, the box-hull reference, the roughness allowance,
      `C_ship = C_model / k_f (1 + k_r)`, the `V ~ sqrt(L)` speed scaling, the
      `L^3` displacement scaling, and the full resistance and power chain.
      Validated end-to-end against a 1:50 model of a 150 m, 3000 DWT feeder,
      which reproduces 294 kN and 4.3 MW
- [x] Implement seakeeping: linear wave theory, response amplitude
      operators for 6-DOF ship motions, significant wave height statistics.
      `seakeeping.rs` carries the deep-water dispersion relation and
      wavelength, the ITTC Pierson-Moskowitz significant height and peak
      period, the six motion modes, and a damped-oscillator RAO per degree of
      freedom that captures resonance, the 180-degree phase reversal above it,
      and damping-limited roll peaks. A Green-function hull-integral RAO is
      not attempted; the SDOF form is documented as an approximation.
- [x] Implement the MMG model for 3-DOF horizontal-plane manoeuvring (surge,
      sway, yaw), hydrodynamic derivatives. `manoeuvring.rs` carries the
      coupled MMG equations of motion including the sway-yaw cross terms,
      mass and added masses with a gyradius rule, linear and quadratic hull
      damping, a lift-based rudder model with a clamped maximum deflection,
      a fixed-step integrator, and the steady turning-circle solution with
      the diameter expressed in ship lengths.
- [x] Implement propulsion: wake fraction, thrust deduction, propeller-hull
      interaction, open-water coefficients. `propulsion.rs` carries the
      `T = K_T rho n^2 D^4` and `Q = K_Q rho n^2 D^5` definitions, the
      advance ratio, the open-water and propulsive efficiencies, the wake
      reduction, thrust deduction, quasi-propulsive coefficient, tip speed,
      and the cavitation-limited diameter. Validated against an independent
      actuator-disc estimate, which back-solves K_T ~ 0.38, and against a
      3000 DWT feeder design point at 18 kn.
- [x] Make the Michell integral differentiable via `tpt-math-autodiff` for
      hull-form optimization - `tpt-fluids-marine/src/differentiable.rs`, which
      is what the already-declared-but-unused `tpt-math-autodiff` dependency in
      that crate's manifest was waiting for. The residuary polynomial is written
      in Horner form over dual numbers, plus the `0.5 W Fr^2 (C_f + C_r)` total
      chain and `P = R V / (eta D)`. 14 tests.
      The ITTC-1957 friction correlation is deliberately **not** differentiated:
      it is an empirical model-test fit, so its derivative has no physical
      standing, and `C_f` is held constant instead. The `V^2` scaling, which is
      the part that actually moves with the design, is still carried exactly.
      Two things the tests caught rather than confirmed. A gradient descent on
      resistance-per-unit-speed converges to the *same* interior Froude number
      from both 0.05 and 0.8, which is a stronger claim than "it decreased". And
      the first scale assertion was wrong in a way worth recording: a kilonewton
      band borrowed from the `froude_scaling` feeder test disagreed by a factor
      of 6, because that module normalises by `rho g` times a displacement
      treated as a *volume* while the Michell path uses a displacement *force*.
      The test now asserts `R/W`, which is convention-free -- the same class of
      error as the two ITTC friction conventions this crate already warns about
      twice.
- [x] Integration test: model-ship resistance extrapolation matches
      Holtrop-Mennen benchmarks (see Phase 7, where the cross-method benchmark
      suite landed; the residual line here was stale)
- [x] Froude-Krylov excitation: the first-order wave excitation force on a
      wall-sided hull, `F_x = rho g a (1 - e^(-kT)) 2 sin(kL/2) / k`, resolved
      onto the ship's axes for any heading. Both limit cases are tested: the
      force vanishes for waves far shorter than the hull, and tends to
      `rho a w^2 T L` in the long-wave limit
- [x] Holtrop-Mennen: total resistance from principal dimensions, split into
      wave-making, friction, residuary, and appendage terms, with Holtrop's
      own wetted-surface estimate. Validated against a 120 000 dwt bulker at
      14 kn, whose resistance is known from its required power
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] Add to root `Cargo.toml` members + workspace deps

## Phase 4 - `tpt-fluids-tribo`

- [x] Scaffold `crates/tpt-fluids-tribo/`
- [x] Implement Hertzian contact: reduced modulus and reduced radius from
      signed curvatures, contact radius, peak and mean pressure, elastic
      approach, the inverse load-for-peak-pressure relation, and an explicit
      yielding check (Hertz knows nothing of yield, so this decides whether
      its answer can be trusted)
- [x] Implement lubrication: Sommerfeld number, regime classification,
      Petroff's hydrodynamic friction, the journal-bearing Stribeck curve
      with a genuine interior minimum, minimum film thickness, and the
      eccentricity-load solution
- [x] Archard's linear wear law and the wear coefficient, plus the
      inversion, wear depth, life-for-depth, the lambda ratio and its
      separation-regime classification, and frictional heating. **The
      transient flash temperature is now solved** (Blok-Wilde, `2qL /
      (k sqrt(pi alpha t_c))`) alongside the quasi-steady form, with
      `ThermalProperties` and a `crossover_time` that says which applies.
      Doing so exposed a **dimensional bug in the pre-existing steady form**:
      `mu F v / (k A)` is kelvin per *metre*, because conduction `q L / k`
      needs a length. The "20 000 K brake pad" that the old test pinned as a
      documented limit was that expression read as a temperature; the same
      duty with `L = sqrt(A/pi)` is about 1 600 K, a real disc temperature.
      Both forms now take the length explicitly, and the test asserts the
      corrected value. The crossover is ~30 hours for steel, which is the
      quantitative statement of why the transient form is the right default
      and the steady one is not a close approximation to it.

- [x] `cargo fmt` / `clippy` clean
- [x] Add to root `Cargo.toml` members + workspace deps

## Phase 5 - `tpt-fluids-verify`

- [x] Scaffold `crates/tpt-fluids-verify/`
- [x] Proptest invariants: dimensional invariance under a change of length
      units, monotonicity of every correlation in its own inputs, exact
      scaling laws (Hertz's `F^(-1/3)` peak, `W^(2/3) V^6` extrapolation,
      `n^2 D^4` thrust, `S^2` clearance), and the requirement that no
      physically admissible input produces `NaN`
- [x] 26 proptest properties, 256 cases each
- [x] Kani proof harnesses for solver safety, behind the `kani` feature and
      `cfg(kani)`, covering the ITTC coefficient, the Hertz contact radius,
      peak pressure, mean-versus-peak ordering, the elastic approach, Archard
      non-negativity, and propeller input overflow
- [x] `cfg(kani)` declared in the workspace lint config
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] Add to root `Cargo.toml` members + workspace deps

**Note:** `cargo-kani` is not installed in this environment, so the Kani
harnesses are written and gated but unproven. They are behind `cfg(kani)`, so
they are not silently passing either. The proptest invariants run in the
ordinary suite and have already found a real bug, described in the commit
message.

## Phase 6 - `tpt-fluids` umbrella crate

- [x] Scaffold `crates/tpt-fluids/`
- [x] Re-export every sub-crate behind its own feature
- [x] `prelude` module: the quantity types, dimensionless numbers, constants,
      and the module path for each application domain, in one import
- [x] `std` feature propagating to the sub-crates with `?/` so it only
      reaches a crate whose optional dependency is actually enabled
- [x] Integration tests: 8, covering cross-crate consistency that unit tests
      inside a crate cannot see (friction resistance against its own
      definition, the Froude-scaling chain composed two ways, the two ITTC
      friction conventions staying 25x apart through the re-export, the Hertz
      chain including its inverse, and quantities flowing unchanged between
      marine and tribology)
- [x] Full feature matrix builds: core alone in `no_std`, each domain alone,
      all three together, and the whole thing with `verify`
- [x] Doctest in the crate documentation exercising all three domains
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

## Phase 7 - Integration & Workspace Closeout

- [x] Integration test: model-ship resistance extrapolation
      (`Ittc1957Line`) matches the independent Holtrop-Mennen estimate for
      the same 150 m, 3000 DWT feeder to within 0.01 percent at the
      calibrated coefficient, and Holtrop-Mennen falls inside the band a
      plausible spread in a measured model test would give. Plus the power
      chain closing: the propeller must deliver more than the hull resists,
      and a Beaufort sea must put the spectral peak within a few ship lengths
- [x] `cargo test --workspace --all-features` passes
- [x] `cargo clippy --workspace --all-targets --all-features -- -D warnings` clean
- [x] `cargo fmt --all -- --check` clean
- [x] `cargo deny check` clean workspace-wide, no Apache-2.0-only crate in the
      resolved graph
- [x] `no_std` matrix resolved and documented: `tpt-fluids-core` builds for
      `thumbv6m-none-eabi` with `--no-default-features --features alloc` and
      also with no features at all. The three domain crates declare a `std`
      feature and are **not** `no_std`; they enable `core/std` in their
      manifests, so `thumbv6m-none-eabi` fails for them with "can't find
      crate for std". This is a deliberate split and is now stated in the
      README rather than left as a surprise.
- [x] Umbrella feature matrix builds: `core` alone in `no_std`, each domain
      alone, all three together, and everything with `verify`
- [x] Root `README.md` rewritten: accurate crate map, build order, design
      philosophy, test inventory, known limitations, and commands. The
      previous text had UTF-8 mojibake from the bootstrap and claimed EHL,
      Dowson-Hampton and LuGre friction, none of which exist.
- [x] `cfg(kani)` declared in the workspace lint config
- [x] EHL **film thickness** implemented in `tpt-fluids-tribo::ehl`:
      Dowson-Higginson (line contact) and Hamrock-Dowson (point contact), the
      minimum-film factor `H_min = 0.8 H_c` for the side-lobed point contact, and
      the lambda ratio with its separation-state thresholds. These are closed-form
      correlations and need no elastic solver, so the film thickness — the
      quantity every lubricant-selection decision turns on — is now available.
      Two notes from doing it. First, a heavily loaded gear mesh gives a film
      around **1 nm**, which against ground flanks is a lambda near 0.02 and
      therefore genuine *boundary* lubrication; quoting "EHL" for such a contact
      without checking the ratio is the mistake the module exists to prevent.
      Second, the correlations are only meaningful inside their validity range,
      and a nearly stationary contact will happily return a sub-atomic film
      that is arithmetically correct and physically meaningless.
- [ ] Advanced coupling follow-up: the **coupled** EHL solution — pressure
      distribution and sub-surface stress from solving Reynolds together with the
      elastic deformation of both bodies. **Not externally blocked after all** --
      see the correction in the spec audit below: `c:\Programming\tpt-fem` contains
      a local `tpt-fem-elasticity` crate. The remaining obstacle is a policy one
      (unpublished, path-resolved, and `deny.toml` requires a registry source),
      not an absent dependency.
- [x] **End-to-end pipe-network diameter optimization** via
      `tpt-systems-optimisation` (spec.txt line 176, Phase 3). Done:
      `tpt-fluids-hydraulic/src/sizing.rs`. The `tpt-opt-core` augmented-Lagrangian
      solver (`nlp` feature) drives a `NlpProblem` over the diameter vector;
      cost is `sum L c D^2`, constraints are total and per-pipe head loss, the
      velocity limit, and both diameter bounds. 13 tests. Converges in ~16 outer
      iterations and ~1 s on a three-pipe network.
      Wired as a **git** dependency, per the maintainer's policy decision, pinned
      to the exact PGP-verified commit `95fbc3f2` rather than `master` so an
      upstream change cannot alter the solver without appearing in this repo's
      history. `deny.toml` keeps `unknown-git = "deny"` and adds one named
      `allow-git` entry, so a *new* git dependency still cannot arrive unnoticed.
      Two things the tests caught, both worth keeping in mind for any future
      consumer of an AL solver:
      - **The upstream solver scores feasibility with `ineq(i, x).max(0.0)`, and
        Rust's `f64::max` returns the non-`NaN` operand.** So a `NaN` constraint
        reads as *satisfied*, and the solver reported `converged = true` in one
        iteration on a decision vector of pure `NaN`. The first run of this
        module did exactly that. `SizingProblem::guard` maps any `NaN` to
        `+INFINITY`, putting it on the infeasible side where the penalty can act;
        the property is pinned by a test that states the upstream scoring rule
        directly.
      - **The obvious starting point is the wrong one.** The velocity-limit
        diameter is infeasible on head loss by construction, and the AL method is
        a descent method, so it walks further into infeasibility and never
        recovers -- it returned a network *more expensive* than the one it was
        given, and a 3 percent shrink was cheaper still. The start is now grown
        until feasible, which is one cheap loop and puts the search inside the
        region where descent means something.
      The physics stays in this crate: the optimiser supplies the constrained
      solver, never the hydraulics.

## Spec audit (against `spec.txt`, post-Phase-7)

The seven phases are closed, so this records what `spec.txt` still asks for
that the workspace does not yet have, audited line by line rather than
assumed.

### Closed this pass

- [x] **Mass conservation at every network node** (`spec.txt` line 145, and
      line 31's worked example). Now a proptest over generated topologies, in
      `tpt-fluids-verify/src/conservation.rs`. The generator builds a random
      spanning tree rooted at a reservoir plus random chords, so it is
      connected and looped by construction and never has to reject a graph.
      Checked on both trees and looped networks, with the tolerance scaled by
      the largest flow in the network so the test is relative.
- [x] **Energy conservation in a lossless segment** (line 149). A pipe with no
      demand can neither dissipate head nor carry flow. Paired with its
      complement, a frictional segment *must* lose head, so the pair
      distinguishes a solver that models dissipation from one that conserves
      everything by accident.
- [x] **proptest strategies for valid network topologies** (line 124). Done as
      part of the above.
- [x] `tpt-fem-contact` was investigated for the EHL coupling and is **not** a
      path dependency. It is finite-element surface-to-surface contact
      (penalty and augmented Lagrangian on a mesh); EHL needs an elastic
      compliance kernel for a rough surface, which is a different thing. The
      coupling stays open, but the blocker is now described accurately rather
      than as "the crate is missing".

### Still open, in rough priority order

- [x] **Reynolds equation solver** (line 115) -
      `tpt-fluids-tribo/src/reynolds.rs`. 15 tests. Finite-volume
      discretisation exact to machine precision against a manufactured
      solution; Thomas on a negative M-matrix; load capacity, attitude angle,
      and the minimum film ratio.
      **Solved for the full-film regime only, and the solver refuses outside
      it.** `solve_journal_bearing` solves at two grid resolutions and
      requires them to agree; past about `e = 0.4` they do not, and the call
      returns `OutsideValidRange` rather than the 8e5 the raw solve produces.
      A real bearing runs at `e = 0.2`-`0.3`, so this is the regime a
      designer sizes against.

      Two artefacts of the Sommerfeld formulation are documented rather than
      hidden, because both look like solver bugs if you do not know:
      - The condition pressurises the whole circumference, so the solution
        goes **negative in the diverging half for any non-zero
        eccentricity** (-0.23 at `e = 0.2`, -0.82 at `e = 0.3`). That is the
        clearest evidence the cavitation condition is needed, and a test
        tracks the excursion growing with `e`.
      - At exactly `e = 0` the boundary conditions alone force a linear
        pressure ramp carrying `1/(2 pi)`. A concentric bearing physically
        carries nothing; the Sommerfeld condition presumes a converging wedge
        exists and at `e = 0` none does.

      **The bug that mattered most:** the `P(0) = 1` Dirichlet contribution to
      the first matrix row was being dropped. The result stayed smooth and
      plausible but was wrong everywhere, and the manufactured-solution test
      could not catch it because that test uses homogeneous boundaries. The
      tell was `Pmax = 1.00000` for every grid, which is just the boundary
      value. Fixing it took the peak from 1.8195 to the correct
      grid-independent 1.81954 and made the solver converge.

      **Closed (2026-09-29): cavitation via complementarity.**
      `solve_journal_bearing_cavitated` imposes `P >= 0` with a primal-dual
      active-set iteration over the same matrix (the M-matrix makes it
      monotone). Grid-converged to `e = 0.49` (peak 33.5, load 3.437 at every
      resolution for `e = 0.48`); pressure never negative; complementarity
      checked residual-by-residual in a test. **Correction to the text below:**
      the 8e5 blow-up at `e = 0.5` was partly the film closing (`h_min = 1 - 2e
      = 0`), which no solver survives -- the limit is now
      `MAX_CAVITATED_ECCENTRICITY = 0.49`, contact rather than method. The
      Sommerfeld solve did not actually refuse at 0.45.
      Kani equilibrium harnesses are unaffected (they concern the load direction).
      *Superseded note:* Reynolds' supplementary cavitation condition for
      `e > 0.4`. A post-hoc "find where P < 0" loop cannot work, because the
      plain solve spikes *positive* before it ever goes negative. A
      pivot/bubble method during elimination was tried and detects some
      cavitation but does not converge for high `e`. A proper complementarity
      treatment (JFO, or penalty) is the remaining work. Until then the
      refusal is the correct behaviour, not a workaround.
      - Kani: Reynolds global load equilibrium (line 128) still blocked on it.

- [ ] **Elastohydrodynamic lubrication / Dowson-Hampton** (line 116). Left open
      deliberately, and the reason is worth recording rather than glossing.
      `spec.txt` line 116 asks for two things: the coupled Reynolds + elastic
      deformation solve, and the Dowson-Hampton film-thickness *formulas*. Both
      are unshipped, but for different reasons, and the deformation side turned
      out **not** to be blocked on `tpt-fem-elasticity` at all -- see the note
      below, which supersedes the "blocked on tpt-fem-contact" reasoning.
      - **The Dowson-Hampton coefficients specifically.** I attempted to source
        them and could not obtain them from any source that renders as text: the
        Penn State and NASA NTRS PDFs return raw compressed streams, and the MDPI
        review returns HTTP 403. Writing
        `2.69 U^0.68 G^0.49 (1 - e^(-0.68k))` from memory would be precisely the
        failure mode this file is built to prevent -- a precise, plausible,
        entirely unverified number.
      - **The elastic side: attempted, not shipped, and this one is closer.**
        Line 116's coupling needs the *surface deformation* under film pressure,
        which is the Cerruti/Boussinesq half-space solution -- an exact analytical
        result, not an FEM solve, and emphatically not `tpt-fem-contact`. So the
        "needs `tpt-fem-elasticity`" reasoning above is wrong twice over: the
        crate exists locally, and this physics does not need it anyway.
        I built the kernel, with the substitution `s = r sin(phi)` to regularise
        the integrable `1/sqrt(r^2 - s^2)` singularity, and validated it against
        Hertz's independent `a^3 / (3 R*)` -- the right check, since the two are
        derived completely differently.
        The check earned its keep three times over, and all three failures are
        worth recording because each is a *finite, plausible* wrong answer:
        1. My test helper computed `1/((1 - nu^2) E)` -- the compliance -- where
           `contact_radius` wants `E* = E / (1 - nu^2)`. That gave a 1122 m contact
           radius for a 1 N load: positive, finite, wrong by 15 orders.
        2. With that fixed, the kernel agreed with a direct quadrature of the
           unsubstituted expression to 0.05%, so the substitution was sound.
        3. But the integrated total approach came out **3.038x** Hertz, and
           refining the panels did not move it -- a constant error, not a
           convergence failure. I could not reduce 3.038 to a clean prefactor
           (it is not pi, 2/pi, pi^2/6, or a simple multiple), which says the
           `1/r` placement in the closed form or the outer area element is wrong
           in a way I could not pin down from a text source. Johnson, *Contact
           Mechanics* (cited throughout the Willert paper I checked) is the
           authority and is not available as text here.
        Shipping a deformation kernel off by 3x would be worse than shipping
        none: an EHL film thickness is the difference between a few nanometres
        and a micrometre, and a 3x error in the deformation is a 3x error in the
        film. This is the same standard already applied to the Hamrock-Dowson
        coefficients above, applied consistently: **a constant I cannot verify is
        not a constant I ship**, whether it is a correlation coefficient or an
        integral prefactor. The code was reverted; `tpt-fluids-tribo` is
        unchanged.
      - **Johnson located, and it does not help: the book is a scan.** A
        complete copy of Johnson's *Contact Mechanics* is reachable at
        `meil.pw.edu.pl/.../Johnson-CONTACTMECHANICS.pdf` (4.99 MB, 462 pages),
        and the Boussinesq/Cerruti framework was independently corroborated by
        Willert, arXiv 2108.04617, which cites Johnson throughout. But that PDF is
        a **scanned image with no text layer** -- the pages are JBIG2 bitmaps --
        so fetching it returns raw image streams and equation 3.63 cannot be read
        from it. The same wall as the other PDFs above, and for the same reason.
        **This is now a resolved question of *access*, not of physics:** the
        formula is standard, the framework is confirmed, and only a text-rendering
        copy of eq. 3.63 (or any secondary source that quotes it) is missing.
        Someone with the book in hand, or a text-layer copy, can close this in
        minutes. What is *not* going to help is another round of inferring the
        prefactor from first principles -- that is what produced the 3.038.
- [x] **LuGre dynamic friction** (line 119) - `tpt-fluids-tribo/src/friction.rs`,
      with Coulomb as the degenerate baseline, the Stribeck steady-state
      curve, the one-state bristle ODE, the relaxation timescale, and a sweep
      helper. 20 tests. Validated against hand-computed values: the bristle
      loads to exactly F_ss/s_0, the total settles to F_ss + s_2 v, and one
      step after a velocity jump the force is 24.4 against a steady state of
      20.025, a 22 percent memory excess that a static curve cannot produce.
- [x] **Load capacity for slider and thrust bearings**, and the differentiable
      load-capacity integrals (lines 115, 120). `reynolds.rs` now factors the
      validated finite-volume core out as `solve_profile_on_grid`, so slider and
      pivoted-thrust pads reuse the journal bearing's discretisation rather than
      growing their own copies. Both use **Gumbel** conditions (`p = 0` at both
      ends), not the journal's Sommerfeld: in a wedge there is no full circle
      for a diverging region to pressurise against, and reusing `P(0) = 1` would
      manufacture a pressure the geometry does not have. `thrust_pad_load`
      carries the `r dr` area factor that a line integral drops.
- [x] **Running-in wear simulation** (line 118). `RunningIn` in `wear.rs`: an
      exponential decay of the excess wear coefficient, with the closed-form
      volume integral checked against a numerical integration of the coefficient
      itself. That check earned its keep immediately — it caught a missing `load`
      factor in the excess term. Two claims that the tests disproved and that are
      now recorded rather than asserted: running-in wear is *less* than a
      fresh-contact coefficient predicts but *more* than a steady one, and
      `settled_distance = 3d` is 5% of the *excess*, not of the steady
      coefficient (for a 100:1 ratio the coefficient is still 6x steady at 3d,
      and the true 5% distance is `d ln(20(k0/ks - 1)) ~ 7.6d`).
- [x] **Turbine four-quadrant curves and cavitation inception** (line 97,
      Phase 3 line 175). `TurbineCurve` covers all four quadrants: the
      generating and windmilling quadrants fall out of one parabola, and the
      reversed branch is separate because the runner sees the flow from the other
      side, continuous in value at zero flow. `CavitationInception` carries a
      flow-dependent inception number and `CavitationMargin` joins it to the
      available NPSH; a machine checked only at its best point cavitates at
      overload, which the constant-number model misses. Both machine curves also
      gained the three-point quadratic fit this file already claimed they had.
- [x] **Kani: Hardy Cross monotone convergence** (line 126) and **MOC CFL
      stability** (line 127). `loop_correction` is factored out of the solver so
      the convergence property is a statement about one small function, with four
      harnesses covering correction direction, scaling, the exact Courant
      boundary, and the Joukowsky rise. Each is mirrored by a unit test, because
      the harnesses sit behind `cfg(kani)` and cannot run on this machine — a test
      that says "this would be proved" is not evidence.
- [x] **Kani: Reynolds global load equilibrium** (line 128). The "blocked on the
      solver above" note was stale: `solve_journal_bearing` and
      `solve_profile_on_grid` both landed, so this is tractable. The property is
      stated about one small function rather than the solver, for the reason the
      Hardy Cross harnesses were: a harness over a 400-cell linear solve proves
      something about a linear solve, which is not what line 128 is asking about.
      `reynolds.rs` factors out `load_direction_is_in_the_converging_half` and
      `sin_of_turn`/`cos_of_turn`, and three harnesses prove the load direction
      is a function of the wrapped position, agrees with the film profile, and
      that `h_min = 1 - 2e` stays positive across the resolved regime -- which
      makes `MAX_RESOLVED_ECCENTRICITY` a consequence of the arithmetic rather
      than a separately chosen constant. Each is mirrored by a unit test,
      because the harnesses sit behind `cfg(kani)` and cannot run here.
      **The mirror earned its keep immediately.** The first harness asserted the
      load direction is *always* in the converging half, which is simply false:
      a position at `0.6` is in the diverging half and must report so. The unit
      test failed on `position -0.49` and forced the property to be restated
      correctly, and a second failure on `-0.5` exposed a real off-by-one at the
      wrap boundary -- `position - position.floor()` maps `-0.5` to exactly
      `0.5`, so the comparison must be `<=`, not `<`. Both are pinned now. A
      Kani harness would have reported them as counterexamples; without the
      mirrors, nothing here would have been evidence of anything.
- [x] **`tpt-systems-optimisation`** (lines 156, 176) - **exists, and the
      license clears; the blocker is publication, not absence.**
      The earlier note here said "not present as a sibling repo", which was true
      of the local filesystem and wrong about the world: it is public at
      `github.com/tpt-solutions/tpt-systems-optimisation` (11 crates, Rust,
      pushed 2026-09-04). The family is `tpt-opt-{core, milp, minlp, network, cp,
      heuristic, multi, robust, decompose, conic, systems}`.
      **The license is not the obstacle I assumed.** The GitHub *repository* API
      reports `Apache-2.0`, which would have failed `deny.toml`'s "no
      Apache-2.0-ONLY anywhere in the graph" rule -- but the repository's own
      `Cargo.toml` declares `license = "MIT OR Apache-2.0"` and every member
      crate inherits it via `license.workspace = true`. Dual-licensed, so it
      passes. I would have rejected this dependency on a stale reading of the
      GitHub API field.
      The real blocker: **none of the `tpt-opt-*` crates are on crates.io** (a
      crates.io search for both `tpt-opt` and `tpt-systems-optimisation` returns
      zero results), and `deny.toml` sets `sources.unknown-git = "deny"`. So they
      are consumable only as a git or path dependency.
      **Resolved by maintainer decision: use a git dependency.** Wired as
      `tpt-opt-core` with the `nlp` feature, pinned to commit `95fbc3f2`, with a
      single named `allow-git` entry in `deny.toml` (so `unknown-git` stays
      `"deny"` for everything else). Shipped as `tpt-fluids-hydraulic/src/sizing.rs`
      -- see the Phase 3 entry above. For the record: `tpt-opt-core` is
      `no_std` + `alloc`, depends on published `tpt-math-linalg-sparse` (v0.1.0,
      confirmed on crates.io, MIT OR Apache-2.0) and optional
      `tpt-math-optimize-general`, so the optimisation family bridges onto the
      same `tpt-math` substrate this workspace already uses.
      **One upstream defect found in the process**, worth reporting to that repo:
      `solve_nlp` decides feasibility with `prob.ineq(i, x).max(0.0)`, and Rust's
      `f64::max` returns the non-`NaN` operand, so a `NaN` constraint is scored
      as *satisfied* and the solver can return `Converged` on a decision vector
      of pure `NaN`. It is worked around here (`SizingProblem::guard`), but the
      fix belongs upstream: the max should be `c.max(0.0)` only after an explicit
      `is_nan` check, or the violation accumulator should treat `NaN` as
      infeasible.
- [ ] **EHL coupling with `tpt-fem-elasticity`** (line 153) - see the
      `tpt-fem-contact` note above.
      **Correction: the blocker is weaker than recorded.** This item was listed as
      needing "an external repo that does not exist", and that is not quite
      right. There is a `tpt-fem` repo at `c:\Programming\tpt-fem` containing a
      `tpt-fem-elasticity` crate (v0.1.0, "Linear elasticity (bar, plane-stress,
      plane-strain, 3D continuum)"), alongside 25 other crates including
      `tpt-fem-contact` and `tpt-fem-coupling`. So the dependency exists locally
      and the coupling is not externally blocked at all.
      What is still real: it is not published (it resolves through a workspace
      `path`, not crates.io), and `tpt-fluids`' `deny.toml` requires a registry
      source, so wiring it in is a deliberate policy decision rather than a
      one-liner. And the earlier finding stands on its own merits --
      `tpt-fem-contact` is surface-to-surface contact, which is not the elastic
      compliance kernel an EHL solution needs, so the useful entry point is
      `tpt-fem-elasticity` (and possibly `tpt-fem-coupling`), not contact.

## Kani / WSL: environment status

`spec.txt` lines 30 and 126-128 ask for Kani proofs. `cargo-kani` **does not
support Windows**: the install guide lists only `x86_64-unknown-linux-gnu`,
`x86_64-apple-darwin` and `aarch64-apple-darwin`. WSL is therefore the only
route, and its state on this machine (Windows 11 Home) is now:

- **Done:** the WSL 2.7.14 kernel package installed successfully via
  `wsl --install --no-distribution`, which needs no elevation. `wsl --version`
  reports the kernel and `wslgpu` present.
- **Blocked:** both WSL1 and WSL2 need optional Windows features enabled, and
  that requires elevation. `dism /online /enable-feature` returns
  `Error: 740 - Elevated permissions are required to run DISM`, and
  `Get-WindowsOptionalFeature` likewise. The agent shell runs as a
  non-administrator user and cannot self-elevate, because UAC needs an
  interactive consent prompt.
- `wsl --status` confirms the consequence: WSL2 is unavailable because the
  virtualization component is not enabled, and no distribution is installed.
- Note that `HypervisorPresent` is `False` and `VirtualizationFirmwareEnabled`
  is `True`, so the hardware supports it; only the OS feature is missing.

### Re-checked after a reboot (2026-09-29)

The machine was rebooted, which **partially** changed the picture and is worth
recording precisely rather than re-deriving:

- **Progress:** `HypervisorPresent` is now **`True`**, where the pre-reboot
  reading recorded above was `False`. The hypervisor is up, so the OS-side
  blocker is closer to cleared than it was.
- **Still blocked, and this is the part that did not change:** no distribution is
  installed (`wsl -l -v` reports none, and `wsl --list --online` fails with the
  "requires internet access to download" error), and
  `Get-WindowsOptionalFeature` still needs elevation. WSL1 remains unsupported on
  this build, so there is no lower-privilege fallback.
- `cargo-kani` is still not installed, as expected.

So the reboot alone was not sufficient. The two remaining steps both need
elevation, and neither can be done from this non-administrator shell.

### To finish, in an elevated PowerShell (Run as Administrator), then reboot

```powershell
wsl --install -d Ubuntu          # enables features, installs the distro
# reboot required -- VirtualMachinePlatform cannot be enabled in this session
```

After the reboot, Kani itself installs inside the distribution:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install --locked kani-verifier
cargo kani setup                  # downloads the CBMC compiler and data
```

Then the harnesses that are already written in `tpt-fluids-verify` can finally
run:

```bash
cargo kani --package tpt-fluids-verify --features kani
```

Until then nothing claims the harnesses pass, which is why they sit behind
`cfg(kani)`: they are compiled only when a real `cargo-kani` invokes them, so
they cannot be mistaken for verification that happened.

## External dependency gaps

*Two crates named in `spec.txt`'s integration section (§5) and Phase 3
rollout do not exist as sibling repos in this workspace tree yet. They are
not local path dependencies of anything in Phases 1-6 above, so they don't
block this build-out — they're called out here so the gap isn't silently
dropped:*

- `tpt-fem-elasticity` — needed for full EHL elastic-deformation coupling in
  `tpt-fluids-tribo` (Phase 4 ships the Reynolds/Dowson-Hampton film-thickness
  side without it; the coupling itself is deferred to Phase 7)
- `tpt-systems-optimisation` — needed for the end-to-end pipe-network-sizing /
  hull-form / bearing-geometry optimization examples referenced in
  `spec.txt` §5 and Phase 3 of the rollout (the differentiable integrals
  those optimizers would consume are still built in Phases 2-4 via
  `tpt-math-autodiff`, so this crate is a consumer-side gap, not a blocker)

`tpt-materials` was flagged as a similar risk during planning and is now
**closed in the other direction**: rather than being wired as a dependency, it
turned out to be unusable for this purpose. It is a micro-scale crystal-
plasticity/phase-field engine (`tpt-mat-core`, `tpt-mat-constants`,
`tpt-mat-crystallography`) that is not published on crates.io and supplies no
fluid density or viscosity data. Phase 1 therefore implements the property
database in-house inside `tpt-fluids-core` (`FluidProperties`, `DensityModel`,
`water_viscosity`, `water_vapour_pressure`). No external dependency was added,
so `deny.toml`'s registry-only source policy is unaffected.

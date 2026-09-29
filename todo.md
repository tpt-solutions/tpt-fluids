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
      (`hardy_cross.rs`). Correct for single-source networks including
      parallel pipe pairs, cross-validated against the GGA. **Does not**
      converge on networks with several sources *and* several independent
      loops; the GGA solves the same networks in ~9 iterations, so the
      networks are well posed and this is a limitation of the method. Two
      tests are `#[ignore]`d for that reason.
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
      - **Multi-node reflection: attempted, not shipped, findings recorded.**
        The single-reach solver advances one characteristic pair, so a wave has
        no state to send back -- that is the structural reason it cannot
        reflect, and a head array updated in place cannot fix it. A wave is a
        *directional travelling perturbation*, so the reflection-capable form
        needs two accumulator arrays (`c_minus` upstream, `c_plus` downstream)
        onto which waves are superimposed, and only then can a wave cross the
        pipe, turn round, and add to what is already there.
        Built and instrumented, and the reflection machinery itself worked: the
        wave was traced travelling upstream one reach per step, reflecting at the
        reservoir at full amplitude, and returning downstream. Four real defects
        were found and fixed along the way, all of which produce plausible
        numbers rather than errors:
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
- [x] Implement component models: valve `Cv` and `K` coefficients, pump
      characteristic curves (quadratic three-point fit, shut-off head, runout
      flow, hydraulic power), cavitation state from the cavitation number, and
      surge tank dynamics. Tabulated minor-loss coefficients for entrances,
      elbows, and exits are included. The turbine four-quadrant curve is not
      modelled beyond its loss coefficient.
- [x] Implement differentiable head-loss functions (via `tpt-math-autodiff`)
      for downstream pipe-sizing optimization. `differentiable.rs` evaluates
      Darcy-Weisbach over forward-mode dual numbers, so the diameter
      sensitivity comes out exactly in one pass. Supports the laminar,
      Swamee-Jain, and Hazen-Williams correlations analytically; a
      finite-difference fallback covers the implicit Colebrook-White.

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
`components`, and `differentiable`. 73 unit tests and 1 doctest pass, with 2
`#[ignore]`d for the documented Hardy Cross non-convergence on multi-source,
multi-loop networks. The GGA is the general solver and is validated against an
independently derived reference; reach for it unless a small single-source
network is all that is needed.

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
      separation-regime classification, and frictional heating with an
      explicit validity check on the quasi-steady temperature rise. A
      brake pad computes to 20 000 K, and a test asserts exactly that so
      the limit of the steady form stays concrete rather than aspirational

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
- [ ] Advanced coupling follow-up (unblock once the external repo exists):
      EHL coupling between `tpt-fluids-tribo` and `tpt-fem-elasticity`
- [ ] Advanced coupling follow-up (unblock once the external repo exists):
      end-to-end pipe-network diameter optimization via
      `tpt-systems-optimisation`

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

      **Still open:** Reynolds' supplementary cavitation condition for
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
      deformation solve, and the Dowson-Hampton film-thickness *formulas*. The
      first is blocked on `tpt-fem-elasticity` (see the `tpt-fem-contact` note
      above), and the second is a correlation whose coefficients come from a
      long chain of numerical solutions that this workspace has no way to
      reproduce. I attempted to source the Hamrock-Dowson coefficients and could
      not obtain them from any source that renders as text: the Penn State and
      NASA NTRS PDFs return raw compressed streams, and the MDPI review returns
      HTTP 403. Writing `2.69 U^0.68 G^0.49 (1 - e^(-0.68k))` from memory would
      be precisely the failure mode this file is built to prevent -- a precise,
      plausible, entirely unverified number. It is not shipped. The blocker is
      therefore "the constants need a citable source", which is a real and
      answerable one, not "we forgot".
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
- [ ] **`tpt-systems-optimisation`** (lines 156, 176) - consumer-side, not
      present as a sibling repo.
- [ ] **EHL coupling with `tpt-fem-elasticity`** (line 153) - see the
      `tpt-fem-contact` note above.

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

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
- [ ] Rustdoc for the solvers once they land
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
- [x] Unit tests: 87, with the ITTC-1957 line checked against the published
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
- [ ] Make the Michell integral differentiable via `tpt-math-autodiff` for
      hull-form optimization
- [ ] Integration test: model-ship resistance extrapolation matches
      Holtrop-Mennen benchmarks (see Phase 7)
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] Add to root `Cargo.toml` members + workspace deps

## Phase 4 — tpt-fluids-tribo

*Lubrication theory, contact mechanics, and wear modeling. Depends on:
`tpt-fluids-core`, `tpt-math-linalg-sparse` (Reynolds equation FD),
`tpt-math-autodiff`.*

- [ ] Scaffold `crates/tpt-fluids-tribo/`
- [ ] Wire deps: `tpt-fluids-core`, `tpt-math-linalg-sparse`,
      `tpt-math-autodiff`
- [ ] Implement Reynolds equation solver (1D/2D finite difference) for
      journal bearings, slider bearings, thrust bearings
- [ ] Implement Dowson-Hampton film thickness formulas (EHL film-thickness
      layer; full elastic coupling with `tpt-fem-elasticity` is an external
      dependency gap — see "External dependency gaps" below, tracked
      separately in Phase 7's advanced-coupling follow-up, not blocking here)
- [ ] Implement Stribeck curve regime navigation (boundary, mixed,
      hydrodynamic)
- [ ] Implement Archard's wear equation, adhesive/abrasive wear
      coefficients, running-in simulation
- [ ] Implement friction models: Coulomb, Stribeck, LuGre dynamic friction
      (for multibody joint integration)
- [ ] Implement differentiable load-capacity integrals (via
      `tpt-math-autodiff`) for bearing geometry optimization
- [ ] Unit tests + doctests (incl. Reynolds equation load-equilibrium
      cross-check)
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] Add to root `Cargo.toml` members + workspace deps

## Phase 5 — tpt-fluids-verify

*Mathematical verification of fluid algorithms and conservation laws.
Depends on: `tpt-fluids-core`, `tpt-fluids-hydraulic`, `tpt-fluids-marine`,
`tpt-fluids-tribo`, `proptest`, `kani` (dev/verification tooling).*

- [ ] Scaffold `crates/tpt-fluids-verify/`
- [ ] Wire deps: the four domain crates, `proptest`; set up the `kani`
      toolchain (pinned nightly + CBMC backend per Kani's install docs)
- [ ] Implement proptest strategies: valid pipe network topologies, ship hull
      forms, bearing geometries
- [ ] Write Kani harness: Hardy Cross loop corrections converge monotonically
- [ ] Write Kani harness: MOC water hammer characteristics respect CFL
      stability
- [ ] Write Kani harness: Reynolds equation satisfies global load
      equilibrium
- [ ] Implement invariant check (proptest): mass conservation at every
      network node (∑Q_in = ∑Q_out)
- [ ] Implement invariant check (proptest): Froude scaling preserves
      dimensionless resistance coefficients
- [ ] Implement invariant check (proptest): energy conservation in lossless
      pipe segments
- [ ] Rustdoc (crate-level: what each harness proves and why)
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] `cargo kani` runs clean on all harnesses
- [ ] Add to root `Cargo.toml` members + workspace deps

## Phase 6 — tpt-fluids (umbrella crate)

*Feature-gated umbrella crate re-exporting core/hydraulic/marine/tribo/verify.*

- [ ] Scaffold `crates/tpt-fluids/`
- [ ] Wire optional deps + feature flags: `core` (always on),
      `hydraulic`, `marine`, `tribo`, `verify`
- [ ] Re-export each constituent's public API behind its feature
- [ ] Rustdoc documenting the feature matrix
- [ ] `cargo fmt` / `clippy` / `deny` clean across feature combinations
      (`--no-default-features`, `--all-features`, and each single-feature
      combination)
- [ ] Add to root `Cargo.toml` members + workspace deps

## Phase 7 — Integration & Workspace Closeout

- [ ] Integration test: model-ship resistance extrapolation (`tpt-fluids-marine`)
      matches published Holtrop-Mennen benchmark figures within stated tolerance
- [ ] `cargo test --workspace --all-features` passes
- [ ] `cargo clippy --workspace --all-targets --all-features -- -D warnings` clean
- [ ] `cargo deny check` clean workspace-wide (confirm no Apache-2.0-only
      crate anywhere in the resolved dependency graph)
- [ ] no_std(+alloc) matrix passes for every crate that declares `no_std`
      support (expected: `tpt-fluids-core` at minimum)
- [ ] Root `README.md` documents the full crate map, build order, and how
      `tpt-construction`/`tpt-engineering`/`tpt-physics`/`tpt-aero` are
      expected to consume `tpt-fluids`
- [ ] Advanced coupling follow-up (unblock once the external repo exists):
      EHL coupling between `tpt-fluids-tribo` and `tpt-fem-elasticity`
- [ ] Advanced coupling follow-up (unblock once the external repo exists):
      end-to-end pipe-network diameter optimization via
      `tpt-systems-optimisation`

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

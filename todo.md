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

- [ ] Scaffold `crates/tpt-fluids-hydraulic/`
- [ ] Wire deps: `tpt-fluids-core`, `tpt-math-graph`, `tpt-math-linalg`,
      `tpt-math-linalg-sparse`, `tpt-math-autodiff`
- [ ] Implement friction models: Darcy-Weisbach, Hazen-Williams,
      Colebrook-White iterative solver, minor loss coefficients
- [ ] Implement network topology on `tpt-math-graph`: directed graph
      representation, node/edge/loop extraction
- [ ] Implement Hardy Cross method (loop-based) network solver
- [ ] Implement Global Gradient Algorithm (node-based Newton-Raphson) network
      solver (via `tpt-math-linalg`/`tpt-math-linalg-sparse`)
- [ ] Implement transient analysis: water hammer via Method of
      Characteristics (MOC), Joukowsky equation, column separation modeling
- [ ] Implement component models: pump/turbine four-quadrant characteristic
      curves, valve Cv models, cavitation inception, surge tank dynamics
- [ ] Implement differentiable head-loss functions (via `tpt-math-autodiff`)
      for downstream pipe-sizing optimization
- [ ] Unit tests + doctests (incl. a hand-verified small network for Hardy
      Cross / GGA cross-check)
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] Add to root `Cargo.toml` members + workspace deps

## Phase 3 — tpt-fluids-marine

*Naval architecture and marine hydrodynamics for ships and underwater
vehicles. Depends on: `tpt-fluids-core`, `tpt-math-linalg`/
`tpt-math-linalg-fixed`, `tpt-math-signal-fft` (wave spectra),
`tpt-math-autodiff` (Michell integral differentiability).*

- [ ] Scaffold `crates/tpt-fluids-marine/`
- [ ] Wire deps: `tpt-fluids-core`, `tpt-math-linalg`,
      `tpt-math-linalg-fixed`, `tpt-math-signal-fft`, `tpt-math-autodiff`
- [ ] Implement frictional resistance: ITTC-1957 model-ship correlation
      line, Granville method
- [ ] Implement residuary resistance: Michell's thin-ship integral,
      Holtrop-Mennen empirical method
- [ ] Implement air drag, appendage drag, correlation allowance
- [ ] Implement Froude scaling: model-to-ship extrapolation with form factor
      and roughness allowance
- [ ] Implement seakeeping: linear wave theory, Response Amplitude Operators
      (RAOs) for 6-DOF ship motions, significant wave height statistics
- [ ] Implement MMG model for 3-DOF horizontal-plane maneuvering (surge,
      sway, yaw), hydrodynamic derivatives
- [ ] Implement propulsion: wake fraction, thrust deduction,
      propeller-hull interaction, open-water diagrams
- [ ] Implement Froude-Krylov excitation for underwater vehicles / floating
      structures
- [ ] Make the Michell integral differentiable via `tpt-math-autodiff` for
      hull-form optimization
- [ ] Unit tests + doctests
- [ ] Integration test: model-ship resistance extrapolation matches
      Holtrop-Mennen benchmarks (see Phase 7)
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] Add to root `Cargo.toml` members + workspace deps

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

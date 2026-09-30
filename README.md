# tpt-fluids

**Applied fluid mechanics for the TPT Solutions engineering stack.** Pure Rust,
no C/C++ FFI, no external solver binaries, no proprietary libraries.

A feature-gated workspace covering three fluid domains:

- **1D hydraulic networks** - Darcy-Weisbach friction with Colebrook-White,
  Swamee-Jain, and Hazen-Williams; network topology with a true fundamental
  cycle basis; the Hardy Cross and Global Gradient Algorithm solvers;
  water-hammer transients; and a library of component models.
- **Marine hydrodynamics** - ITTC-1957 and ITTC-57 friction, Granville's
  roughness and appendage extension, the ITTC-1957 model-ship line, Michell's
  thin-ship residuary integral, Holtrop-Mennen resistance, Froude scaling,
  seakeeping response amplitude operators, MMG manoeuvring, propulsion, and
  Froude-Krylov wave excitation. The Michell integral and the total resistance
  chain are also available over forward-mode dual numbers, so a gradient-based
  hull-form optimiser gets `dR/dFr` exactly.
- **Tribology** - Hertzian contact, Petroff and Stribeck lubrication, the
  hydrodynamic journal bearing, Archard's wear law, the lambda ratio, and
  frictional heating.

Every solver and correlation is implemented from scratch on top of the
[`tpt-math`](https://github.com/tpt-solutions/tpt-math) substrate. There is no
EPANET, WAMIT, or tribology-solver FFI anywhere in the graph.

## The crate map

| Crate | What it holds | `no_std` |
|-------|---------------|----------|
| `tpt-fluids-core` | Unit-safe quantities, dimensionless numbers, equations of state, fluid properties, viscosity, surface tension, constants, `libm` shims | **yes**, `thumbv6m-none-eabi` verified |
| `tpt-fluids-hydraulic` | Network topology, friction factors, Hardy Cross, GGA, water hammer, components, differentiable head loss | no |
| `tpt-fluids-marine` | Ship resistance, Froude scaling, seakeeping, propulsion, MMG manoeuvring, wave excitation | no |
| `tpt-fluids-tribo` | Hertzian contact, lubrication, wear, frictional heating, EHL film | no |
| `tpt-fluids-verify` | Proptest invariants and Kani proof harnesses | no |
| `tpt-fluids` | Umbrella: feature-gated re-exports and a `prelude` | partially, `core` only |

**Build order** is the dependency order: `core` first, then the three domain
crates in parallel, then `verify` (which needs all three), then the umbrella.
Cargo resolves this from the manifests; nothing has to be built in a
particular order by hand.

The three domain crates declare a `std` feature and are not `no_std`. Only
`core` is, and that is a deliberate split: the correlations are pure
arithmetic that belongs on a microcontroller, while the solvers allocate
network topologies and are not going to.

## Design philosophy

- **Unit-safe.** Phantom-typed quantities make `Pa + m` a compile error, and
  Reynolds, Froude, Weber, and cavitation numbers are distinct types that
  cannot be confused with each other or with dimensioned quantities.
- **Dimensional checks before reference values.** The most dangerous bugs in
  this domain are not wrong correlations, they are correlations that are
  internally consistent and dimensionally wrong by a factor of 1000. Several
  real ones were caught this way, including Petroff's friction law missing a
  length, the ITTC-1957 hull-form figure being 25x the skin-friction
  coefficient it is often mistaken for, and the quasi-steady temperature rise
  returning kelvin per metre because it too was missing a conduction length.
  Each is now a test. The temperature one is the clearest argument for the
  practice: it produced a plausible-looking 20 000 K for a brake pad, a number
  no material survives, and a reference-value test would have *endorsed* it.
- **Documented limits.** Where a model is an approximation, or is only valid
  in a regime, the module says so at the point of use. The transient flash
  temperature is now solved rather than merely flagged, and what remains
  unmodelled — the layered pad-on-disc composite — is named.
- **Documentation is a build gate.** `cargo doc` runs in CI under
  `RUSTDOCFLAGS: -D warnings`, so an unresolved intra-doc link fails the build.
  A dead cross-reference is otherwise invisible: rustdoc drops the link and
  emits the text as plain, with no compile error and no test failure. That gate
  caught four on its first run, two of them written while fixing unrelated
  bugs — the exact failure a test suite cannot see.
- **Declared dimensions are enforced at compile time.** Every product and
  quotient in `core` asserts its dimension in a `const` context, so a wrong
  unit is a build failure rather than a plausible wrong answer. The check found
  three relations the previous macro had emitted without verifying any of them,
  including `Length * Time = Velocity`.

## Usage

```rust
use tpt_fluids::prelude::*;

let speed = Velocity::new(12.5 * 0.514_444);
let friction = tpt_fluids_marine::resistance::friction_resistance(
    Length::new(300.0),
    Length::new(45.0),
    Length::new(14.0),
    speed,
    Density::new(1025.0),
    1.05e-6,
);
```

Features on the umbrella: `core` (default), `hydraulic`, `marine`, `tribo`,
`verify`, `std`. `--no-default-features --features core` builds `no_std`.

## Testing

548 tests across the workspace, plus 8 doctests. **No test is `#[ignore]`d.**

| Suite | Count | What it checks |
|-------|-------|----------------|
| `core` | 44 | Quantity algebra, **compile-time dimension checks**, EOS, viscosity |
| `hydraulic` | 142 | Correlations against published values; solver convergence; network sizing; MOC wave propagation, reflection, and friction |
| `marine` | 140 | Hull form, resistance, seakeeping, manoeuvring, propulsion, Froude scaling |
| `tribo` | 155 | Contacts, lubrication, wear, friction, flash temperature, Reynolds film solve, EHL film thickness |
| `verify` | 41 | Property sweeps over every correlation, including the MOC solver and EHL |
| `umbrella` | 9 | Cross-crate consistency through the re-exports |
| `marine` benchmarks | 8 | Two independent methods agreeing about the same ship |

The two hydraulic tests that used to be `#[ignore]`d are now enabled. Both
covered Hardy Cross on a multi-source, multi-loop network, which the solver
previously could not close; it now solves that case and the tests assert it
directly. One of the fixes that made them pass was a real bug rather than a
missing feature: the module's `continuity_error` metric counted fixed-head
reservoirs as continuity failures, so a converged two-reservoir solution
reported an "error" equal to the network's own net demand. The GGA has always
excluded them, which is why the two solvers disagreed about what the number
meant.

Kani proof harnesses are written and gated behind the `kani` feature and
`cfg(kani)`. `cargo-kani` is not part of a normal toolchain, so they are not
claimed to pass; they are behind a cfg so they are not silently passing
either.

## Build and check

```text
cargo build --workspace
cargo test  --workspace --all-features
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo fmt   --all -- --check
cargo deny  check
cargo build -p tpt-fluids-core --no-default-features --features alloc \
    --target thumbv6m-none-eabi
```

## Known limitations

These are stated rather than hidden, and tracked in `todo.md`:

- **Water hammer** has a multi-node MOC network solver (`MocNetwork`): each
  branch is discretised into reaches on a Courant-common time step, a wavefront
  travels one reach per step, and the node boundaries are a fixed-head reservoir,
  a prescribed-flow valve, a reflecting dead end, and a continuity junction. A
  dead end reflects at full amplitude -- the head disturbance doubles on arrival,
  which is the form needed to quote a Joukowsky head at a closed end -- and
  `solve_with_heads` lets a reservoir whose level moves **launch** a wave rather
  than only absorb one. **Unsteady friction** is available opt-in through
  `with_friction_lag`, which lags the friction force and damps the surge the
  quasi-steady law over-predicts; it is off by default and provably cannot
  perturb a steady result. `with_derived_friction_lag` derives the time constant
  per branch from each pipe's own `2 L / (g |V|)`.
- **Tribology `contact`** carried a silent zero: `contact_radius` and
  `approach` both answered `0.0` for an infinite reduced radius, reporting a
  loaded flat joint as having no contact patch and no deformation at all. The
  point-contact formula simply does not apply there — it is a cylinder problem —
  so the answer is now `NaN`, and `flat_approach` is the determinable route.
  Both branches were untested, and the existing test asserted the wrong answer.
- **Seakeeping RAOs** are damped-oscillator approximations, not
  Green-function hull-integral solutions. Their peak, phase and monotonicity
  are verified against the closed-form oscillator, but they are not a hull
  integral.
- **The flash temperature** treats the contact as one homogeneous semi-infinite
  body. The transient (Blok-Wilde) form is implemented, so the quasi-steady
  absurdity is resolved, but a layered pad-on-disc composite is not solved.
- **The coupled EHL problem is not solved**, but the **film thickness is**.
  `tpt-fluids-tribo::ehl` delivers the Dowson-Higginson line contact and the
  Hamrock-Dowson point contact, the minimum (side-lobe) film, and the resulting
  lambda ratio — the quantity every lubricant-selection decision turns on, and
  which needs no elastic solver because these are closed-form correlations.
  What is missing is the *coupled* solution: the pressure distribution and
  sub-surface stress from solving Reynolds together with the elastic
  deformation of both bodies, which does need the external
  `tpt-fem-elasticity` crate.
- **Hardy Cross** closes every loop, but it cannot choose between the several
  flow fields its own equations admit, so on a multi-source network it can
  converge to a loop-consistent answer that is not the physical one. The GGA is
  the solver to reach for; see `hardy_cross` for the details.

Two items previously listed here no longer apply:

- ~~**Turbines** are a loss coefficient, not a four-quadrant curve.~~ The
  four-quadrant model exists: `components::TurbineCurve` covers generating,
  windmilling and pumping, with `runout_flow` and the extracted hydraulic
  power.
- ~~**Colebrook** falls back to finite differences for the derivative.~~ The
  derivative over dual numbers is now **analytic and exact** for every
  correlation including Colebrook-White, by running the iteration itself over
  forward-mode duals rather than approximating its derivative.
- ~~**EHL and LuGre friction** are not implemented.~~ LuGre dynamic friction is
  implemented in `tpt-fluids-tribo::friction`, with Coulomb and the Stribeck
  baseline. EHL **film thickness** is now implemented too, in
  `tpt-fluids-tribo::ehl`: Dowson-Higginson for line contacts, Hamrock-Dowson for
  point contacts, the minimum-film and lambda-ratio consequences, and property
  tests pinning every exponent. Only the *coupled* pressure solution is
  outstanding, and that one is listed above.
- ~~**The flash temperature** is the quasi-steady form only.~~ The transient
  Blok-Wilde solution is implemented, and it exposed a dimensional bug in the
  quasi-steady one: `mu F v / (k A)` is kelvin *per metre*, not kelvin, because
  conduction needs a length. Both forms now take a conduction length, and the
  documented "20 000 K brake pad" was a consequence of that missing factor —
  the same pad computes to about 1 600 K, which is a real disc temperature.

## Licence

MIT OR Apache-2.0.

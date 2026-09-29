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
| `tpt-fluids-tribo` | Hertzian contact, lubrication, wear, frictional heating | no |
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
  length, and the ITTC-1957 hull-form figure being 25x the skin-friction
  coefficient it is often mistaken for. Each is now a test.
- **Documented limits.** Where a model is an approximation, or is only valid
  in a regime, the module says so at the point of use. The quasi-steady flash
  temperature computes to 20 000 K for a brake pad, and there is a test
  asserting exactly that number so nobody mistakes it for a result.

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

446 tests across the workspace, plus 8 doctests.

| Suite | Count | What it checks |
|-------|-------|----------------|
| `core` | 42 | Quantity algebra, dimensional consistency, EOS, viscosity |
| `hydraulic` | 94 + 2 ignored | Correlations against published values; solver convergence |
| `marine` | 139 | Every correlation against a reference value, plus limit cases |
| `tribo` | 124 | Hertz against its closed forms; Stribeck, Archard, lambda |
| `verify` | 31 | Proptest invariants: dimensional invariance, monotonicity, exact scalings |
| `umbrella` | 8 | Cross-crate consistency through the re-exports |
| `marine` benchmarks | 8 | Two independent methods agreeing about the same ship |

The two ignored hydraulic tests are the Hardy Cross multi-loop case that is
documented in `todo.md` as not converging; they are left in place rather than
deleted so the limitation stays visible.

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

- **Water hammer** handles one reservoir-fed reach, not multi-node wave
  reflection.
- **Turbines** are a loss coefficient, not a four-quadrant curve.
- **Hardy Cross** does not converge on multi-source, multi-loop networks. The
  GGA solver covers that case and is the recommended one.
- **Colebrook** falls back to finite differences for the derivative rather
  than an analytic gradient.
- **Seakeeping RAOs** are damped-oscillator approximations, not
  Green-function hull-integral solutions.
- **The flash temperature** is the quasi-steady form only; a transient
  solution is not attempted.
- **EHL and LuGre friction** are not implemented, despite earlier drafts of
  this README claiming otherwise.

## Licence

MIT OR Apache-2.0.

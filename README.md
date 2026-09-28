# tpt-fluids

**Applied fluid mechanics for the TPT Solutions engineering stack.** Pure Rust,
no C/C++ FFI, no external solver binaries.

`tpt-fluids` consolidates three specialized fluid domains into one
feature-gated workspace:

- **1D hydraulic networks** — Darcy-Weisbach / Colebrook-White friction, Hardy
  Cross and Global Gradient Algorithm network solvers, water-hammer transients.
- **Marine hydrodynamics** — ITTC-1957 / Granville friction, Michell's thin-ship
  integral, Holtrop-Mennen residuary resistance, Froude scaling, seakeeping
  RAOs, MMG maneuvering, propulsion.
- **Tribology** — Reynolds-equation lubrication, Dowson-Hampton EHL film
  thickness, Stribeck regimes, Archard wear, LuGre friction.

Every solver is implemented **from scratch** on top of the
[`tpt-math`](https://github.com/tpt-solutions/tpt-math) substrate. There is no
EPANET, WAMIT, or tribology-solver FFI anywhere in the graph.

## Design philosophy

- **Unit-safe.** Phantom-typed quantities make `Pa + m` a compile error, and
  Reynolds/Froude/Weber/Mach/cavitation numbers are distinct types that cannot
  be confused with each other or with dimensioned quantities.
- **`no_std` + `alloc` at the core.** `tpt-fluids-core` builds for bare-metal
  targets (`thumbv6m-none-eabi`).
- **Differentiable.** Head-loss integrals and load-capacity integrals are
  expressed so that `tpt-math-autodiff` can propagate gradients, enabling
  gradient-based sizing and hull/bearing optimization.
- **Verified.** `tpt-fluids-verify` carries Kani harnesses (CFL stability,
  Hardy Cross monotonic convergence, Reynolds load equilibrium) and proptest
  invariants (node mass conservation, Froude invariance, lossless energy
  conservation).

## License posture

Every crate is `MIT OR Apache-2.0`. **No Apache-2.0-only dependency is permitted
anywhere in the graph** — the allow-list is enforced by `deny.toml` and CI.

## Crate map

| Crate | Purpose | `no_std` |
|-------|---------|----------|
| `tpt-fluids-core` | Unit-safe types, equations of state, viscosity models, fluid property database | yes (alloc) |
| `tpt-fluids-hydraulic` | 1D pipe networks, friction, Hardy Cross, GGA, water hammer | no |
| `tpt-fluids-marine` | Ship resistance, Froude scaling, seakeeping, MMG maneuvering, propulsion | no |
| `tpt-fluids-tribo` | Reynolds equation, EHL, Stribeck, Archard wear, friction | no |
| `tpt-fluids-verify` | Kani harnesses + proptest invariant strategies | no |
| `tpt-fluids` | Feature-gated umbrella crate re-exporting the above | no |

## Consumers

`tpt-fluids` is the applied fluid-mechanics layer beneath `tpt-construction` and
`tpt-engineering` (HVAC, plumbing, municipal water networks), and it hands
hydrodynamic and friction models to `tpt-multibody-dynamics`. `tpt-physics` (CFD)
is reserved for full 3D Navier-Stokes where 1D and semi-empirical specialization
stops being accurate.

## Building

```sh
cargo build --workspace
cargo test  --workspace --all-features
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo deny check
```

`spec.txt` holds the full design specification; `todo.md` tracks the phased
build-out.

Licensed under MIT OR Apache-2.0. Copyright (c) 2026 TPT Solutions.

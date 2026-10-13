# tpt-fluids-multiphase-bev

Gas-liquid physics for carbonated beverages, part of the
[tpt-fluids](https://github.com/tpt-solutions/tpt-fluids) workspace.

It answers two formulation questions for a drink containing xylitol, dissolved
minerals, microencapsulated particles and pH buffers:

1. **How much CO2 does it hold, and how is it shared?** Henry's law with
   Sechenov salting-out, carbonic-acid and buffer speciation, and the
   bicarbonate reservoir that releases CO2 slowly.
2. **Will it foam over or go flat?** Heterogeneous classical nucleation with
   contact angles, pre-existing gas pockets, diffusion-limited growth, and a
   time-integrated release model that gives a flat / stable / foam-over verdict.

## Status

The physics is standard and its reference values are tested. The formulation
inputs are **not** measured: the xylitol Sechenov constant
(`Solute::xylitol_estimate`, tagged `Provenance::Estimate`), particle contact
angles, cavity size distributions and the foam-over thresholds
(`foaming::Thresholds::default`, marked uncalibrated) are placeholders. Use the
model to rank formulations; calibrate against bench data before reading absolute
foam-over or shelf-life numbers.

## Licence

MIT OR Apache-2.0.

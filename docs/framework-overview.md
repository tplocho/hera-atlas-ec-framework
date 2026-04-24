# Framework Overview

## High-level split
The repository is organized around a two-layer idea:

1. **HERA** — generation of resonance-related quantities from tokamak geometry, plasma-profile assumptions, wave propagation choices, and resonance physics.
2. **ATLAS** — post-processing, target evaluation, ranking, and visualization of actuator configurations from generated datasets.

## Practical interpretation
- HERA is the physics-generation layer.
- ATLAS is the decision / ranking / visualization layer.
- `ATLAS_int.m` is the interactive entry layer.

## Public-safe release principle
This repository deliberately separates the reusable code from manuscript-specific packaging. That keeps the public release truthful and avoids exposing unstable submission material prematurely.

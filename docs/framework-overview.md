# Framework Overview

## High-level structure
The repository is organized around a two-layer workflow for electron-cyclotron resonance analysis and actuator evaluation.

### HERA
`HERA.m` is the resonance-physics generation layer. It defines tokamak and plasma-profile inputs, evaluates propagation and resonance-related quantities, and writes tabulated outputs for downstream analysis.

### ATLAS
`ATLAS.m` is the actuator-analysis layer. It operates on generated datasets to evaluate targeting scenarios, compare actuator settings, and produce ranking and visualization outputs.

### ATLAS_int
`ATLAS_int.m` is the interactive configuration layer for ATLAS.

## Practical flow
A typical use pattern is:
1. run or prepare datasets associated with the desired propagation / mode assumptions
2. configure the ATLAS analysis interactively or through a scripted configuration structure
3. run ATLAS for ranking, landscape mapping, or heating / targeting analysis

## Scope of the repository
The repository is centered on the reusable software components of the framework rather than manuscript packaging or publication-layout assets.

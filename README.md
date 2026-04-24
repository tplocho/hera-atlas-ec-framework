# HERA–ATLAS EC Framework

MATLAB research code for electron-cyclotron actuator evaluation and resonance mapping in ITER-class tokamak scenarios.

## Overview
The repository contains a two-layer workflow for studying resonance structure and actuator relevance under prescribed tokamak geometry, plasma-profile, and propagation assumptions:

- **HERA** generates resonance-related quantities from the underlying physics model.
- **ATLAS** post-processes generated datasets to evaluate, rank, and visualize actuator-target relationships.
- **ATLAS_int** provides an interactive configuration path for ATLAS.

## Repository contents
- `src/HERA.m` — resonance-physics generation layer
- `src/ATLAS.m` — actuator-target analysis and visualization layer
- `src/ATLAS_int.m` — interactive configuration layer
- `docs/framework-overview.md` — framework structure and workflow summary
- `docs/repository-scope.md` — repository boundaries and included material
- `docs/usage-and-data.md` — notes on execution flow and generated datasets

## Typical workflow
1. Generate or prepare resonance datasets with `HERA.m`.
2. Configure an analysis through `ATLAS_int.m` or a scripted configuration structure.
3. Run `ATLAS.m` to analyze, rank, and visualize actuator-target configurations.

## Quick start
```matlab
cfg = ATLAS_int();
results = ATLAS(cfg);
```

## Repository scope
This repository contains the code and supporting documentation for the HERA–ATLAS framework. Manuscript drafts, internal review material, and selected publication assets are maintained separately from the repository.

## Affiliation
Theotokis Plochoros  
National Technical University of Athens (NTUA)

## Contact
`tplocho@gmail.com`

## Citation
Citation metadata is provided in `CITATION.cff`.

## License
MIT License.

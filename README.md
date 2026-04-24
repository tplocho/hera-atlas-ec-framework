# HERA–ATLAS EC Framework

MATLAB research code for a physics-informed framework for electron-cyclotron actuator evaluation and resonance mapping in ITER-class tokamak scenarios.

## Status
This repository is a **public-safe code cut** of an active research project. It is intended to expose the software layer cleanly while manuscript drafts, submission-specific figures, and internal review artifacts remain outside the public repository.

## What is included
- `src/HERA.m` — resonance-physics generation layer
- `src/ATLAS.m` — actuator-target post-processing and ranking layer
- `src/ATLAS_int.m` — interactive configuration layer
- lightweight documentation describing scope and current release posture

## What is intentionally excluded in this cut
- manuscript drafts
- internal audit PDFs
- submission-sensitive figures
- unstable publication packaging assets
- any materials not yet cleared for stable public release

## Purpose
The framework sits between overly simple geometric resonance criteria and much heavier ray-/beam-tracing workflows. Its role is to provide a computationally lighter, interpretable, control-oriented layer for exploring EC actuator relevance across tokamak geometry and operating conditions.

## Current repository posture
This public-safe cut is designed for:
- technical transparency
- portfolio presentation
- early code visibility
- future reproducibility hardening

It is **not yet** a finalized publication companion package.

## Quick start
A minimal entry sequence is:

```matlab
cfg = ATLAS_int();
results = ATLAS(cfg);
```

## Current limitations
Before treating this repository as a mature public research release, the following should be tightened:
- documented end-to-end example workflow
- clearer dataset-generation story for ATLAS inputs
- post-submission decision on which figures/results can be public

## Suggested starting points
Begin by reading the code headers in:
- `src/HERA.m`
- `src/ATLAS.m`
- `src/ATLAS_int.m`

Then see:
- `docs/framework-overview.md`
- `docs/release-gate.md`
- `docs/reproducibility-status.md`

## Affiliation
Theotokis Plochoros  
National Technical University of Athens (NTUA)

## Contact
`tplocho@gmail.com`

## License
MIT License.

## Citation
A machine-readable citation draft is provided in `CITATION.cff`. Once the project's publication state stabilizes, the preferred citation can be updated to the final paper and/or archived software release.

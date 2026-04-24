# Usage and Data Notes

## Core files
- `src/HERA.m`
- `src/ATLAS.m`
- `src/ATLAS_int.m`

## Execution pattern
`HERA.m` is used to generate resonance-related datasets under selected model and scan assumptions.
`ATLAS.m` then operates on those datasets for actuator evaluation, ranking, and visualization.

## Interactive path
A minimal interactive entry sequence is:

```matlab
cfg = ATLAS_int();
results = ATLAS(cfg);
```

## Dataset expectations
ATLAS is designed to operate on generated CSV datasets associated with the chosen propagation mode or analysis path. The file naming conventions and expected inputs are reflected in the source headers and interactive prompts.

## Working style
The framework can be used either:
- interactively through `ATLAS_int.m`, or
- programmatically by constructing a configuration structure and passing it to `ATLAS.m`

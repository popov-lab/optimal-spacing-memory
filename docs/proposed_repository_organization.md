# Proposed repository organization

Status: proposal only.  No repository-wide moves or renames have been made.

## Design goals

The repository should support additional learning/forgetting models, multiple
versions of source-of-activation-confusion models, further paper/data
extractions, and both deterministic and Monte Carlo studies without tying the
project to the layout of any one paper release.

The key boundaries should be:

- immutable third-party material versus project-authored code;
- raw/imported data versus deterministic derivatives versus study results;
- general model implementations versus paper-specific compatibility adapters;
- model definitions versus experimental schedules and fitting procedures;
- reference results versus new scientific results.

## Proposed target tree

```text
optimal-spacing-memory/
├── pyproject.toml
├── README.md
├── src/optimal_spacing_memory/
│   ├── models/
│   │   ├── base.py
│   │   ├── anderson_2023/
│   │   ├── generic_learning_forgetting/
│   │   └── source_activation_confusion/
│   ├── data/
│   │   ├── schemas.py
│   │   ├── loaders.py
│   │   └── validation.py
│   ├── simulation/
│   │   ├── engines.py
│   │   ├── random_streams.py
│   │   └── conditioning.py
│   ├── fitting/
│   │   ├── objectives.py
│   │   └── optimizers.py
│   └── metrics.py
├── external/
│   └── anderson_2023/
│       ├── README.md
│       ├── SHA256SUMS
│       └── matlab/
├── data/
│   ├── interim/
│   ├── derived/
│   │   └── anderson_2023/
│   └── schemas/
├── references/
│   ├── papers/
│   ├── model_notes/
│   └── bibliography/
├── studies/
│   ├── replications/
│   │   └── anderson_2023/
│   │       ├── configs/
│   │       ├── scripts/
│   │       └── reports/
│   └── comparisons/
├── results/
│   ├── manifests/
│   ├── predictions/
│   ├── fits/
│   └── figures/
├── docs/
│   ├── models/
│   ├── data/
│   ├── replications/
│   └── decisions/
├── tests/
│   ├── unit/
│   ├── integration/
│   └── regression/
└── scripts/
    └── one-off migration or release utilities only
```

## Responsibilities

### `src/optimal_spacing_memory/models`

Every model should implement a small common interface over an explicit
occurrence history, such as activation/odds/probability prediction plus
parameter validation.  Paper-specific prototype-history rules and known code
quirks should live in a named compatibility module, not in the common model.

Suggested semantic modes are explicit enums or configuration values, never
hidden booleans:

- `paper_intended`
- `release_compat`
- `corrected_common_engine`

The generic learning/forgetting model and future source-of-activation-confusion
versions would each receive their own subpackage and notation/reference page.

### `external`

Byte-for-byte third-party releases, accompanied by provenance, checksums,
access dates, and rights notes.  Project code must never import and then mutate
files in place.  Any adapter belongs under `src`; any derivative belongs under
`data/derived`.

### `data`

- `interim`: expensive, reproducible intermediate representations that are not
  stable public interfaces;
- `derived`: schema-versioned, analysis-ready tables generated from imported
  sources;
- `schemas`: machine-readable column definitions and constraints shared by
  loaders and validators.

Generated tables should carry a manifest with source hashes and pipeline
version.  Data loaders should select columns and stream/chunk large tables;
model code should not parse MATLAB or know CSV filenames.

### `references`

Human-facing mathematical transcriptions, paper PDFs where redistribution is
appropriate, bibliographic metadata, and notes about notation.  This is
separate from executable implementations so a model's source and derivation
remain reviewable even when code evolves.

### `studies` and `results`

A study is an executable scientific specification: inputs, model mode,
parameter bounds, seeds, optimizer/simulator settings, and acceptance tests.
Results are generated artifacts keyed by study/config hash.  Reference-cache
replays and new comparisons must be labeled separately.

Small canonical regression fixtures can be committed.  Large predictions,
optimizer traces, and Monte Carlo shards should use the selected large-data
strategy and should never be ordinary Git blobs.

### `tests`

- `unit`: hand-computed equations and validation rules;
- `integration`: MAT → CSV → loader → model paths;
- `regression`: small golden arrays and published summary anchors.

Stochastic tests should assert seeded summaries and invariants rather than
fragile full-array identity unless the random engine itself is the test.

## Practices to establish early

1. One installable Python package and one declared environment, with versions
   locked for replication runs.
2. Type-checked parameter objects with keyword-only construction and explicit
   conversion from released vector order.
3. Stable model, dataset, study, and condition IDs.
4. Configuration files for every fit or simulation; no scientific constants
   embedded only in notebooks.
5. Deterministic root seeds, child-stream derivation, and recorded shard IDs.
6. Input/output manifests with hashes and schema versions.
7. Automated fast tests on every change; expensive replication jobs as
   separately invoked workflows.
8. Decision records for notation changes, compatibility modes, data-policy
   choices, and deviations from paper code.
9. Notebooks only for exploration and exposition.  Any result used in a report
   must be reproducible through package code or a study runner.
10. Separate scientific correction from faithful reproduction: preserve both
    outputs when they answer different questions.

## Decisions needed before applying the reorganization

1. **Package/API name.** `optimal_spacing_memory` is descriptive, but a broader
   name may be preferable if the repository will become a general memory-model
   comparison framework.
2. **Large-data policy.** Choose Git LFS, DVC/object storage, or generated-local
   data with only manifests committed.  The current author release and full
   CSV derivatives are too large for ordinary Git history.
3. **Environment policy.** Choose a lockfile/workflow that also supports any R
   code that remains scientifically relevant.
4. **Result-retention policy.** Decide which fitted predictions and Monte Carlo
   summaries are canonical enough to version and which are disposable build
   products.
5. **Reference-material policy.** Decide which papers/code releases may be
   stored, linked, or redistributed, and record permission/provenance per item.

## Safe migration sequence after approval

1. Record the approved decisions in `docs/decisions`.
2. Add the package skeleton and compatibility-preserving imports; do not move
   scientific code yet.
3. Move the independent Anderson implementations with regression tests.
4. Move pipeline code behind stable loaders and schemas.
5. Introduce study configuration/runner infrastructure.
6. Migrate existing figures and older scripts, retaining provenance in an
   archive where they are not worth modernizing.
7. Add the user's new model families one at a time through the common API.
8. Only then remove compatibility shims and obsolete top-level paths.

This sequence keeps every step reviewable and reversible while the scientific
scope grows.

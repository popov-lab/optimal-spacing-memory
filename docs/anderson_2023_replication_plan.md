# Anderson et al. (2023) replication plan

## Objective

The replication should answer three different questions without conflating
them:

1. **Artifact replay:** Can the released MATLAB inputs and cached predictions
   reproduce the paper's reported summaries?
2. **Implementation replay:** Can an independently written implementation
   reproduce the released code's predictions, including documented quirks?
3. **Model test:** How do the paper-intended equations behave under a coherent,
   reproducible simulator and robust fitting procedure?

A discrepancy can be a release defect, an implementation defect, Monte Carlo
error, optimizer variability, or a substantive difference between the printed
model and the code.  Every study must therefore name its semantics explicitly:
`release_compat`, `paper_intended`, or `corrected_common_engine`.

## Stage 0 — data integrity and lossless conversion

Status: conversion implemented and verified; raw-to-summary reconstruction
pending.

Acceptance criteria:

- all source hashes match the vendored release manifest;
- all MATLAB objects used by a study have a documented CSV schema and lineage;
- normalized tables reconstruct the original array/cell relationships;
- all released totals and special mappings pass automated checks;
- a second validation pass reproduces every generated CSV hash.

Remaining Stage 0 task: independently reconstruct the environmental summaries
from the raw event tables.  The release omits the code that creates
`Combined.mat/counts225`, so this is a replication rather than a direct replay.
Both the released and independently rebuilt tables should remain available;
small discrepancies must not be overwritten.

## Stage 1 — deterministic environmental models

Models: GPE, ACT-R, P&A, PPE, MCM, and AMPE.

Run two evaluations for every model:

- `release_compat`: reproduce the MATLAB prototype-history construction,
  probability mapping, truncations, and parameter-vector order;
- `paper_intended`: evaluate the transcribed equations with common notation and
  explicit odds-to-probability conversion.

Tests:

1. Hand-calculated unit cases for one, two, and at least three occurrences.
2. Golden fixtures on cells where the printed equation and release agree.
3. Full 513-cell replay against `published_environmental_surfaces.csv`.
4. Separate attribution tests for known differences: ACT-R's integral
   approximation, P&A's released recursion, MCM's zero truncation, AMPE's
   inclusive range, and prototype schedules for `N>2`.  The MCM replay must
   also update all traces simultaneously from their shared pre-update state.

Primary metrics are RMSE in log probability, squared Pearson correlation in
log probability, and the exact number of compared cells.  A deterministic
release replay should agree at floating-point tolerance (target `1e-10` for
direct array values where algorithms are identical).  Rounded paper tables
are checked only to their printed precision.

## Stage 2 — behavioral AMPE studies

Fit the 14 experiments using the normalized condition/schedule relations.
Keep the source `range_input` separate from presentation ages.

Required modes:

- `release_compat`: reproduce `Appendix.m`, including the hard-coded first
  Rumelhart row and equal averaging of the two Pavlik–Anderson schedules;
- `paper_intended`: additionally apply the stated Rumelhart guessing rule,
  \(p_\mathrm{correct}=1/3+(2/3)p_\mathrm{recall}\);
- `robust_fit`: use deterministic multistart optimization with declared bounds
  and seeds.

Report parameter estimates, RMSE, squared correlation, optimizer convergence,
start sensitivity, and prediction residuals for every experiment.  The
released `starts14` values are optimizer starts, not final results.  Because
exact final predictions were not released, Table 2 is a rounded validation
target rather than a bitwise target.

Additional checks:

- the eight chance-only Rumelhart conditions never receive invented schedules;
- the 64 Pavlik–Anderson condition predictions each average exactly two
  schedule variants of weight 0.5;
- zero, repeated, and fractional ages remain unchanged;
- plotted subsets never replace the full modeling dataset.

## Stage 3 — Hick-law analysis

First replay Equation 16 from the ten *cached* microenvironment odds and the
response times in the CSVs.  Preserve the released single-start Nelder–Mead
result as `matlab_single_start`.  The parameters in `hick_micro_fit.csv`
belong to this empirical microenvironment fit—not AMPE.

Reconstructing those microenvironment odds is a separate extraction study.
They use 12/18/24/30/36-text windows with frequency fixed at six and a
repeated/nonrepeated split.  The release provides the resulting ten odds but
not the code that constructs them, so raw-corpus agreement is an end-to-end
target rather than an artifact-replay prerequisite.

Separately use the 32 × 32 environmental surface, history weights, released
interpolation, and alternative-set frequencies to replay the paper's full
1,000-text environment row, then replay each model row in Table 3.  This
prevents exact recovery of either cached input from being mistaken for a
raw-stream reconstruction.

Then run a robust bounded multistart fit.  This is essential because the
released GPE comparison is a poor local/boundary solution: alternative starts
reduce its latency RMSE substantially.  Report the released and robust results
side by side; never replace one with the other.

Acceptance anchors:

- exact Equation-16 replay from the ten released microenvironment odds;
- independently reconstructed micro-window odds, with any deviation from the
  cache reported rather than silently calibrated away;
- exact recovery of the cached empirical-microenvironment transform parameters
  and selected 2/4/6 predictions;
- explicit reproduction status for the full 1,000-text environment row;
- Table 3 agreement to printed precision for `matlab_single_start`;
- an optimizer-stability report across starts for every model.

## Stage 4 — A&M simulator smoke studies

Run small, fast simulations before any fit claim.

Test both exponential and power decay for:

- valid occurrence probabilities or released raw prediction scores, depending
  on mode;
- no NaNs or negative ages;
- sparse-history conditioning and empty-cell handling;
- increasing revival probability increasing recent activity;
- deterministic results for a named seed and shard count;
- invariance to worker count and shard scheduling.

This stage validates mechanics only.  It should not be described as a
replication of Table 1.

Use the normalized `FigureC1.mat` tables as additional diagnostics: frequency
distributions, post-revival decay, and cross-revival occurrence pairs.  First
replay the cached figure transformations exactly; regeneration of every cache
from its original upstream files is impossible because several named
preprocessing artifacts were not released.

## Stage 5 — released A&M simulation replay

Implement the released simulation semantics exactly, including:

- direct desirability-times-retention Bernoulli draws;
- raw, potentially above-one target scores retained for averaging;
- \(p_\mathrm{revival}=\lambda e^{-\lambda}\);
- variant-specific initialization and boundary behavior;
- the released prefilter and 2,000-window convention;
- count-at-least-1,000 simulation masks;
- post-hoc geometric-mean scaling;
- all known exponential/power implementation differences.

Compare cell predictions and Table 1 metrics with the final cache in
`display461.mat`.  Exact random-number replay is impossible because the
authors did not release seeds or search code.  The criterion is statistical
agreement with uncertainty, not bit identity.

## Stage 6 — corrected common-engine A&M study

Hold the simulation engine constant and vary only the declared model feature:

- odds converted with \(p=O/(1+O)\);
- revival probability \(1-e^{-\lambda}\);
- one documented discrete-time observation convention;
- the same initialization, window boundaries, conditioning, and output scale
  for exponential and power decay.

This study estimates how much of the apparent decay-family comparison is due
to the decay law versus release-specific implementation differences.  Its
results are new analyses, not failed replications.

## Stage 7 — parameter-search replication and uncertainty

Use a reproducible random-stream design:

- a named root seed stored in configuration;
- `SeedSequence` with one child stream per model × stage × shard;
- an explicitly splittable stream design; use a counter-based generator such
  as Philox when schedule-independent counter addressing is required (or
  document jump/split semantics for alternatives such as PCG64DXSM);
- common random numbers during parameter comparisons;
- independent held-out streams for final evaluation;
- results invariant to worker count.

Use deterministic multistart for non-simulation models and retain every start,
termination code, objective value, and parameter vector.  For simulation
models, stream sufficient counts to disk instead of retaining full histories.

Begin with approximately 20 independent pilot shards.  Estimate cell-level
Monte Carlo standard errors and the simulation contribution to overall
log-RMSE.  Increase computation adaptively until a prespecified precision
criterion is met—for example, Monte Carlo contribution below 0.005 RMSE and
key-cell confidence-interval widths below a declared tolerance.  Whether those
intervals cover the released cache is a replication outcome, not a stopping
rule; systematic model differences must be allowed to remain visible.  Only
then decide whether the authors' very large original scale is necessary.

## Study outputs

Each executable study should eventually produce:

- one versioned configuration file;
- an input-manifest hash;
- environment and git revision metadata;
- tidy condition/cell predictions;
- tidy parameter and optimizer traces;
- validation metrics with cell counts;
- seed/shard metadata and Monte Carlo uncertainty where applicable;
- a machine-readable status (`pass`, `partial`, `expected_difference`, or
  `unresolved`);
- a short human-readable report linking every discrepancy to evidence.

No study result should be committed as an unexplained binary workspace or
overwritten in place.  The proposed repository locations for these artifacts
are described separately and will not be applied until approved.

# Minimal Anderson et al. (2023) replication plan

## Goal

The goal is only to establish that the released data are read correctly and
that the main environmental models are correctly understood before further
model comparisons are done in Python or R. This is not a full audit or a
replication of every analysis in the paper.

## 1. Validate the imported data

Already implemented in this PR:

- verify the vendored files against `external/anderson_2023/SHA256SUMS`;
- convert released MATLAB arrays to documented CSV tables without changing
  values or relationships;
- validate table shapes, keys, counts, and manifest hashes;
- keep generated large CSVs out of Git while retaining a reproducible manifest.

One later data task remains: independently rebuild `Combined.mat/counts225`
from the released raw event streams. The producer and original tokenizer were
not released, so this is a validation of our reconstruction, not part of the
lossless MAT-to-CSV conversion.

## 2. Validate the main deterministic models

Models: GPE, ACT-R, P&A, PPE, MCM, and AMPE.

For each model:

1. Test hand-calculable histories with one, two, and three occurrences.
2. Reproduce the released MATLAB calculation on representative cells, including
   its history approximation where the fitter constructs a prototype history.
3. Compare the Python release calculation with all published environmental
   cells used for Table 1.
4. Keep a corrected model separate only where the source check found a material
   issue:
   - full P&A decay recursion and odds conversion;
   - PPE activation mapped to $\alpha\exp(B)$ rather than treated as odds;
   - MCM's corrected Equation A8 dummy index while retaining its released
     simultaneous truncated update.
5. Report RMSE and $r^2$ in log probability for the release comparison. Do not
   reuse the released fitted scale for a corrected equation; refit that equation
   before comparing model quality.

Success means that differences from MATLAB are either numerical noise or a
single, named correction with a direct source citation. No general compatibility
framework is needed.

## 3. Smoke-test the stochastic A&M models

Before any expensive simulation:

- run small exponential and power simulations with declared seeds;
- check finite probabilities, valid event histories, and revival behavior;
- test both the paper's odds mapping and the released direct-value mapping;
- verify conditional history aggregation on a small hand-built fixture.

Only if later work needs the A&M surfaces should we run enough seeded Monte
Carlo batches to compare with the released cache and report sampling error. The
authors' multibillion-history run is not a prerequisite for using or extending
the deterministic models.

## Out of scope

- reaction-time and Hick-law fits;
- the 14 behavioral-study refits;
- SAM reconstruction;
- a paper-wide error audit;
- speculative repository or workflow infrastructure.

## Commands for the current checkpoint

```bash
python3 -m unittest discover -s tests -v
python3 scripts/sanity_check_anderson_2023.py
python3 scripts/prepare_anderson_2023_data.py --validate-only
git lfs fsck
```

The CSV validation command requires a complete generated snapshot. Partial
development builds must remain in a separate output directory.

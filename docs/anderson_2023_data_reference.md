# Anderson et al. (2023) data reference

## Purpose

This document defines the analysis-ready representation of the authors' MATLAB release.  The conversion is deterministic and does not alter the vendored source files.  It separates four things that the release stores together:

1. corpus event streams;
2. environmental aggregates and compressed history cells;
3. behavioral condition means and model-encoded schedules; and
4. cached predictions and fitted parameter vectors.

These distinctions are necessary for modeling.  In particular, a binned probability is not a raw observation, a weighted history cell is not a full occurrence schedule, and a cached prediction is not input data.

The normative column types, keys, nullability rules, controlled values, and variable-level lineage are in `anderson_2023_csv_schema.md`.

## Rebuilding the CSVs

The current converter requires Python, NumPy, and SciPy.  The Reddit v7.3 reader itself uses only the standard library and NumPy; it deliberately does not require MATLAB or `h5py`.  A project-wide locked environment is part of the repository proposal and has not yet been imposed.

From the repository root:

```bash
python scripts/prepare_anderson_2023_data.py
python scripts/prepare_anderson_2023_data.py --validate-only
```

The exporter first verifies every input against `external/anderson_2023/SHA256SUMS`.  It writes files atomically, checks semantic invariants during conversion, and records row counts, byte sizes, source lineage, and SHA-256 hashes in `data/derived/anderson_2023/manifest.csv`.

For a quicker development pass, the two largest layers can be omitted:

```bash
python scripts/prepare_anderson_2023_data.py \
  --skip-corpus --skip-exact-history --skip-appendix-c \
  --output-dir /tmp/anderson_2023_smoke
```

Skip modes refuse to run into a destination that already contains one of the excluded tables.  This prevents a partial manifest from silently validating a mixed full/partial snapshot; use a fresh or separate directory.

## CSV conventions

- UTF-8, comma-delimited, one header row, and LF line endings.
- `snake_case` field names.
- Missing values are empty fields, never the strings `NA`, `NaN`, or `null`.
- Booleans are `true` or `false`.
- Floating-point values use enough significant digits for a lossless binary64 round trip.
- Source, text, condition, bin, and presentation indices are 1-based when they correspond to MATLAB or paper indexing.
- IDs are stable strings; numeric word IDs are always scoped by `corpus_id`.
- Tables contain scalar columns only.  MATLAB cell arrays and list-valued schedules are normalized into related tables.
- Vocabulary strings are literal data.  Ten released tokens are named `na`, `nan`, `null`, or `none`; readers that apply default missing-value word lists can corrupt them.  For example, read `corpus_vocabulary.csv` with pandas using `keep_default_na=False` and an explicit string type for `token`.

## Table catalog

### Corpus event layer

| Table | Rows | Unit | Primary relationship |
|---|---:|---|---|
| `corpus_sources.csv` | 1,502 | source stream | one row per Twitter account history or Reddit subreddit-day |
| `corpus_vocabulary.csv` | 60,000 | retained term | 20,000 terms for each of Twitter, Reddit April, and Reddit May |
| `corpus_texts.csv` | 2,162,837 | text event | preserves empty Reddit comments because they contribute event distance |
| `corpus_occurrences.csv` | 24,342,574 | retained unique term in a text | joins to vocabulary within `corpus_id` |

Use `(corpus_id, source_id, source_text_index)` for temporal operations. `corpus_text_index` is only a surrogate key and must not be used to construct lags across source boundaries.  `position_in_text` preserves the column order of the released matrix; because the authors removed function words and collapsed duplicates, it is not a character offset or an index into raw text.

The two Reddit days have independent vocabularies, so they are represented as `reddit_apr23` and `reddit_may5`.  Only 65 shared strings retain the same numeric ID across the two dictionaries.  The 501 subreddit names occur in the same order on both dates, but the paper treats them as 1,002 separate source-day streams.

### Environmental modeling layer

| Table | Rows | Contents |
|---|---:|---|
| `bin_definitions.csv` | 37 | the 32 square-number bins and five grouped spacing bins |
| `frequency_recency_bins.csv` | 3,072 | Twitter, combined, and exactly inferred Reddit frequency-by-recency cells |
| `spacing_recency_bins.csv` | 3,072 | Twitter, combined, and inferred Reddit twice-occurring spacing cells |
| `spacing_figure_cells.csv` | 576 | the six Figure 2-style curves for each corpus |
| `environmental_fit_targets.csv` | 672 | all 32 × 21 target positions; 513 have `fit_included=true` |
| `exact_history_cells.csv` | 10,404,456 | sparse nonzero cells from `counts225` with their integer weights |
| `three_occurrence_range_bins.csv` | 2,048 | raw recency × range aggregates for Twitter and Reddit |
| `three_occurrence_range_groups.csv` | 320 | the five released grouped range curves |
| `range_spacing_histories.csv` | 128,332 | histories with 3–5 occurrences and their later outcome counts |
| `range_frequency_histories.csv` | 582,736 | source-level range/frequency histories for Twitter and Reddit |

The range/frequency table carries the same `corpus_id` and `source_id` keys as `corpus_sources.csv`, so source metadata can be joined without reconstructing IDs from cell order.

### Appendix-C diagnostic layer

| Table | Rows | Contents |
|---|---:|---|
| `appendix_c_frequency_distributions.csv` | 675 | environmental, released power-A&M, and fitted negative-binomial frequency distributions |
| `appendix_c_decay_replicates.csv` | 2,311,902 | released empirical/exponential/power normalization and delay-count arrays |
| `appendix_c_revival_pairs.csv` | 9,364,833 | paired extra-occurrence counts for the environmental and old power-A&M curves |
| `appendix_c_no_revival_trials.csv` | 2,219,920 | released trials underlying the new/no-revival comparison |

These tables are a lossless scalar expansion of `FigureC1.mat`, which supports the paper's Appendix-C diagnostics.  Several upstream files named inside that MAT file's notes were not released, so the cached arrays can be replayed but not all can be regenerated from their original preprocessing code.  Only the first two columns of `resultsNoRevival` are used by `FigureC1.m`; the release does not define the meanings of columns 3–8, so the CSV preserves them under neutral `released_column_*` names rather than inventing semantics.

`exact_history_cells.csv` represents the axes

```text
(most-recent age, second-most-recent age, frequency N)
```

and a `history_count` weight.  For `N=1`, the released second-age index of 1 is a sentinel and the semantic second age and spacing are left blank.  For `N>2`, occurrences older than the two most recent are absent.  This table is therefore not a collection of complete schedules and contains no cell-level future-hit count.  The released fitting code supplies model-specific prototype histories for the missing occurrences.

In `environmental_fit_targets.csv`, `matlab_linear_index` retains the original 32 × 21 array indexing.  The spacing panel comes from MATLAB columns 16–21, whereas its modeling-friendly `condition_index` deliberately restarts at 1–6. Use `(panel, recency_bin, condition_index)` as the semantic key.

The three corpus labels in the binned tables have distinct provenance:

- `twitter`: directly released in `twitterData.mat`;
- `combined`: directly released in `Combined.mat`;
- `reddit_derived`: reconstructed exactly as `combined - twitter`, because `redditData.mat` does not contain the advertised summary variables.

### Behavioral schedule layer

| Table | Rows | Contents |
|---|---:|---|
| `behavioral_experiments.csv` | 14 | experiment metadata and source labels |
| `behavioral_conditions.csv` | 353 | released aggregate condition probabilities |
| `behavioral_schedule_variants.csv` | 409 | one or more model schedules associated with a condition |
| `behavioral_presentations.csv` | 2,010 | one occurrence age per schedule event |
| `behavioral_plot_points.csv` | 352 | the many-to-many mapping used by the released plotting code |
| `behavioral_fit_starts.csv` | 14 | five optimizer starting values per experiment |

The behavioral release contains aggregate condition probabilities only.  It does not contain participant-level responses, successful recalls, trial denominators, or standard errors.  `observed_probability` must not be treated as a binomial proportion with an invented denominator.

`age_events` and `range_input` are separate released model inputs.  They can use different clocks in mixed within-day/between-day designs and must not be derived from one another.  Ages preserve zeros, repeated pseudo-exposures, fractional values, and source order.

Two mappings require special attention:

- Rumelhart has eight chance-level conditions with no schedule.  The release fixes their predictions at `1/3`.
- Pavlik and Anderson has two schedules per observed condition.  The release evaluates both and averages them with weights `0.5, 0.5`.

`starts14` contains optimizer starts, not final fitted estimates.  Exact final behavioral parameters and predictions were not released.

### Reference-result layer

| Table | Rows | Contents |
|---|---:|---|
| `environmental_model_parameters.csv` | 30 | released parameter vectors in named, code-order form |
| `published_environmental_surfaces.csv` | 6,048 | canonical data/model probability surfaces |
| `hick_observations.csv` | 6 | response-time observations |
| `hick_microenvironment_odds.csv` | 10 | released repeated/nonrepeated odds for 2–6 alternatives |
| `hick_frequency_recency_surface.csv` | 1,024 | released high-frequency surface and history weights |
| `hick_micro_fit.csv` | 3 | released latency-transform parameters |

The canonical environmental cache uses `display461.mat` for exponential A&M, power A&M, and AMPE, because that undocumented file reproduces Table 1.  It uses `display46.mat` for the deterministic models and spacing arrays missing from `display461.mat`.  The stale A&M arrays in `display46.mat` are retained in the immutable release but deliberately excluded from the canonical CSV.

## Verified invariants

The exporter or test suite checks all of the following:

- source-file SHA-256 checksums;
- 1,029,655 Twitter texts and 7,812,969 retained incidences;
- 1,133,182 Reddit comments and 16,529,605 retained incidences;
- exact per-term reconciliation of every streamed Reddit occurrence against the two released count vectors;
- no duplicate retained term ID within an event and zero padding only at row ends;
- `counts225` shape, 10,404,456 nonzero cells, and total weight 2,384,082,074;
- bit-exact reconstruction of the released 32 × 21 `lags` array;
- 513 finite environmental fitting targets;
- the published Table 1 log-RMSE, squared-correlation, and comparison-cell count for all eight models;
- the released row and outcome totals for the 3–5-occurrence range analysis;
- 353 behavioral conditions, 409 schedules, and 2,010 presentations;
- CSV hashes and byte counts on a second `--validate-only` pass.

## Release limitations and discrepancies

These are preserved rather than silently repaired:

1. The stored Twitter arrays contain 1,029,655 tweets, 8,977 fewer than the paper's reported 1,038,632.
2. Released Reddit comments average 14.587 retained unique terms, not the reported 14.85.
3. Forty-seven Reddit comments fill all 100 storage columns and may have been truncated.  The release provides no overflow flag.
4. Reddit timestamps, comment IDs, thread boundaries, and raw text are absent. Twitter's serial dates are exported as both the original MATLAB number and an ISO timestamp using the UTC convention; the release itself does not store an explicit timezone field.
5. The tokenizer and function-word exclusion resources are absent, so the original-text preprocessing stage cannot be reproduced.
6. `redditData.mat` lacks the four summary arrays described in the readme.
7. No released code constructs `Combined.mat/counts225`; rebuilding that representation from the event streams is a new replication task.
8. Rebinning `counts225` differs slightly from the supplied aggregate cells: 191 frequency-recency cells differ (absolute difference sum 65,958), and the `N=2` slice has 18,139 more histories than `counts2`.
9. Range is exclusive (`last-first`) in the 3–5-occurrence and three-occurrence analyses, but inclusive (`last-first+1`) in `allNs1000` and AMPE.
10. `future_count` in both range-history tables is a count in the observation window, not a binary probability.  `future_hit` is an explicit derived indicator.
11. Figure 4's range-bin labels and its MATLAB bin edges disagree by one for several boundaries.
12. `getNs.m` and `freqRange.m` omit the final logically valid 2,000-event window for every source.
13. `fit4s.mat` is intentionally not treated as canonical data. It is a saved exploratory four-parameter/alternative-currency fit whose active equation differs from the final AMPE model and remains available in the immutable vendor snapshot.

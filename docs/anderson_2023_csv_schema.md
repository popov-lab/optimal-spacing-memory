# Anderson et al. (2023) CSV schema, version 1.0.0

This is the normative column-level companion to `anderson_2023_data_reference.md`.  Types are `string`, `integer`, `number`, and `boolean`; a trailing `?` means the CSV field may be empty.  All other fields are required.  Primary keys (PK), foreign keys (FK), controlled values, and transformation lineage are stated explicitly.

## Shared conventions and controlled values

- All source-derived indices are 1-based integers.
- `corpus_id`: `twitter`, `reddit_apr23`, or `reddit_may5` in raw relational tables.  Aggregate tables can instead use `combined`, `reddit_derived`, or the family-level labels documented for that table.
- `source_id`: `{corpus_id}_source_{source_index:03d}`.
- `probability`: number in `[0,1]`; empty only when its denominator is zero or the released cache has a missing cell.
- `hits`, `trials`, `count`, `history_count`, and `future_count`: nonnegative integers.
- A blank field is the only missing-value representation.
- String values that resemble missing markers are not missing.  In particular, `na`, `nan`, `null`, and `none` are legitimate vocabulary tokens; disable default NA-word inference when reading that table.

## Corpus event tables

### `corpus_sources.csv`

PK: `(corpus_id, source_id)`.  Alternate unique key: `(corpus_id, source_index)`.

| Column | Type | Meaning |
|---|---|---|
| `corpus_id` | string | vocabulary/event-stream namespace |
| `source_id` | string | stable source-stream identifier |
| `source_index` | integer | cell index in the released array |
| `source_name` | string | Twitter handle or subreddit name |
| `text_count` | integer | number of stored texts for the source |
| `eligible_1001` | boolean | whether at least one 1,001-event history/target window exists |

Lineage: `twitterData.mat::{allTweets,tweeters}` and `redditData.mat::{arraysApr23,arraysMay5,subredditsApr23,subredditsMay5}`.

### `corpus_vocabulary.csv`

PK: `(corpus_id, word_id)`.

| Column | Type | Meaning |
|---|---|---|
| `corpus_id` | string | vocabulary namespace |
| `word_id` | integer | released 1-based term ID |
| `token` | string | released retained string |
| `retained_occurrence_count` | integer | texts containing that retained term |

Lineage: `twitterData.mat::{twitter20000,countTweets}` and the two `redditData.mat::{vocab*,counts*}` pairs.  Reddit count vectors are also used as full-stream validation targets.

### `corpus_texts.csv`

PK: `(corpus_id, source_id, source_text_index)`.  `corpus_text_index` is a globally unique surrogate key, not a temporal key across sources.

| Column | Type | Meaning |
|---|---|---|
| `corpus_id` | string | event-stream namespace |
| `source_id` | string | FK to `corpus_sources` |
| `source_text_index` | integer | temporal index within source |
| `corpus_text_index` | integer | exporter-assigned global surrogate |
| `matlab_datenum` | number? | original Twitter serial date; blank for Reddit |
| `timestamp_utc` | string? | ISO conversion of Twitter date; blank for Reddit |
| `unique_token_count` | integer | nonzero retained IDs in the released row |

Lineage: one row per row of `allTweets`, `arraysApr23`, or `arraysMay5`; empty Reddit rows are retained.

### `corpus_occurrences.csv`

PK: `(corpus_id, source_id, source_text_index, position_in_text)`.  FKs: the first three columns to `corpus_texts`, and `(corpus_id, word_id)` to `corpus_vocabulary`.

| Column | Type | Meaning |
|---|---|---|
| `corpus_id` | string | event-stream namespace |
| `source_id` | string | source FK |
| `source_text_index` | integer | text FK |
| `corpus_text_index` | integer | redundant fast join to `corpus_texts` |
| `position_in_text` | integer | nonzero released matrix-column order, not raw-text offset |
| `word_id` | integer | vocabulary FK |

## Environmental tables

### `bin_definitions.csv`

PK: `(bin_scheme, bin_index)`.  `bin_scheme` is `square32` or `spacing5`.

| Column | Type | Meaning |
|---|---|---|
| `bin_scheme` | string | bin family |
| `bin_index` | integer | 1-based bin |
| `lower_inclusive` | integer | smallest represented value |
| `upper_inclusive` | integer | largest represented value |
| `representative_value` | number | released midpoint/plot value |

Lineage: `base32.mat::{bounds,means,bounds5}`.

### `frequency_recency_bins.csv`

PK: `(corpus, recency_bin, frequency_bin)`.

| Column | Type | Meaning |
|---|---|---|
| `corpus` | string | `twitter`, `combined`, or `reddit_derived` |
| `recency_bin`, `frequency_bin` | integer | square-bin indices |
| `recency_lower`, `recency_upper`, `frequency_lower`, `frequency_upper` | integer | inclusive bounds |
| `recency_midpoint`, `frequency_midpoint` | number | released representatives |
| `hits`, `trials` | integer | target occurrences and eligible histories |
| `probability` | number? | `hits/trials` |
| `in_fitted_frequency_domain` | boolean | frequency bin 1–15 |
| `meets_fit_count_threshold` | boolean | fitted domain and at least 5,000 trials |

Lineage: `Combined.mat::{hitsA,countsA}`, `twitterData.mat::{hitsA,countsA}`; `reddit_derived` is the exact elementwise difference.

### `spacing_recency_bins.csv`

PK: `(corpus, recency_bin, spacing_bin)`.  Columns are `corpus`, the five recency fields and four spacing fields analogous to the table above, then `hits: integer`, `trials: integer`, and `probability: number?`.

Lineage: released `hits2/counts2` arrays, with Reddit again computed as combined minus Twitter.

### `spacing_figure_cells.csv`

PK: `(corpus, recency_bin, condition_index)`.

Columns: `corpus: string`; recency bin/bounds/midpoint; `condition_index: integer`; `condition_type: string` (`frequency_one_reference` or `twice_occurring_spacing`); `condition_lower: integer`; `condition_upper: integer`; `hits: integer`; `trials: integer`; `probability: number?`; and `fit_included: boolean`.

Lineage: exact pooling specified by `displayPair.m`; condition 1 copies the frequency-one curve and conditions 2–6 pool the five spacing ranges.

### `environmental_fit_targets.csv`

PK: `target_id`; alternate unique key `(panel, recency_bin, condition_index)`.

Columns: `target_id: string`; `matlab_linear_index: integer`; `panel: string` (`frequency_recency` or `spacing_recency`); recency bin/bounds/midpoint; condition index/type/bounds; `probability: number?`; `fit_included: boolean`.

Lineage: `Combined.mat::lags`, independently reconstructed bit-exactly from the aggregate hit/count arrays using `displayPair.m`.

### `exact_history_cells.csv`

PK: `(frequency, most_recent_age, raw_second_age_index)`.

| Column | Type | Meaning |
|---|---|---|
| `frequency` | integer | occurrences in the preceding 1,000 events, 1–225 |
| `most_recent_age` | integer | age of latest occurrence, 1–1,000 |
| `second_most_recent_age` | integer? | age of previous occurrence; blank for singleton |
| `raw_second_age_index` | integer | released dimension-2 index; singleton sentinel is 1 |
| `spacing` | integer? | second age minus most-recent age |
| `history_count` | integer | number of histories represented by the cell |

Lineage: nonzero cells of `Combined.mat::counts225`.  Older ages and cell-level future outcomes were not released.

### Three-occurrence range tables

`three_occurrence_range_bins.csv` has PK `(corpus, recency_bin, range_bin)` and columns: `corpus`; recency bin/bounds/midpoint; range bin/bounds/midpoint; `range_convention` (always `exclusive_last_minus_first`); `hits`; `trials`; and nullable `probability`.

`three_occurrence_range_groups.csv` has PK `(corpus, recency_bin, range_group)` and replaces range-bin fields with `range_group`, `range_lower`, and `range_upper`; it also contains exact pooled `hits`, `trials`, nullable `probability`, and `fit_included`.

Lineage: `threesTwitter.mat` and `threesReddit.mat`; grouped denominators and probabilities are reconstructed from the released 32-bin arrays and checked against `groups5`.

### `range_spacing_histories.csv`

PK: `history_id`.

Columns: `history_id: string`; `frequency: integer`; `gap_1` and `gap_2: integer`; nullable `gap_3` and `gap_4`; `range_exclusive: integer` (sum of gaps); `future_count: integer`; `future_hit: boolean` (derived as count > 0); and `corpus: string` (always `unknown_combined`).

Lineage: row-preserving normalization of `results3to5.mat::{results3,results4,results5}`.  Source/corpus identity was discarded before release.

### `range_frequency_histories.csv`

PK: `history_id`.  FK: `(corpus_id, source_id)` to `corpus_sources`.

Columns: `history_id`; `corpus_id`; `corpus_family` (`twitter` or `reddit`); `snapshot` (`release`, `apr23`, or `may5`); `source_id`; `source_index`; `range_inclusive`; `frequency`; `future_count`; and derived `future_hit`.

Lineage: row-preserving normalization of `allNs1000.mat::{resultsTW,resultsRE}`; the `notes` variable establishes April-then-May Reddit cell order.

## Behavioral tables

### `behavioral_experiments.csv`

PK: `experiment_id`.

Columns: `experiment_id: string`; `source_order: integer`; `released_name: string`; `paper_group: string` (`one_day`, `between_days`, or `mixed`); `response_measure: string` (always `aggregate_probability`); `x_axis_label: string?`; `age_clock: string`; `range_clock: string`; `source_file: string`; and `notes: string?`.

### `behavioral_conditions.csv`

PK: `condition_id`.  FK: `experiment_id`.

Columns: condition/experiment IDs; MATLAB linear/row/column indices; `observed_probability: number`; `fixed_prediction: number?`; `schedule_aggregation: string` (`fixed`, `single`, or `equal_weight_mean`); and `n_schedule_variants: integer`.

### `behavioral_schedule_variants.csv`

PK: `schedule_id`.  FK: `condition_id`.

Columns: schedule/condition IDs; `variant_index: integer`; `mixture_weight: number`; `source_pattern_index: integer`; `range_input: number`; `n_presentations: integer`; and `encoding_provenance: string`.

### `behavioral_presentations.csv`

PK: `(schedule_id, presentation_index)`.  FK: `schedule_id`.

Columns: `schedule_id: string`, `presentation_index: integer`, and `age_events: number`.  Zero, repeated, fractional, and source-ordered values are valid.

### `behavioral_plot_points.csv`

PK: `plot_point_id`; FKs: `experiment_id`, `condition_id`.

Columns: the three IDs; `display_row: integer`; `display_series_index: integer`; `display_series_label: string?`; `display_x: number?`; and `provenance: string`.  The relation is many-to-many because the released P&A display repeats some conditions and omits others.

### `behavioral_fit_starts.csv`

PK/FK: `experiment_id`.  Remaining numeric columns are `b`, `t_prior`, `range_prior`, `threshold`, and `noise_scale`; `parameter_status` is always `optimizer_start`.

Lineage for all behavioral tables: `experiments14.mat::{names14,data14,patterns14,gaps14,labels14,xAxis14, xLabels14,starts14}` plus the special mapping/aggregation logic in `Appendix.m`.

## Appendix-C diagnostic tables

### `appendix_c_frequency_distributions.csv`

PK: `(frequency, series)`.  Columns: `frequency: integer` (1–225), `series: string` (`environmental_data`, `power_am_release_simulation`, or `fitted_negative_binomial`), and `probability: number`.

Lineage: long-form expansion of `FigureC1.mat::probs`.

### `appendix_c_decay_replicates.csv`

PK: `(series, replicate_index, statistic, delay_events)` with the blank delay forming the unique baseline row.  Columns: `series` (`environmental_data`, `exponential_am_release_simulation`, or `power_am_release_simulation`); `replicate_index`; `statistic` (`baseline_normalization_count` or `future_count`); `delay_events: integer?`; and `count: integer`.

Lineage: long-form expansion of `FigureC1.mat::{decayData,decayExp,decayPower}`.  Each replicate has one baseline and delays 1–200.

### `appendix_c_revival_pairs.csv`

PK: `(series, corpus_id, source_cell_index, within_cell_index)`.

Columns: key fields plus `extra_occurrences_period_1: integer` and `extra_occurrences_period_2: integer`.  `series` is `environmental_data` or `power_am_old_release_simulation`; `corpus_id` identifies Twitter, either Reddit snapshot, or simulation.

Lineage: `FigureC1.mat::{revivalsTwitter,revivalsApr23,revivalsMay5, revivalsModel}`.

### `appendix_c_no_revival_trials.csv`

PK: `(simulation_item, within_item_index)`.  The first two outcome columns are `extra_occurrences_period_1` and `extra_occurrences_period_2`.  Six additional `number` fields are deliberately named `released_column_3` through `released_column_8` because neither their semantics nor their generator was released.

Lineage: lossless expansion of `FigureC1.mat::resultsNoRevival`.

## Reference-result and Hick tables

### `environmental_model_parameters.csv`

PK: `(model_id, parameter_index)`; alternate unique key `(model_id, parameter_name)`.  Columns: `model_id`, `parameter_index`, `parameter_name`, `value: number`, `status` (always `released_cached_fit`), and `source_variable`.  Lineage: named vectors in `modelParams.mat`.

### `published_environmental_surfaces.csv`

PK: `(panel, model_id, recency_bin, condition_index)`.  Columns: `panel`, `model_id`, `recency_bin`, `condition_index`, `probability: number?`, and `cache_source`.  Lineage: deterministic arrays from `display46.mat`; final A&M/AMPE arrays from `display461.mat`.

### Hick tables

- `hick_observations.csv`: PK `(n_alternatives, repetition_status)`; fields `n_alternatives: integer`, `repetition_status: string` (`repeated` or `nonrepeated`), `rt_ms: integer`.
- `hick_microenvironment_odds.csv`: the same PK and first two fields, with `odds: number`.
- `hick_frequency_recency_surface.csv`: PK `(corpus, recency_bin, frequency_bin)`; `corpus` is `combined`; it contains bin bounds/midpoints, `probability: number`, and `history_count: integer`.
- `hick_micro_fit.csv`: PK `parameter_name`; values are `intercept_ms`, `scale_ms`, and `power`, with `value: number` and `status: released_cached_fit`.

Lineage: `Schneider.mat::{schneiderData,odds,fullprobs,params6}` plus `Combined.mat::countsA` and `base32.mat`.  `params6` is the empirical microenvironment latency fit, not AMPE.

## Manifest

`manifest.csv` has one row per data table and is intentionally not self-listed.  Columns are `schema_version`, `table_id`, `filename`, `rows`, `bytes`, `sha256`, `source_files`, and `notes`.  `filename` and `table_id` are unique.  Validation requires the output directory's complete CSV inventory to equal the listed filenames plus `manifest.csv`.  It also rejects unexpected regular files other than the directory's `README.md`, so interrupted-write artifacts cannot hide outside the `.csv` inventory.  The validator checks the exact manifest header and schema version, unique table/file IDs, declared row and byte counts, LF termination, and every SHA-256 digest.

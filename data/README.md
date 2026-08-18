# Optimal spacing gap by test delay — Cepeda et al. meta-analysis data

## Contents

- `cepeda2006_optimal_gaps.csv` — one row per experiment: the test delay (retention
  interval) and the observed optimal gap (interstudy interval), both in days.

## What the data are

Each point is one spacing experiment in which the gap between study sessions was
varied at a fixed test delay, and the optimal gap — the gap producing the best
final-test performance, flanked by shorter and longer tested gaps — could be
identified. The core of the dataset is the study-level database of the Cepeda,
Pashler, Vul, Wixted, & Rohrer (2006) meta-analysis, and it also includes the
same lab's later experiments: Cepeda et al. (2009) Experiments 1 and 2 and the
four retention-interval conditions of Cepeda et al. (2008).

This is the dataset plotted (on log–log axes) in Figure 5 of Cepeda et al.
(2009) and in the meta-analysis figure of Mozer et al. (2009) (Figure 3 in the
NeurIPS proceedings; Figure 4 in some manuscript versions). The repo's
`cepeda.png` shows the Cepeda et al. (2009) version.

## Provenance and extraction

The values were recovered directly from the embedded vector graphics of the
Mozer et al. (2009) figure in the proceedings PDF — i.e., they are the exact
plotted coordinates, not a raster digitization. All 46 circle markers were
decomposed from the figure's drawing commands and calibrated against its
log10 gridlines (calibration residuals < 0.001 decades); remaining uncertainty
from PDF coordinate rounding is about ±0.2%. Recovered values fall on round
experimental units (e.g., 5.787e-5 days = 5 s; 0.998 ≈ 1 day), confirming they
reproduce the underlying data values. The `source` column attributes points to
the meta-analysis or to the later experiments; among coincident points the
assignment is by composition (the coordinates are identical).

## Data dictionary

| Column             | Type    | Description |
|--------------------|---------|-------------|
| `point`            | integer | Data point identifier, 1–46, ordered by increasing test delay. |
| `test_delay_days`  | numeric | Retention interval between the final study session and the final test, in days. |
| `optimal_gap_days` | numeric | Gap (interstudy interval) between learning sessions that produced the best observed final-test performance, in days. |
| `test_delay_readable`  | string | `test_delay_days` re-expressed as the nearest round value in seconds/minutes/hours/days, for readability (approximate; within the ±0.2% recovery noise for nearly all points). |
| `optimal_gap_readable` | string | `optimal_gap_days` re-expressed the same way. |
| `source`           | string  | Origin of the data point: the Cepeda et al. (2006) meta-analysis database, or the lab's later experiments (Cepeda et al., 2008; Cepeda et al., 2009). |

## Notes

- Several points coincide exactly, so the figures show 43 distinct markers for
  46 points: (7, 1) is plotted three times and one sub-minute point twice.
- Cepeda et al. (2009) Figure 5 shows this dataset minus the three
  Cepeda et al. (2008) points at test delays 35, 70, and 350 days. Its text
  reports n = 48 points, which cannot be fully reconciled with the recoverable
  markers; the difference is presumably additional exact overplotting.
- Cepeda et al. (2009) Experiments 2a and 2b both had an observed optimal gap
  of 28 days at a 168-day delay; the source figure plots a single circle there.
- The gaps plotted for the Cepeda et al. (2008) delays do not all match that
  paper's reported observed recall optima (1, 11, 21, 21 days at delays 7, 35,
  70, 350); e.g., the plotted 7-day gap at the 35-day delay matches its
  recognition optimum.

## References

> Cepeda, N. J., Pashler, H., Vul, E., Wixted, J. T., & Rohrer, D. (2006).
> Distributed practice in verbal recall tasks: A review and quantitative
> synthesis. *Psychological Bulletin*, *132*(3), 354–380.
> https://doi.org/10.1037/0033-2909.132.3.354

> Cepeda, N. J., Vul, E., Rohrer, D., Wixted, J. T., & Pashler, H. (2008).
> Spacing effects in learning: A temporal ridgeline of optimal retention.
> *Psychological Science*, *19*(11), 1095–1102.
> https://doi.org/10.1111/j.1467-9280.2008.02209.x

> Cepeda, N. J., Coburn, N., Rohrer, D., Wixted, J. T., Mozer, M. C., &
> Pashler, H. (2009). Optimizing distributed practice: theoretical analysis
> and practical implications. *Experimental Psychology*, *56*(4), 236–246.
> https://doi.org/10.1027/1618-3169.56.4.236

> Mozer, M. C., Pashler, H., Cepeda, N. J., Lindsey, R. V., & Vul, E. (2009).
> Predicting the optimal spacing of study: A multiscale context model of
> memory. In *Advances in Neural Information Processing Systems 22*
> (pp. 1321–1329).

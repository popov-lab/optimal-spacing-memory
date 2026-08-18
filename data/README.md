# Spacing-effect data recovered from the Cepeda / Mozer figures

## Contents

- `cepeda2006_optimal_gaps.csv` — one row per experiment: the test delay (retention
  interval) and the observed optimal gap (interstudy interval), both in days.
  Plotted by `src/plot_optimal_gaps.py` → `figures/optimal_gap_by_test_delay.png`.
- `mozer2009_figure2_recall.csv` — one row per recall observation from panels a–d
  of Figure 2 in Mozer et al. (2009): the forgetting and spacing functions of
  Cepeda et al. (2009) Experiments 1, 2a, and 2b, and the spacing functions of
  the four retention-interval conditions of Cepeda et al. (2008).
  Plotted by `src/plot_mozer2009_figure2.py` → `figures/mozer2009_fig2{a,b,c,d}.png`.

Both datasets were recovered directly from the embedded vector graphics of the
Mozer et al. (2009) proceedings PDF — i.e., they are the exact plotted
coordinates, not raster digitizations. Marker centers were decomposed from the
figures' drawing commands and calibrated against the axes (calibration
residuals < 0.03 units), leaving only PDF coordinate-rounding noise. Where the
recovered values coincide with round units (integers, round day counts) the
round values are reported; the remaining values are reported to one decimal.

## `cepeda2006_optimal_gaps.csv`

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
`cepeda.png` shows the Cepeda et al. (2009) version. All 46 circle markers were
recovered from the vector figure; the coordinates fall on round experimental
units (5 s, 90 s, 1 day, 168 days, ...), so the round values are the underlying
data and the numeric columns give their exact conversion to days. The single
value that does not correspond to a round unit (the 32.2 s test delay of
point 9) is kept as recovered. The `source` column attributes points to the
meta-analysis or to the later experiments; among coincident points the
assignment is by composition (the coordinates are identical).

### Data dictionary

| Column             | Type    | Description |
|--------------------|---------|-------------|
| `point`            | integer | Data point identifier, 1–46, ordered by increasing test delay. |
| `test_delay_days`  | numeric | Retention interval between the final study session and the final test, in days — the exact day conversion of `test_delay_readable`. |
| `optimal_gap_days` | numeric | Gap (interstudy interval) between learning sessions that produced the best observed final-test performance, in days — the exact day conversion of `optimal_gap_readable`. |
| `test_delay_readable`  | string | The same value in natural units (seconds/minutes/hours/days). |
| `optimal_gap_readable` | string | The same value in natural units. |
| `source`           | string  | Origin of the data point: the Cepeda et al. (2006) meta-analysis database, or the lab's later experiments (Cepeda et al., 2008; Cepeda et al., 2009). |

### Notes

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

## `mozer2009_figure2_recall.csv`

One recall percentage per row, from panels a–d of Figure 2 in Mozer et al.
(2009) (empirical points only; the panels' model fits are not part of the
dataset, and panel e repeats panel d's data with a different model). In each
experiment, subjects studied the material in two sessions separated by a gap
(the interstudy interval, ISI) and took a final test after a retention interval
(RI) following the second session:

- **Spacing function** — recall on the final test, as a function of ISI at
  fixed RI.
- **Forgetting function** — recall at the start of the second session, i.e.,
  after a single study exposure; the study–test lag of that measurement is the
  ISI itself, so `ri_days` is left empty for these rows.

Panels a–c are Cepeda et al. (2009) Experiment 1 (Swahili–English word pairs,
RI = 10 days), Experiment 2a (obscure facts, RI = 168 days), and Experiment 2b
(object names, RI = 168 days). Panel d is Cepeda et al. (2008): all 26
gap × RI conditions of the study (RIs 7, 35, 70, and 350 days), spacing
functions only.

### Data dictionary

| Column       | Type    | Description |
|--------------|---------|-------------|
| `panel`      | string  | Source panel in Mozer et al. (2009) Figure 2: `a`–`d`. |
| `experiment` | string  | Experiment the observation comes from. |
| `material`   | string  | Study material of that experiment. |
| `function`   | string  | `spacing` (final-test recall) or `forgetting` (recall at the start of session 2). |
| `isi_days`   | numeric | Interstudy interval (gap between the two study sessions), in days, as plotted. |
| `ri_days`    | numeric | Retention interval between the second session and the final test, in days; empty for forgetting rows (the measurement's delay equals `isi_days`). |
| `recall_pct` | numeric | Percent correct recall. |

### Notes

- ISI 0 denotes the massed condition; the actual gap was ~5 min in
  Cepeda et al. (2009) Experiment 1, ~20 min in Experiment 2, and within-session
  in Cepeda et al. (2008).
- Recall values for panels b, c, and the panel-a spacing function are exact
  integers in the source vectors; the remaining values carry one decimal.
- The panel-a spacing values are consistent with the pairwise differences in
  Table 1 of Cepeda et al. (2009) except for two small discrepancies (gaps
  2→4: plotted 1.0 vs. reported 0.4; gaps 1→14: plotted 9.0 vs. reported 8.4).

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

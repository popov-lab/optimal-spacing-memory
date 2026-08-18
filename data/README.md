# Spacing-effect data recovered from published figures

## Contents

- `cepeda2006_optimal_gaps.csv` — one row per experiment: the test delay (retention
  interval) and the observed optimal gap (interstudy interval), both in days.
  Plotted by `src/plot_optimal_gaps.py` → `figures/optimal_gap_by_test_delay.png`.
- `cepeda_spacing_recall.csv` — one row per recall observation: the forgetting
  and spacing functions of Cepeda et al. (2009) Experiments 1, 2a, and 2b, and
  the spacing functions of the four retention-interval conditions of
  Cepeda et al. (2008).
  Plotted by `src/plot_spacing_recall.py` → `figures/spacing_recall_{a,b,c,d}.png`
  (linear x axis) and `figures/spacing_recall_{a,b,c,d}_log.png`
  (logarithmic x spacing; linear below ISI = 1 day so the massed condition stays on the axis).

`cepeda2006_optimal_gaps.csv` was recovered directly from the embedded vector
graphics of the Mozer et al. (2009) proceedings PDF — i.e., it holds the exact
plotted coordinates, not raster digitizations. Marker centers were decomposed
from the figure's drawing commands and calibrated against the axes (calibration
residuals < 0.03 units), leaving only PDF coordinate-rounding noise. Where the
recovered values coincide with round units (integers, round day counts) the
round values are reported; the remaining values are reported to one decimal.

`cepeda_spacing_recall.csv` is read from the original Cepeda papers; its
provenance is described in its own section below. An earlier version of this
dataset, taken from Figure 2 of Mozer et al. (2009), is kept in `archive/`
together with the figures made from it; `archive/README.md` records why it was
retired.

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

## `cepeda_spacing_recall.csv`

One recall percentage per row, from the two studies. In each experiment,
subjects studied the material in two sessions separated by a gap (the
interstudy interval, ISI) and took a final test after a retention interval
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
functions only. The panel letters are retained from the earlier version of
this dataset so that the two are directly comparable.

### Provenance

The two studies are read from different sources, because their own figures
differ in quality.

**Cepeda et al. (2009) — panels a–c.** Read directly from Figures 3 and 4 of
that paper. Those figures are clean: markers are large, the axes are finely
ticked, and every series is drawn at its true gap position. Marker centres
were located from the raster at native resolution — circles and squares by the
midpoint of their upper and lower edges, triangles by their base row less a
fixed base-to-centre offset calibrated on the markers of known value — and
converted with axis calibrations fitted to the frame and tick marks
(residuals < 0.1 percentage points). Extraction precision is about
±0.2 points. Experiment 1 values carry one decimal; every clean Experiment 2
marker fell within 0.16 of an integer, so those are reported as integers.

Two caveats. At ISI = 2 days in Experiment 1 the forgetting and spacing
markers coincide, so the forgetting value (68.9) is the centre of the merged
pair and is uncertain by about ±0.5. At ISI = 0 and 1 day in Experiment 2 the
markers of adjacent conditions overlap; they were separated by their
individual edges.

**Cepeda et al. (2008) — panel d.** *Not* read from that paper's own figure.
Its Figure 3a plots the massed and 1-day conditions with a horizontal jitter
that is not applied to the error bars, overlaps several markers, and draws the
series as cubic splines, so marker centres there cannot be located reliably.
Instead each condition is the mean of three independent published replots of
the same study:

> Carpenter, S. K., Cepeda, N. J., Rohrer, D., Kang, S. H. K., & Pashler, H.
> (2012). Figure 3.
> Toppino, T. C., & Gerbier, E. (2014). Figure 4.4.
> Carpenter, S. K. (2020). Figure 2.

Each was digitized independently at native raster resolution with its own axis
calibration. The three agree closely — the standard deviation across sources
is below 0.3 points at 23 of the 26 conditions and never exceeds 0.4 — so the
mean is reported to one decimal. Two further replots (Walsh et al., 2018,
Figure 5; Antony et al., 2024, Figure 3A) were digitized as well and agree
with the same values; they are not included in the mean because the Walsh
figure carries a uniform offset of about +1 point and the Antony figure plots
integer-rounded values.

### Data dictionary

| Column       | Type    | Description |
|--------------|---------|-------------|
| `panel`      | string  | Panel of the figure set: `a`–`d` (see above). |
| `experiment` | string  | Experiment the observation comes from. |
| `material`   | string  | Study material of that experiment. |
| `function`   | string  | `spacing` (final-test recall) or `forgetting` (recall at the start of session 2). |
| `isi_days`   | numeric | Interstudy interval (gap between the two study sessions), in days. |
| `ri_days`    | numeric | Retention interval between the second session and the final test, in days; empty for forgetting rows (the measurement's delay equals `isi_days`). |
| `recall_pct` | numeric | Percent correct recall. |

### Notes

- ISI 0 denotes the massed condition; the actual gap was ~5 min in
  Cepeda et al. (2009) Experiment 1, ~20 min in Experiment 2, and within-session
  in Cepeda et al. (2008).
- Error bars are not included. Both papers plot ±1 SEM, but the bars overlap
  the markers throughout the crowded short-gap conditions and could not be
  measured to a useful precision.
- The panel-a spacing values are consistent with the pairwise differences in
  Table 1 of Cepeda et al. (2009) except for two small discrepancies (gaps
  2→4: plotted 1.0 vs. reported 0.4; gaps 1→14: plotted 9.0 vs. reported 8.4).
- Figure 2 of Pashler, Rohrer, Cepeda, & Carpenter (2007) replots the three
  Cepeda et al. (2009) spacing functions as vector art, and its exact
  coordinates match the values recorded here at all 18 points — an independent
  confirmation of the spacing values, and evidence that they predate both
  2009 publications.
- Where this dataset differs from the Mozer et al. (2009) version in
  `archive/`, see `archive/README.md`.

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

> Carpenter, S. K., Cepeda, N. J., Rohrer, D., Kang, S. H. K., & Pashler, H.
> (2012). Using spacing to enhance diverse forms of learning: Review of recent
> research and implications for instruction. *Educational Psychology Review*,
> *24*(3), 369–378. https://doi.org/10.1007/s10648-012-9205-z

> Toppino, T. C., & Gerbier, E. (2014). About practice: Repetition, spacing,
> and abstraction. In B. H. Ross (Ed.), *The Psychology of Learning and
> Motivation* (Vol. 60, pp. 113–189). Academic Press.
> https://doi.org/10.1016/B978-0-12-800090-8.00004-4

> Carpenter, S. K. (2020). Distributed practice or spacing effect. In
> *Oxford Research Encyclopedia of Education*. Oxford University Press.
> https://doi.org/10.1093/acrefore/9780190264093.013.859

> Pashler, H., Rohrer, D., Cepeda, N. J., & Carpenter, S. K. (2007). Enhancing
> learning and retarding forgetting: Choices and consequences.
> *Psychonomic Bulletin & Review*, *14*(2), 187–193.
> https://doi.org/10.3758/BF03194050

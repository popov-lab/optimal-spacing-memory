# Data extracted from Figure 5 of Cepeda et al. (2009)

## Source

Data points were digitized from Figure 5 of:

> Cepeda, N. J., Coburn, N., Rohrer, D., Wixted, J. T., Mozer, M. C., & Pashler, H. (2009).
> Optimizing distributed practice: theoretical analysis and practical implications.
> *Experimental Psychology*, *56*(4), 236–246. https://doi.org/10.1027/1618-3169.56.4.236

Figure 5 is a log–log scatter plot of the optimal gap (interstudy interval that
produced the best final-test performance) against the test delay (retention
interval), for studies in the Cepeda, Pashler, Vul, Wixted, & Rohrer (2006)
meta-analysis in which the optimal gap was flanked by shorter and longer gaps,
plus the two experiments reported in Cepeda et al. (2009).

## Files

- `cepeda2009_figure5.csv` — extracted data, one row per plotted point.

## Data dictionary (`cepeda2009_figure5.csv`)

| Column             | Type    | Description |
|--------------------|---------|-------------|
| `point`            | integer | Identifier for a plotted data point (one study/condition), numbered 1–40 in order of increasing test delay. |
| `test_delay_days`  | numeric | Retention interval between the final study session and the final test, in days. Digitized from the figure and rounded to three significant figures. |
| `optimal_gap_days` | numeric | Spacing gap (interstudy interval) between learning sessions that produced the best final-test performance, in days. Digitized from the figure and rounded to three significant figures. |

## Extraction method and precision

- The figure image (827 × 567 px) was extracted losslessly from the article PDF.
  Markers were located with a Hough circle transform and refined to subpixel
  precision by ring-template correlation (all refinement shifts < 0.5 px).
- Axes were calibrated by least-squares mapping of pixel position to log10
  value using the plot frame and the interior x = 1 and y = 1 axis lines
  (calibration residuals < 0.003 decades on both axes).
- Digitization precision is limited by marker localization (~0.5 px), i.e.
  roughly ±1.3% in test delay and ±1.7% in optimal gap per coordinate; values
  are therefore reported to three significant figures.
- Five coordinates whose marker centers lay within 0.8 px of a drawn gridline
  were snapped to that gridline's exact value: points 31 (delay = 1, gap = 1),
  32 (delay = 1, gap = 0.001), and 33, 36, 39 (gap = 1).
- Validity checks: the two experiments reported in the paper's text are
  recovered accurately — Experiment 1 (10-day test delay, 1-day optimal gap)
  digitizes as (10.2, 1) and Experiment 2 (168-day test delay, 28-day optimal
  gap) digitizes as (170, 27.5). An ordinary least-squares power regression on
  the extracted points gives optimal gap ≈ 0.050 × delay^0.72 (R² = 0.75 in
  log–log space), which reproduces the figure's dashed best-fit line.

## Caveats

- The text of the paper states that Figure 5 contains n = 48 data points, but
  only 40 distinct markers are visible in the figure. The remaining points
  evidently coincide with (are overplotted by) visible markers and cannot be
  recovered from the image; each CSV `point` therefore represents at least one
  study/condition.
- Values are estimates read from a published figure, not the original study
  values. Anyone needing exact values should consult the Cepeda et al. (2006)
  meta-analysis ( https://doi.org/10.1037/0033-2909.132.3.354 ) and the original
  studies it reviews.

# Archive — data and figures derived from Mozer et al. (2009)

The files in this folder were recovered from Figure 2 of

> Mozer, M. C., Pashler, H., Cepeda, N. J., Lindsey, R. V., & Vul, E. (2009). Predicting the optimal spacing of study: A multiscale context model of memory. In *Advances in Neural Information Processing Systems 22* (pp. 1321–1329).

and were the basis of `data/mozer2009_figure2_recall.csv` and `figures/mozer2009_fig2*.png` until they were replaced by `data/cepeda_spacing_recall.csv`, which is read from the original Cepeda papers instead. They are kept here for reference, not for use.

## Why they were retired

Mozer et al. (2009) Figure 2 does not reproduce the values published in the experiments it cites. Comparing it against the original figures and against five later papers that replot the same studies:

- **Cepeda et al. (2008)** (Mozer's panel d): the NeurIPS values differ from the published Figure 3a at roughly half of the 26 conditions, by up to 7 percentage points. Five independent replots — Carpenter et al. (2012), Toppino & Gerbier (2014), Walsh et al. (2018), Carpenter (2020), and Antony et al. (2024) — all agree with the published figure and not with Mozer: the mean absolute difference between any two of those replots is 0.2–1.2 points, whereas each of them differs from Mozer by 2.5–2.7 points on average. Most visibly, the 35-day retention function peaks at a gap of 11 days in every other source (matching the optimum stated in the 2008 paper's own text) but at 7 days in Mozer's panel.
- **Cepeda et al. (2009)** (Mozer's panels a–c): the spacing functions match the source exactly, and Experiment 1 matches at every point. The Experiment 2 forgetting functions do not: at a gap of 28 days Mozer plots the forgetting value equal to the spacing value in both panels (56 for facts, 26 for objects), where the published Figure 4 shows 48 and 25. Three further points differ by one percentage point.

Neither paper reports any exclusions, reanalysis, or preprocessing that would explain the differences, and none of the discrepant values appears anywhere else in the literature.

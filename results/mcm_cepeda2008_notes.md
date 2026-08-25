# MCM replication: Cepeda et al. (2008)

## Why this analysis is post-hoc

Mozer et al. (2009) state that, for each experiment, the four MCM parameters
$\{\mu,\nu,\omega,\xi\}$ were estimated by least-squares fit to the
single-session forgetting function, after which the spacing functions were
predicted with $N=100$ and $\epsilon_r=9$ fixed.

For Cepeda et al. (2008), the published article reports the 26 final spacing
conditions but does not report the Session-2 first-test forgetting function
used by Mozer et al. This analysis was run before that forgetting function was
recovered from Mozer and Lindsey (2016). The recovered observations are now in
`data/cepeda_spacing_recall.csv`, but the results below remain the original
spacing-only analysis and have not been refitted.

As a fallback, we fit one common MCM parameterization directly and jointly to
all 26 spacing observations across retention intervals of 7, 35, 70, and 350
days. The four RI curves are **not** fitted separately; doing so would remove the
cross-RI constraint that is central to the model.

## Joint fit

With $N=100$ and $\epsilon_r=9$ fixed, the least-squares fit gives approximately

$$
\mu = 23.6907,\qquad
\nu = 1.07361,\qquad
\omega = 0.80668,\qquad
\xi = 0.93531.
$$

Across all 26 conditions,

$$
\mathrm{RMSE}=0.04830,
$$

or about **4.83 percentage points**.

The largest residual is the 70-day-RI condition at ISI = 14 days. The empirical
curve dips at 14 days and rises again at 21 days, whereas the fitted MCM spacing
function is smooth and peaks between these observations.

## Interpretation

This analysis answers a different question from the Cepeda et al. (2009)
replications:

- **Cepeda et al. (2009):** fit forgetting, freeze parameters, predict spacing.
- **Cepeda et al. (2008), analysis recorded here:** forgetting data was then
  unavailable, so the spacing surface itself was fitted.

The 2008 result is therefore useful as a descriptive/post-hoc test of whether a
single MCM parameterization can capture the four-RI spacing surface, but it must
not be presented as a replication of Mozer et al.'s parameter-free prediction.

## Reproduction

Run

```bash
python src/fit_mcm_cepeda2008.py
python src/plot_mcm_replication.py
```

The first command writes `results/mcm_cepeda2008_fit.csv` and
`results/mcm_cepeda2008_predictions.csv`. The second writes the corresponding
figures under `figures/`.

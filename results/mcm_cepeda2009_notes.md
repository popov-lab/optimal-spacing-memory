# MCM replication: initial status

This is a first independent implementation of the Multiscale Context Model (MCM) of Mozer et al. (2009), using the Cepeda et al. (2009) data in `data/cepeda_spacing_recall.csv`.

## Implementation

`src/mcm.py` implements the leaky-integrator formulation and exact marginalization described in the MCM paper:

- `N = 100`
- `tau_i = mu * nu^i`
- `gamma_i = xi^i / sum_j xi^j`
- `s_i = sum_{j<=i} gamma_j x_j / sum_{j<=i} gamma_j`
- between study episodes, `x_i` decays exponentially with time constant `tau_i`
- on an encoded study episode, `Delta x_i = epsilon * (1 - s_i)`
- `epsilon = 1` after retrieval failure and `epsilon = epsilon_r = 9` after retrieval success
- encoding succeeds with probability `omega`
- expected recall is computed by explicitly marginalizing over retrieval and encoding outcomes, rather than by Monte Carlo simulation

For the Cepeda experiments, one experimental learning session is treated as one MCM study episode. The within-session test-with-feedback passes are not modeled as separate MCM updates. Reverse-engineering the vector model curves in Mozer et al. Figure 2 confirms that this interpretation reproduces the published MCM spacing curves to plotting precision when the corresponding fitted forgetting curve is used.

## Independent fits to the original Cepeda et al. (2009) data

`src/fit_mcm_cepeda2009.py` fits `{mu, nu, omega, xi}` by multistart least squares to the *forgetting function only*, then freezes the parameters and predicts the spacing function. `epsilon_r = 9` and `N = 100` are fixed.

Current best fits:

| Panel | Material | mu | nu | omega | xi | forgetting SSE | spacing RMSE |
|---|---|---:|---:|---:|---:|---:|---:|
| a | Swahili-English | 0.979869 | 3.718688 | 0.937351 | 0.491356 | 0.00272149 | 9.16 pp |
| b | obscure facts | 1.781635 | 12.845737 | 0.956422 | 0.329904 | 0.000163441 | 6.30 pp |
| c | object names | 5.138786 | 1.051299 | 0.920765 | 0.966589 | 0.000422430 | 5.78 pp |

The full point predictions are in `results/mcm_cepeda2009_predictions.csv`.

## Important reproduction issue

The published MCM Figure 2 does not use exactly the forgetting-function observations reported in the original Cepeda et al. (2009) paper. This is already documented in `archive/README.md`.

For Experiment 2a, the clearest discrepancy is the 28-day forgetting observation: Mozer et al. plot approximately 56%, whereas the original Cepeda figure gives approximately 48%. Fitting the MCM to the Mozer-plotted forgetting values reproduces the published panel-b spacing prediction essentially exactly. Fitting to the original Cepeda values, as done here, necessarily gives a different prediction.

Panel c has an additional unresolved issue. The solid blue forgetting curve plotted by Mozer et al. is substantially worse, in raw least-squares error, than the optimum obtainable from the plotted forgetting points with the stated four-parameter model. This may reflect unreported fitting constraints, a local optimizer solution, or a mismatch between the data used for fitting and the data shown in the figure. It should not yet be treated as a model discrepancy.

## Next checks

1. Generate overlays of the independent fits/predictions against the original Cepeda data and against the Mozer Figure 2 curves.
2. Recover the exact Mozer Figure 2 solid curves from the PDF vector graphics and store them as validation targets.
3. Test sensitivity to the nominal `ISI = 0` approximation versus the actual approximately 5-minute and 20-minute gaps.
4. Then turn to the 2008 experiment; its experimental protocol is known, but the single-session forgetting observations used to fit the MCM are still missing from the published 2008 article.

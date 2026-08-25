# MCM history ablation summary

The strongest available test is the Cepeda et al. (2009) forgetting-constrained prediction.

Pooled spacing RMSE (percentage points):

- Full stochastic MCM: 7.23
- No retrieval branching (use expected learning rate): 16.47
- No encoding branching (use deterministic omega-scaled update): 7.71
- No branching: 19.55

Thus encoding-history branching is not needed for the two-study predictions, whereas replacing retrieval success/failure by its mean update badly damages Experiments 1 and 2a. The reason is specific: with only two study episodes, retrieval branching matters through the nonlinear capped recall readout `min(1, s_N)` when the successful-retrieval update (`epsilon_r = 9`) saturates.

For the Cepeda et al. (2008) data, where only a post-hoc direct spacing fit is possible, the same comparison is much less diagnostic: RMSEs are 4.66, 6.57, 5.09, and 4.93 pp respectively. The deterministic no-branching model can absorb most differences by changing its multiscale parameters.

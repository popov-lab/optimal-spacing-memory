# MCM stochastic-history ablations

## Question

Mozer et al. (2009) marginalize over two stochastic outcomes at study: retrieval success/failure and encoding success/failure. Are those discrete history branches actually necessary for the spacing predictions?

The cleanest test removes the *branching* while preserving its one-step expectation.

## Mean-field ablations

The full model has

- retrieval success with probability $p=s_N$, producing $\epsilon=\epsilon_r=9$; retrieval failure produces $\epsilon=1$;
- encoding with probability $\omega$, otherwise the study event has no effect.

We compare it with three deterministic/mean-field replacements:

1. **No retrieval branching.** Replace the Bernoulli retrieval outcome by its expected learning rate

   $$
   \bar\epsilon = (1-p)\,1 + p\epsilon_r = 1+p(\epsilon_r-1).
   $$

2. **No encoding branching.** Replace encode/no-encode branching by a deterministic expected update

   $$
   \Delta x_i = \omega\epsilon(1-s_i).
   $$

   This preserves the single-study expectation $\omega F(t)$ exactly, so the original forgetting-function fits can be reused without refitting.

3. **No branching.** Apply both mean-field replacements simultaneously.

For Cepeda et al. (2009), all four variants therefore use exactly the same forgetting-constrained parameters. Only the spacing prediction changes.

## Cepeda et al. (2009): predictive test

Spacing RMSE in percentage points:

| Variant | Exp. 1 | Exp. 2a | Exp. 2b | Pooled |
|---|---:|---:|---:|---:|
| Full stochastic | 9.16 | 6.30 | 5.78 | **7.23** |
| No retrieval branching | 21.06 | 18.35 | 5.78 | **16.47** |
| No encoding branching | 8.86 | 8.10 | 5.87 | **7.71** |
| No branching | 27.18 | 19.32 | 5.87 | **19.55** |

The result is asymmetric. Replacing stochastic encoding by its deterministic expectation barely changes predictive accuracy. Replacing the retrieval-success branch by its expected learning rate damages Experiments 1 and 2a severely.

For these two-study experiments, this has a useful algebraic interpretation. Retrieval at Session 2 is the final stochastic branch before the test. If the final readout were linear, averaging the two retrieval-contingent updates before versus after the branch would give the same expectation. The difference arises because MCM uses

$$
P(\mathrm{recall})=\min(1,s_N).
$$

With $\epsilon_r=9$, the successful-retrieval branch can saturate at 1 while the unsuccessful branch does not. Therefore

$$
E[\min(1,S)] \neq \min(1,E[S]),
$$

and the explicit retrieval mixture matters. In Experiment 2b the predicted strengths remain low enough that this saturation is essentially absent, and the full and no-retrieval-branching predictions are numerically identical.

## Hard mechanism removals

A stronger ablation removes the mechanisms rather than only their stochastic branching:

- **No retrieval-dependent update:** set $\epsilon=1$ regardless of retrieval success.
- **Always encode:** set $\omega=1$ and refit the forgetting function with the remaining decay parameters.

Pooled Cepeda et al. (2009) spacing RMSE:

| Variant | RMSE (pp) |
|---|---:|
| Full stochastic | **7.23** |
| No retrieval-dependent update | **16.54** |
| Always encode ($\omega=1$) | **9.95** |
| Both hard removals | **17.08** |

The retrieval-dependent update is therefore important for these predictions. The poorer $\omega=1$ model is less clean evidence for stochastic encoding itself, because fixing $\omega=1$ also removes the amplitude parameter needed to match imperfect immediate encoding. The mean-field encoding ablation above is the more diagnostic test of whether encoding *history branching* is necessary.

## Cepeda et al. (2008): direct post-hoc fit

Because the forgetting observations used by Mozer et al. are unavailable, these variants were also fit directly to all 26 spacing observations. This is a post-hoc flexibility check, not a predictive test.

Using the same broad bounds for all variants ($10^{-4}\leq\mu\leq10^3$, $10^{-3}\leq\nu-1\leq10^{1.5}$, $.45\leq\omega\leq.999$, $.05\leq\xi\leq.999$ where applicable), the best RMSEs found were:

| Variant | RMSE (pp) |
|---|---:|
| Full stochastic | **4.66** |
| No retrieval branching | **6.57** |
| No encoding branching | **5.09** |
| No branching | **4.93** |

Because all parameters are fit directly to the spacing surface and the parameterization is weakly identified, these differences are much less diagnostic. In particular, the fully deterministic mean-field version can absorb much of the difference by changing its multiscale parameters.

## Interpretation

For the available two-study data, the evidence does **not** support treating stochastic encoding histories as an essential part of MCM: deterministic expected encoding performs almost identically when the forgetting fit is held fixed.

The retrieval-contingent component matters much more. But the present result should be stated precisely: with only two study episodes, the need for explicit retrieval branching is tied to the nonlinear capped recall readout and the large $\epsilon_r=9$ update. It is not yet evidence that long stochastic retrieval histories are intrinsically necessary for producing a spacing effect.

A natural next comparison is therefore the deterministic multiscale core with a continuous strength-dependent learning rule, which would put MCM on the same footing as the stripped SAC model.

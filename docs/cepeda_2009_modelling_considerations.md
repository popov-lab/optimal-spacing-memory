# Modelling considerations for Cepeda et al. (2009)

## The complication

The forgetting curves in Cepeda et al. (2009) are easy to misread as memory after one study event. They are actually measured at the beginning of Session 2, after participants completed a multi-event learning procedure in Session 1. The spacing curves are measured after further test-feedback-restudy events in Session 2.

This distinction matters especially for SAC. New learning at event $j$ depends on the strength retrieved immediately before that event:

$$
u_j=\delta[1-B(t_j^-)],
\qquad
B(t)=\sum_{k:t_k<t}u_k f(t-t_k).
$$

Several events within a session therefore do not behave like one event with a larger scale. Each event changes the learning produced by the next, generating higher-order interactions in $\delta$. Once the sessions contain repeated events, $\delta$ can no longer be removed as a common multiplicative scale and retained only as the coefficient on a single two-event interaction.

The first response in Session 2 is particularly informative: it supplies the forgetting observation before feedback, and the same pre-feedback strength determines the subsequent SAC update.

## What is known about the procedure

The details below come from the method section of [Cepeda et al. (2009)](https://doi.org/10.1027/1618-3169.56.4.236).

| | Session 1 | Session 2 | Final test |
|---|---|---|---|
| **Experiment 1: Swahili-English pairs** | One passive presentation of all 40 pairs, 7 s each. Participants then completed test-with-feedback cycles until every item had been translated correctly twice, not necessarily consecutively. Feedback displayed both words for 5 s regardless of accuracy. An item was removed after reaching criterion. | Exactly two complete test-with-feedback cycles. The first response in the first cycle is the forgetting observation and occurs before its feedback. | Session 3 occurred 10 days after Session 2. The nominal zero-day Session 1–2 gap was approximately 5 minutes. |
| **Experiment 2: obscure facts and object names** | After a no-feedback pretest, each 23-item set received one initial exposure followed by three fixed test-with-feedback blocks: four known learning events per item. Initial statements were shown for 13 s; test cues were shown for 13 s and then the answer for 5 s. Facts always preceded objects. | Exactly two test-with-feedback blocks. Again, the first response is the forgetting observation before feedback. | Session 3 occurred 168 days after Session 2. The nominal zero-day Session 1–2 gap was approximately 20 minutes. |

The Experiment 2 pretest is not treated as a learning event: it had no feedback, and correctly answered items were excluded.

## What cannot be reconstructed from the paper

Experiment 2 has known pass counts, and its mean within-session timing can be approximated from the fixed presentation durations. Experiment 1 is less recoverable. The paper does not report:

- the distribution of test trials required to reach criterion;
- item-level response times or cycle durations;
- total Session 1 duration;
- when each item reached criterion relative to the end of the session.

The last point matters because items were removed after criterion. An easy item may have had fewer learning events and a longer effective delay before Session 2 than a difficult item. Event count and event timing are therefore correlated.

Further approximations would also be needed to distinguish passive presentation, retrieval, and feedback-restudy as learning events. The current analyses use the same SAC update rule for all events and treat the pre-feedback response as a measurement rather than an additional update.

## Session-batch approximation already tested

As a first diagnostic, all events within a session were collapsed to the same time while retaining their nonlinear SAC updates. For $m$ effectively massed events,

$$
A_m=1-(1-\delta)^m.
$$

With $m_1$ events in Session 1, $m_2$ events in Session 2, intersession gap $a$, and final retention interval $b$,

$$
B_F(a)=A_{m_1}f(a),
$$

$$
B_S(a,b)=A_{m_1}f(a+b)+A_{m_2}[1-A_{m_1}f(a)]f(b).
$$

The analysis fixed $m_1=4$ for Experiment 2 and $m_2=2$ for every experiment. Because the Experiment 1 count is unknown, it evaluated $m_1=3,\ldots,7$. The three materials had separate $d$ and $\tau$, while $\delta$, the logistic threshold, and the logistic scale were shared.

Under the central $m_1=4$ assumption:

| Fit | Spacing RMSE | Predicted Experiment 1 optimum |
|---|---:|---:|
| Forgetting fit plus Experiment 1 $\delta$ calibration | 5.31 pp | 3.54 days |
| Joint forgetting and spacing fit | 4.27 pp | 5.19 days |
| Spacing-only fit | 2.64 pp | 1.72 days |

Relative to the earlier single-event approximation, the forgetting-trained optimum moved from 12.55 to 3.54 days and its spacing RMSE improved from 7.93 to 5.31 percentage points. The joint optimum moved from 6.53 to 5.19 days. The spacing-only result was essentially unchanged.

Thus repeated learning explains a substantial part of the original mismatch, but not all of it. Even with seven assumed Session 1 events, the joint-fit Experiment 1 optimum was 3.59 days, whereas the spacing-only fit remained near 1.5–1.7 days.

Full assumptions, sensitivity results, parameters, predictions, and figures are in the [session-batch SAC report](../results/sac_cepeda2009_batch_report.md). The reproducible implementation is in [`src/fit_sac_cepeda2009_batch.py`](../src/fit_sac_cepeda2009_batch.py).

## Next steps

1. Check published models of these data, supplementary material, or raw data for an estimate of the Experiment 1 trial-count distribution.
2. Replace the fixed Experiment 1 count with a mixture over event counts, averaging response probabilities rather than latent strengths.
3. Add finite within-session timing. Experiment 2 can use the approximately reconstructable pass schedule; Experiment 1 requires a sensitivity analysis over plausible pass intervals and post-criterion delays.
4. Refit the forgetting-trained, joint, and spacing-only protocols. If the discrepancy remains, then consider separate learning gains for passive study and test-feedback-restudy, or treat retrieval itself as an update.

The immediate unresolved input is therefore not another SAC parameter. It is a defensible approximation to the distribution of Experiment 1 learning histories.

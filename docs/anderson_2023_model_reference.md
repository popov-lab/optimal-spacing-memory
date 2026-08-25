# Anderson et al. (2023) model reference

This page states the environmental models in one notation and distinguishes the
equations we intend to use from behavior specific to the released MATLAB code.
The source comparison is in the [implementation notes](anderson_2023_implementation_notes.md).

## Common notation

For a history of $n$ occurrences, set the most recent occurrence to the origin:

$$
t_1 < t_2 < \cdots < t_n = 0.
$$

The prediction time $t>0$ is therefore the time since the most recent
occurrence. When a model needs the lag between consecutive occurrences, use

$$
\Delta t_k=t_{k+1}-t_k.
$$

The odds after $n$ occurrences are $O_n(t)$ and the corresponding probability is

$$
p_n(t)=\frac{O_n(t)}{1+O_n(t)}.
$$

The same symbol, $\alpha>0$, is used for every model's multiplicative odds
scale. It is not an activation. If activation $B$ is mapped to recall by

$$
p=\frac{1}{1+\exp[(\mu-B)/\sigma]},
$$

then at unit logistic scale, $\sigma=1$,

$$
O=\exp(-\mu)\exp(B)=\alpha\exp(B).
$$

Equivalently, with intercept $\eta=-\mu$, $\alpha=e^\eta$. When $\sigma\ne1$,
retain $\mu$ and $\sigma$ in the response mapping rather than folding them into
the common unit-scale $\alpha$.

The Python functions currently receive the equivalent positive ages
$(t-t_1,\ldots,t-t_n)$. This is an input representation only; the equations
below use the normalized event times throughout.

For the environmental data, time is discrete text position. An occurrence in
the immediately preceding text has age $t-t_n=1$. Multiple appearances within
one text were collapsed.

## 1. General Performance Equation (GPE)

$$
O_n(t)=\alpha n^c t^{-d}.
$$

Here $c$ controls the effect of frequency and $d\ge0$ controls forgetting.
Anderson et al.'s Table 1 reports $\alpha=.021$, $c=.575$, and $d=.608$.

Code: `gpe_odds`.

## 2. ACT-R base-level model

$$
O_n(t)=\alpha\sum_{j=1}^{n}(t-t_j)^{-d}.
$$

Every occurrence contributes a component with the same decay exponent $d$.
Table 1 prints a negative odds multiplier, but the released parameter file and
`ACTRFit.m` use a positive value near $.040$; the implementation therefore uses
$\alpha>0$.

Code: `actr_odds`.

## 3. Pavlik and Anderson (P&A)

Each occurrence has its own decay exponent. Let the first be

$$
d_1=a.
$$

For occurrence $k\ge2$, its exponent depends on activation from all earlier
components at that occurrence:

$$
B_{k-1}(t_k)
=\sum_{j=1}^{k-1}(t_k-t_j)^{-d_j},
$$

$$
d_k=a+cB_{k-1}(t_k).
$$

Prediction then uses

$$
O_n(t)=\alpha\sum_{j=1}^{n}(t-t_j)^{-d_j}.
$$

The released parameter file stores the positive multiplier
$\alpha=.05296011$; this is approximately $e^{-2.94}$, the value printed as
$B=-2.94$ in Table 1. The general implementation evaluates the full sum above.
The released environmental prototype has a narrower calculation for the newest
component, documented in the implementation notes.

Code: `pavlik_anderson_component_decays` and `pavlik_anderson_odds`.

## 4. Predictive Performance Equation (PPE)

Define normalized recency weights and effective elapsed time:

$$
w_j(t)=\frac{(t-t_j)^{-x}}
{\sum_{k=1}^{n}(t-t_k)^{-x}},
$$

$$
T_n(t)=\sum_{j=1}^{n}w_j(t)(t-t_j).
$$

For $n>1$, spacing determines the decay exponent:

$$
d_n=b+\frac{m}{n-1}\sum_{k=1}^{n-1}
\frac{1}{\log(\Delta t_k+e)}.
$$

For $n=1$, set $d_1=b$. Walsh et al. define the power product as memory
activation:

$$
B_n(t)=n^c T_n(t)^{-d_n}.
$$

With unit-scale logistic noise, the corrected odds are therefore

$$
O_n(t)=\alpha\exp\{B_n(t)\}.
$$

Anderson et al. instead label activation as odds in Equation 6 and implement

$$
O_n^{\mathrm{Anderson}}(t)=\alpha B_n(t).
$$

These are kept separate in code: `ppe_activation`, corrected `ppe_odds`, and
`anderson_ppe_odds` for the Anderson equation and MATLAB comparison. Anderson
et al.'s fitted $\alpha=.018$ belongs to their mapping and must not be treated
as a fitted scale for the corrected odds equation.

## 5. Multiscale Context Model (MCM)

This is Anderson et al.'s environmental adaptation of Mozer et al. (2009). It
uses $K=100$ traces with

$$
\tau_i=\mu\nu^i,
\qquad
\gamma_i=\frac{\omega\xi^i}{\sum_{k=1}^{K}\xi^k},
\qquad i=1,\ldots,K.
$$

Thus $\sum_i\gamma_i=\omega$. Between two events separated by $\Delta t$,

$$
x_i\leftarrow x_i\exp(-\Delta t/\tau_i).
$$

At a later occurrence, the release computes every increment from the same
pre-update state:

$$
\bar x_i=\frac{\sum_{k=1}^{i}\gamma_kx_k}
{\sum_{k=1}^{i}\gamma_k},
\qquad
\Delta x_i=\max(0,1-\bar x_i),
$$

$$
x_i\leftarrow x_i+\Delta x_i.
$$

After decay from $t_n=0$ to prediction time $t$, define

$$
q_n(t)=\min\left(\sum_i\gamma_ix_i(t),.999999\right),
$$

and

$$
O_n(t)=\alpha\frac{q_n(t)}{1-q_n(t)}.
$$

The `max` truncation is present in `MCMFit.m` but absent from printed Equation
A6. Anderson et al.'s Table 1 reports $\mu=.032$, $\nu=1.111$,
$\omega=.704$, $\xi=.978$, and $\alpha=.029$.

Code: `mcm_state` and `mcm_odds`.

## 6. Anderson-Milson (A&M) environmental process

For each item $i$, draw a persistent desirability $\pi_i$ from a Gamma
distribution and a persistent decay $d_i$ from an exponential distribution.
Revivals reset the item's latent elapsed time $u$.

The two retention functions are

$$
r_{\exp}(u;d_i)=e^{-d_i u},
\qquad
r_{\mathrm{pow}}(u;d_i)=u^{-d_i},
\qquad u=1,2,\ldots.
$$

The coherent latent odds and occurrence probability are

$$
Q_i(u)=\pi_i r(u;d_i),
\qquad
p_i(u)=\frac{Q_i(u)}{1+Q_i(u)}.
$$

For a supplied history summary $h$, prediction is the Monte Carlo conditional
mean of $p_i(u)$ among simulated targets with that history. Only after this
conditioning is the common output scale $\alpha$ applied to the conditional
odds.

The released MATLAB scripts instead compare a uniform draw directly with
$Q_i(u)$ and average the raw $Q_i(u)$ values within history cells. The simulator
therefore retains an explicit `released_probability` option for checking the
release while using the odds mapping by default.

The fitted output scales reported in Table 1 are $\alpha=.704$ for exponential
decay and $\alpha=.724$ for power decay. They are post-hoc scales, not
parameters of the latent generator.

Code: `simulate_anderson_milson`,
`anderson_milson_conditional_predictions`, and
`scale_anderson_milson_predictions`.

## 7. AMPE

Include a prior age $t_P>0$ in the currency estimate:

$$
T_n(t)=\mathrm{HM}(t-t_1,\ldots,t-t_n,t_P)+1.
$$

The released environmental implementation uses the inclusive discrete range

$$
G_n=t_n-t_1+1=1-t_1.
$$

Define the effective interval and decay:

$$
M_n=\frac{G_n+g_P}{2},
\qquad
d_n=\frac{b}{M_n}.
$$

The odds are

$$
O_n(t)=\alpha\frac{n}{M_n}T_n(t)^{-d_n}.
$$

Anderson et al.'s Table 1 reports $\alpha=214$, $b=1401$,
$t_P=15.18$, and $g_P=1565$. AMPE is undefined for $n=0$.

Code: `inclusive_range`, `ampe_components`, and `ampe_odds`.

### Human-recall mapping

For human recall, remove the environmental scale from the activation:

$$
B_n(t)=\log n-\log M_n-\frac{b}{M_n}\log T_n(t).
$$

With threshold $\mu$ and logistic scale $\sigma$,

$$
p_n(t)=\frac{1}{1+\exp[(\mu-B_n(t))/\sigma]}.
$$

At $\sigma=1$, this gives $\alpha=e^{-\mu}=e^\eta$ with $\eta=-\mu$, as
in the common mapping above. The released behavioral schedules sometimes use
age zero; `ampe_recall_probability` preserves that released convention without
changing the environmental time convention.

## Checks in this PR

Run

```bash
python3 -m unittest discover -s tests -v
python3 scripts/sanity_check_anderson_2023.py
```

The unit tests include hand-checkable histories for the deterministic models,
the corrected and Anderson PPE mappings, AMPE's spacing crossover, and small
seeded A&M simulations. No reaction-time or SAM replication is in scope.

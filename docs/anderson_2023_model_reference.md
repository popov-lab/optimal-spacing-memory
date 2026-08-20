# Anderson et al. (2023) model reference

Implementation-complete transcription of the models compared in:

> Anderson, J. R., Betts, S., Byrne, M. D., Schooler, L. J., & Stanley, C. (2023). The environmental basis of memory. *Psychological Review, 130*, 1137-1166. https://doi.org/10.1037/rev0000409

This reference states the models in common notation. It excludes the paper's binning, optimization, and corpus preprocessing because those are not part of the models. Discrepancies between the equations and the authors' MATLAB release are recorded in [the release audit](anderson_2023_release_audit.md).

## Common notation and output

Let occurrence times and the prediction time satisfy

$$
u_1 < \cdots < u_N < t.
$$

Define

$$
r_j=t-u_j>0
$$

as the age of occurrence $j$ at prediction, and

$$
\ell_j=u_{j+1}-u_j>0
$$

as the lag between consecutive occurrences. The most recent age is

$$
R=\min_j r_j.
$$

For the environmental data, time is discrete text position: an occurrence in the immediately preceding text has age 1, and occurrences in adjacent texts have lag 1. Multiple appearances within one text were collapsed, so ages are distinct.

All deterministic environmental models first produce odds $O$, then

$$
p=\frac{O}{1+O}.
$$

The implementation exposes the odds directly and uses one shared `odds_to_probability` function. All scales and times must be positive. Decay and spacing-sensitivity parameters are nonnegative unless stated otherwise. The released human-AMPE schedules are the one exception: they sometimes encode a just-presented event with age zero, as specified in Section 8.

## 1. General Performance Equation (GPE)

Paper Equation 3:

$$
O=A N^c R^{-d}.
$$

| Quantity | Meaning |
|---|---|
| $A>0$ | odds scale |
| $c$ | frequency exponent |
| $d\ge0$ | recency-decay exponent |
| $N$ | number of occurrences |
| $R$ | age of the most recent occurrence |

Printed Table 1 fit: $A=.021, c=.575, d=.608$.

Code: `gpe_odds(ages, params)`.

## 2. ACT-R base-level model

Paper Equation 4:

$$
O=\kappa\sum_{j=1}^{N}r_j^{-d}.
$$

Here $d\ge0$ is a common component-decay exponent and $\kappa>0$ is an odds scale. The paper calls the latter $b$, but that symbol is overloaded in later models.

Printed Table 1 gives $d=.792$ and $b=-.040$. A negative multiplier would give negative odds. The released code resolves this as a typographical sign error: it uses a positive scale near $.040$.

Code: `actr_odds(ages, params)`.

## 3. Pavlik and Anderson (P&A)

Each occurrence creates a component whose decay exponent depends on the strength of all earlier components at that occurrence.

First occurrence (Equation A1):

$$
d_1=a.
$$

For occurrence $n\ge2$, compute prior strength at $u_n$:

$$
H_n=\sum_{j=1}^{n-1}(u_n-u_j)^{-d_j},
$$

then assign the new component (Equation A2):

$$
d_n=a+cH_n.
$$

At prediction (Equation 5):

$$
O=\kappa\sum_{j=1}^{N}r_j^{-d_j}.
$$

| Parameter | Meaning |
|---|---|
| $a\ge0$ | minimum component decay |
| $c\ge0$ | increase in a new component's decay per unit current strength |
| $\kappa>0$ | odds scale |

Printed Table 1 gives $a=.758, c=.444, B=-2.94$. The only coherent reading of $B$ is an additive log-odds intercept, so

$$
\kappa=e^B=e^{-2.94}\approx .0529.
$$

These P&A estimates were fitted with the release's direct-probability output and its defective $N>2$ prototype recursion. They are reference values, not fitted estimates for the corrected equation-level implementation below.

Code: `pavlik_anderson_component_decays` and `pavlik_anderson_odds`.

## 4. Predictive Performance Equation (PPE)

Define recency weights and effective elapsed time (Equation A3):

$$
w_j=\frac{r_j^{-x}}{\sum_{k=1}^{N}r_k^{-x}},
\qquad
T_{\mathrm{eff}}=\sum_{j=1}^{N}w_jr_j.
$$

Equivalently,

$$
T_{\mathrm{eff}}
=\frac{\sum_j r_j^{1-x}}{\sum_j r_j^{-x}}.
$$

At $x=1$, this is the harmonic mean; as $x\to\infty$, it approaches the most recent age.

For $N>1$, spacing determines the decay exponent (Equation A4):

$$
d_{\mathrm{PPE}}
=b+\frac{m}{N-1}\sum_{j=1}^{N-1}
\frac{1}{\log(\ell_j+e)}.
$$

For $N=1$, the necessary special case is $d_{\mathrm{PPE}}=b$.

The odds equation is Equation 6:

$$
O=A N^c T_{\mathrm{eff}}^{-d_{\mathrm{PPE}}}.
$$

| Parameter | Meaning |
|---|---|
| $A>0$ | odds scale |
| $c$ | frequency exponent |
| $x$ | recency-weighting exponent |
| $b\ge0$ | minimum decay |
| $m\ge0$ | spacing sensitivity |

Printed Table 1 fit:

$$
x=8.699,\quad c=.612,\quad b=.549,\quad m=.186,\quad A=.018.
$$

Code: `ppe_components` and `ppe_odds`.

## 5. Multiscale Context Model (MCM), environmental variant

This is Anderson et al.'s simplified environmental adaptation of Mozer et al. (2009), not the original full recall model. It uses $K=100$ traces. Trace index $i$ is unrelated to occurrence index $j$.

Time constants and trace weights (Equations A7-A8):

$$
\tau_i=\mu\nu^i,
\qquad
\gamma_i=\frac{\omega\xi^i}{\sum_{k=1}^{K}\xi^k},
\qquad i=1,\ldots,K.
$$

Thus $\sum_i\gamma_i=\omega$. Intended constraints are

$$
\mu>0,\quad \nu>1,\quad 0<\omega<1,\quad 0<\xi<1.
$$

All trace strengths are initialized to one at the first occurrence. Between events (Equation A5):

$$
x_i(t+\Delta t)=x_i(t)e^{-\Delta t/\tau_i}.
$$

At every later occurrence, calculate simultaneously from the pre-update state

$$
\bar x_i
=\frac{\sum_{k=1}^{i}\gamma_kx_k}{\sum_{k=1}^{i}\gamma_k},
\qquad
\Delta x_i=\max(0,1-\bar x_i),
\qquad
x_i\leftarrow x_i+\Delta x_i.
$$

The `max` is present in the released MATLAB but omitted from printed Equation A6. After decay to prediction, define

$$
q=\min\left(\sum_i\gamma_ix_i,.999999\right).
$$

Equation 7 maps this recall-like strength to environmental odds:

$$
O=A\frac{q}{1-q}.
$$

Printed Table 1 fit:

$$
\mu=.032,\quad \nu=1.111,\quad \omega=.704,\quad \xi=.978,\quad A=.029.
$$

Code: `mcm_state` and `mcm_odds`.

## 6. Anderson-Milson (A&M) environmental process

For each item $i$, draw once

$$
\pi_i\sim\operatorname{Gamma}(k_\pi,\theta_\pi),
\qquad
d_i\sim\operatorname{Exponential}(\text{mean }\mu_d).
$$

The Gamma uses shape $k_\pi$ and scale $\theta_\pi$, so $E[\pi_i]=k_\pi\theta_\pi$. The item retains both $\pi_i$ and $d_i$ across revivals.

Revivals follow a homogeneous Poisson process with mean interval $\mu_R$, equivalently rate $\lambda_R=1/\mu_R$. A revival resets latent elapsed time but is not caused by an observed occurrence. At discrete event boundaries the probability of at least one revival in the preceding unit interval is

$$
p_R=1-e^{-1/\mu_R}.
$$

Use age $u=1$ for the first sampled text after a revival. This removes the printed power law's singularity at zero while using one timing convention for both variants:

$$
r_{\exp}(u;d)=e^{-du},
\qquad
r_{\mathrm{pow}}(u;d)=u^{-d},
\qquad u=1,2,\ldots.
$$

The intended latent need odds and occurrence probability are

$$
O_i(u)=\pi_i r(u;d_i),
\qquad
p_i(u)=\frac{O_i(u)}{1+O_i(u)},
\qquad
Y_{i,t}\sim\operatorname{Bernoulli}(p_i(t)).
$$

The Bernoulli and odds-conversion steps are required by the paper's interpretation but are not printed. The released code instead compares a uniform draw directly with $z_i=\pi_i r(u;d_i)$, which implicitly caps the Bernoulli probability at one, but it retains the raw, potentially above-one $z_i$ as the target prediction. The simulator supports both explicitly; `occurrence_mapping="odds"` is the coherent default and `"released_probability"` is the exact MATLAB scoring behavior.

The model's prediction is not the unconditional $p_i(t)$. For a supplied history summary $h$, it is the Monte Carlo conditional mean

$$
\widehat p(h)
=\frac{1}{|C_h|}\sum_{(i,t)\in C_h}p_i(t),
\qquad
C_h=\{(i,t):H_{i,t}=h\}.
$$

The paper's exact cells use frequency plus the ages of the two most recent occurrences in the preceding 1,000 events; later binning is an analysis step, not part of the model. Code: `anderson_milson_conditional_predictions`.

In released mode the same operation averages raw $z_i(t)$, not bounded probabilities, so a conditional cell can in principle exceed one.

Only after that conditional averaging is output scale $A$ applied. The paper says it scales odds,

$$
\widehat p_A(h)
=\frac{A\widehat p(h)/(1-\widehat p(h))}
{1+A\widehat p(h)/(1-\widehat p(h))}.
$$

The MATLAB release instead multiplies its raw conditional cell value directly and chooses

$$
A_{\mathrm{release}}
=\exp\!\left[
\operatorname{mean}\log p_{\mathrm{observed}}
-\operatorname{mean}\log\widehat z
\right].
$$

These alternatives are implemented by `scale_anderson_milson_predictions` and `released_geometric_mean_scale`. Output scale never changes which histories are generated.

| Model | $k_\pi$ | $\theta_\pi$ | $\mu_d$ | $\mu_R$ | $A$ |
|---|---:|---:|---:|---:|---:|
| Exponential A&M | .164 | .139 | .035 | 333 | .704 |
| Power A&M | .199 | .482 | 4.076 | 800 | .724 |

The values of $A$ in this table are post-hoc released prediction scales, not generative parameters for the corrected odds-based simulator. Likewise, the other values were fitted under the released direct-probability generator and should be treated as audit/reference values when using the corrected mapping.

The paper calls $\beta$ a revival *rate* but reports 333 and 800. The released code uses their reciprocals as per-event rates, confirming that the reported values are mean intervals.

Code: `simulate_anderson_milson`. Its returned arrays include latent odds, probabilities, Bernoulli occurrences, revival indicators, elapsed times, and item parameters. This is the latent generator; the conditional function above provides the operational prediction stage.

## 7. AMPE environmental model

For an observed history with $N\ge1$, include one prior pseudo-age $t_P>0$ in the currency estimate (Equation 10):

$$
T
=\operatorname{HM}(r_1,\ldots,r_N,t_P)+1
=\frac{N+1}{\sum_{j=1}^{N}1/r_j+1/t_P}+1.
$$

The released implementation uses the inclusive range

$$
G=\max_j r_j-\min_j r_j+1,
$$

so a singleton has $G=1$ and adjacent occurrences have $G=2$. Define the effective interval (Equation 11), inferred decay (Equation 12), and initial desirability (Equation 13):

$$
M=\frac{G+g_P}{2},
\qquad
d=\frac{b}{M},
\qquad
\pi=\frac{aN}{M}.
$$

Substitution into Equation 9 gives

$$
O=\frac{aN}{M}T^{-b/M}.
$$

For stable computation,

$$
\log O=\log a+\log N-\log M-\frac{b}{M}\log T.
$$

Parameters are $a,b,t_P,g_P>0$. Printed Table 1 fit:

$$
a=214,\quad b=1401,\quad t_P=15.18,\quad g_P=1565.
$$

AMPE is undefined for $N=0$; the paper supplies no unseen-item base rate.

Code: `inclusive_range`, `ampe_components`, and `ampe_odds`. `ampe_components` accepts an explicit range because the behavioral applications sometimes use a different range calendar from the ages entering $T$.

## 8. AMPE human-recall mapping

Let $A_{\mathrm{mem}}=\log O$. With logistic activation noise of scale $s>0$ and retrieval threshold $\tau$, Equation 14 is

$$
P(\text{recall})
=\frac{1}{1+\exp[(\tau-A_{\mathrm{mem}})/s]}.
$$

The desirability scale $a$ is absorbed into the threshold. Define

$$
\alpha
=\log\left(\frac{NT^{-d}}{M}\right)
=\log N-\log M-\frac{b}{M}\log T
$$

and $\eta=\tau-\log a$. Equation 15 becomes

$$
P(\text{recall})
=\frac{1}{1+\exp[(\eta-\alpha)/s]}.
$$

The five fitted parameters per experiment are

$$
b>0,\quad t_P>0,\quad g_P>0,\quad s>0,\quad \eta\in\mathbb R.
$$

Code: `ampe_recall_probability`.

Unlike the environmental window convention, the released behavioral schedules may assign age zero to a just-presented item and may repeat that zero. The behavioral implementation therefore accepts $r_j\ge0$. If any age is zero, MATLAB's harmonic mean is zero and Equation 10 gives $T=1$. Environmental AMPE continues to require positive, distinct text ages.

Schedule adapters needed to reproduce individual behavioral datasets:

- Single-day studies: each study or test opportunity is one event.
- Between-day studies: the paper uses a 3-1-2 pseudo-exposure rule (first session = 3 observations, same-day repeat = 1, later-day repeat = 2), but does not fully specify pseudo-event placement or aging.
- Mixed studies: one day is 500 age events, while range is the sum of within-day ranges and ignores intervening days.
- Rumelhart's three-alternative task adds guessing:

  $$
  P(\text{correct})=P(\text{recall})+\frac{1-P(\text{recall})}{3}.
  $$

These dataset-specific encodings are not silently built into the core model.

## 9. Reaction-time mapping

Paper Equation 16 maps any positive model odds to latency:

$$
RT=I+S O^{-q},
$$

with $S,q\ge0$. Code: `odds_to_reaction_time`.

## 10. SAM status

Search of Associative Memory (SAM) appears only in the behavioral comparison, whose statistics are copied from Walsh et al. (2018). Anderson et al. provide no equations, parameters, initialization, output mapping, or simulation code for SAM, and their release contains no SAM implementation. It is therefore not implemented here. Substituting an arbitrary version of Raaijmakers's SAM would not be a faithful implementation of the model underlying Table 2.

## Confirmatory checks only

Run

```bash
python -m unittest discover -s tests -v
python scripts/sanity_check_anderson_2023.py
```

The sanity script performs no fitting and reads no author data. It checks only:

1. finite probabilities and a recency decline for all six deterministic environmental models;
2. AMPE's intended spacing crossover using the released environmental fit;
3. valid sparse outputs and a post-revival boost in short simulations of both A&M decay variants.

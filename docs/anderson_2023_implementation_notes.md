# Anderson et al. (2023) implementation notes

This is a focused source comparison for implementing the environmental models
and data pipeline in Python. It is not a general audit of the paper.

Claims below were checked against all three primary artifacts:

- [the published paper](https://act-r.psy.cmu.edu/wordpress/wp-content/uploads/2023/02/EnvironmentalBassisofMemory.pdf),
  identified by equation or appendix number;
- the immutable MATLAB files in `external/anderson_2023/matlab/`, linked to
  exact lines below;
- the released `modelParams.mat`, read directly as a MATLAB v5 file.

## Differences that affect model implementation

| Model | Verified paper-to-code result | Python treatment |
|---|---|---|
| GPE | Paper Equation 3 and [`GPEFit.m` lines 3-8](../external/anderson_2023/matlab/GPEFit.m#L3-L8) both use a frequency power, recency power, positive multiplier, and odds-to-probability conversion. | Implement Equation 3 directly. |
| ACT-R | Paper Equation 4 and Table 1 print the scale as $b=-.040$, which would make odds negative. The saved parameter is positive $.04010094$, and [`ACTRFit.m` lines 3 and 14-15](../external/anderson_2023/matlab/ACTRFit.m#L3-L15) multiply by that positive value before converting odds to probability. | Use $\alpha>0$; the released executable value resolves the printed sign inconsistency. |
| P&A | Paper Equations 5 and A1-A2 require the newest decay to depend on all earlier components. In [`PAFit.m` lines 64-76](../external/anderson_2023/matlab/PAFit.m#L64-L76), the prototype decays are reversed and the newest component uses only `decays(end)`. [`PAFit.m` line 16](../external/anderson_2023/matlab/PAFit.m#L16) then uses the scaled sum directly, with no odds conversion. The saved multiplier is positive $.05296011$, approximately $e^{-2.94}$. | Evaluate the full A2 sum and return odds. Use the MATLAB function only as a release-comparison target. |
| PPE | Walsh et al.'s PPE defines $n^cT^{-d}$ as activation and maps activation through a logistic response. Anderson et al. relabel that activation as odds in Equation 6. [`PPEFit.m` lines 13-15](../external/anderson_2023/matlab/PPEFit.m#L13-L15) apply the multiplier and odds conversion; [`PPEFit.m` lines 58-75](../external/anderson_2023/matlab/PPEFit.m#L58-L75) calculate $n^cT^{-d}$. | `ppe_activation` returns the Walsh activation; `ppe_odds` returns corrected $\alpha\exp(B)$ for unit logistic scale; `anderson_ppe_odds` retains Anderson's mapping for comparison. |
| MCM | Paper Equation A6 omits truncation. [`MCMFit.m` lines 74-85](../external/anderson_2023/matlab/MCMFit.m#L74-L85) computes all increments from the pre-update state and applies `max(0, 1-strength)`. [`MCMFit.m` lines 10-12](../external/anderson_2023/matlab/MCMFit.m#L10-L12) also confirms that Equation A8's denominator must sum over its dummy index. | Follow the simultaneous, truncated released update and use the corrected denominator index. |
| AMPE | Paper Equations 9-13 define the currency, effective interval, decay, and odds. [`AMPEFit.m` lines 65-76](../external/anderson_2023/matlab/AMPEFit.m#L65-L76) uses the inclusive discrete range: for two occurrences it is the spacing plus one. | Use the inclusive range by default and permit an explicit range input for behavioral schedules. |
| A&M | The paper describes `desirability * retention` as odds. Both [`expFit.m` lines 51-66](../external/anderson_2023/matlab/expFit.m#L51-L66) and [`powerFit.m` lines 56-71](../external/anderson_2023/matlab/powerFit.m#L56-L71) compare that value directly with a uniform draw and later average the unbounded raw value. Both scripts use $e^{-\lambda}\lambda$, the probability of exactly one Poisson event, for a revival ([`expFit.m` lines 32-35](../external/anderson_2023/matlab/expFit.m#L32-L35), [`powerFit.m` lines 41-45](../external/anderson_2023/matlab/powerFit.m#L41-L45)). | Default to odds-to-probability conversion and the probability of at least one revival. Keep the released behavior explicit for comparison. |

For PPE, the activation/logistic interpretation is also stated in the original
model literature: Walsh et al. (2018), *Cognitive Science*,
[doi:10.1111/cogs.12602](https://doi.org/10.1111/cogs.12602).

## Released parameter order

These values come from `modelParams.mat`, not from reverse-engineering rounded
Table 1 entries:

| Model | Released vector in code order |
|---|---|
| GPE | $[c,d,\alpha]=[.58161261,.61681564,.02081011]$ |
| ACT-R | $[d,\alpha]=[.79718902,.04010094]$ |
| P&A | $[c,a,\alpha]=[.45391362,.75811324,.05296011]$ |
| PPE | $[x,c,b,m,\alpha]=[8.6986,.6178,.5358,.1862,.0180]$ |
| MCM | $[\mu,\nu,\omega,\xi,\alpha]=[.0316,1.1112,.7041,.9784,.0288]$ |
| AMPE | $[\alpha,b,t_P,g_P]=[214.107907,1401.134954,15.174529,1564.979129]$ |
| Exponential A&M | $[k_\pi,\theta_\pi,\mu_d,\mu_R]=[.1637,.1390,.0346,333]$ |
| Power A&M | $[k_\pi,\theta_\pi,\mu_d,\mu_R]=[.1986,.4820,4.0757,800]$ |

The PPE value $\alpha=.018$ was fitted under Anderson's activation-as-odds
mapping. It is not a fit for the corrected exponential mapping.

## Data limitations relevant to the rewrite

- The released DOCX readme says that model fitting starts from
  `Combined.mat/counts225`, but the archive contains no function that produces
  `counts225`. The raw-event-to-`counts225` step therefore cannot be replayed
  exactly from the release.
- The deterministic CSV conversion is exact for released arrays and records
  source hashes in its manifest. Reconstructing `counts225` from raw events is
  a separate validation task, not part of conversion.
- The A&M scripts do not set or record random seeds. Their exact cached Monte
  Carlo surfaces cannot be required as bitwise targets for the Python rewrite.

No reaction-time fits, Hick-law analysis, SAM reconstruction, behavioral-study
refits, or general paper-error catalogue is part of this implementation check.

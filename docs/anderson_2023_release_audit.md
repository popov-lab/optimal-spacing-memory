# Anderson et al. (2023) public-release audit

- Source page: https://act-r.psy.cmu.edu/?p=32939&post_type=publications
- Archive: https://act-r.psy.cmu.edu/wordpress/wp-content/uploads/2021/10/paper_matlab.zip
- Author-hosted PDF: https://act-r.psy.cmu.edu/wordpress/wp-content/uploads/2023/02/EnvironmentalBassisofMemory.pdf

The current archive is 168,168,179 bytes (357,071,370 bytes unpacked), and its
server `Last-Modified` date is 2023-03-02. No GitHub, OSF, Dataverse, or other
public mirror was found. There is no license file. At the user's direction,
the author-hosted release has been imported byte-for-byte under
`external/anderson_2023/matlab/`; its provenance, checksums, omissions, and the
authors' stated redistribution terms are recorded in that directory's README.

## Release inventory

The archive is MATLAB-only: 36 `.m` files, 15 `.mat` files, and a DOCX readme.
It contains:

- processed environmental arrays (`Combined.mat`, `twitterData.mat`,
  `redditData.mat`, `FigureC1.mat`, and smaller derived arrays);
- environmental model functions for GPE, ACT-R, P&A, PPE, MCM, exponential
  A&M, power A&M, and AMPE;
- scripts for the human AMPE fits, Hick-law analysis, Appendix C analyses, and
  paper displays;
- a DOCX readme;
- no package manifest, license file, test suite, fixed environment, or
  one-command reproduction entry point;
- no SAM code.

`Combined.mat` contains the model-fitting target and a
`1000 x 1000 x 225` `int32` history-count array. The array occupies about 900
MB when loaded, contains 10,404,456 nonzero cells, and summarizes
2,384,082,074 histories. It is sufficient to evaluate the deterministic
models, but no released source constructs it. Thus the raw-to-fit pipeline
cannot be rerun end to end.

The release is a working-directory snapshot rather than a reproducible
software package: functions load files by relative name, several display and
fit variants coexist, and intermediate `.mat` outputs determine which result
is shown.

## Equation-to-code comparison

| Model | Agreement | Material difference |
|---|---|---|
| GPE | Core odds equation agrees. | Released parameter vector differs slightly from printed Table 1. |
| ACT-R | Uses \(\sum r_j^{-d}\) and a positive output scale. | Confirms that Table 1's negative scale is a sign error. Earlier ages for \(N>2\) are approximated analytically rather than represented exactly. |
| P&A | Implements presentation-specific decay. | Scale is stored as \(e^B\), but the environmental fit treats scaled strength directly as probability rather than converting odds. Its prototype-history code does not faithfully evaluate all prior components for the newest occurrence. |
| PPE | Equations A3-A4 and the \(N=1\) special case agree. | Earlier occurrences for \(N>2\) are spread over an undocumented prototype ending near the 1,000-text boundary. |
| MCM | Confirms simultaneous trace updates, \(K=100\), and the Equation 7 odds transform. | Truncates negative presentation increments with `max(0, ...)`, omitted from Equation A6. |
| AMPE | Confirms harmonic mean plus one, inclusive range, and the closed-form odds. | Confirms a singleton range of 1 and adjacent-pair range of 2, contradicting the prose equation of range with positional lag; its environmental fitter also imposes an undocumented lower bound greater than 1 on every parameter. |
| A&M | Confirms Gamma desirability, exponential item decay distribution, persistent item parameters, and inverse revival-interval parameterization. | Treats `desirability * retention` directly as probability, not odds; uses a nonstandard per-step revival probability; the exponential and power scripts implement different initial-history conventions. |

Published Equation A8 also has a dummy-index error in the MCM weight
denominator: it must sum \(\xi^j\), as the code does, rather than repeat
\(\xi^i\). The following sentence says \(\omega<1\) makes weights decrease;
the relevant condition is \(\xi<1\). The released fitter enforces \(\xi<1\),
not \(\omega<1\).

## Parameters stored in `modelParams.mat`

These are not identical to all rounded values in Table 1 and appear to come
from a nearby optimization or saved analysis state.

| Model | Released values, in code order |
|---|---|
| GPE | \(c=.58161261, d=.61681564, A=.02081011\) |
| ACT-R | \(d=.79718902, \kappa=.04010094\) |
| P&A | \(c=.45391362, a=.75811324, \kappa=.05296011\) |
| PPE | \(x=8.6986, c=.6178, b=.5358, m=.1862, A=.0180\) |
| MCM | \(\mu=.0316, \nu=1.1112, \omega=.7041, \xi=.9784, A=.0288\) |
| AMPE | \(a=214.107907, b=1401.134954, t_P=15.174529, g_P=1564.979129\) |
| Exponential A&M | \(k_\pi=.1637, \theta_\pi=.1390, \mu_d=.0346, \mu_R=333\) |
| Power A&M | \(k_\pi=.1986, \theta_\pi=.4820, \mu_d=4.0757, \mu_R=800\) |

## Critical discrepancies and errors

### 1. Odds are not handled consistently

The paper defines desirability and the deterministic outputs as odds. GPE,
ACT-R, PPE, MCM, and AMPE generally apply

\[
p=\frac{O}{1+O}.
\]

The released P&A fit instead uses \(\kappa\sum r_j^{-d_j}\) directly as a
probability. Both A&M scripts similarly generate occurrences using

```matlab
history = rand(...) < desirability .* retention
```

without an odds-to-probability conversion. Values above one are therefore
silently treated as probability one for occurrence draws, but the unbounded
raw values are retained and averaged as target predictions. The independent
Python implementation defaults to the stated odds interpretation and exposes
the released scoring rule only as an explicit comparison mode.

### 2. Revival probability is misnamed and approximated unusually

The fitted 333 and 800 are mean inter-revival intervals, not Poisson rates;
the MATLAB code takes their reciprocals. It then uses

\[
p_{\mathrm{revival}}=\lambda e^{-\lambda},
\]

the probability of exactly one Poisson event, rather than the probability of
at least one event,

\[
1-e^{-\lambda}.
\]

The difference is tiny for \(\lambda=1/333\) or \(1/800\), but the implemented
law is not the one described.

### 3. The two A&M variants do not share one simulation engine

The exponential script samples a backward pre-window revival age separately
for each item. The power script constructs elapsed times through one flattened
vector and does not implement the same explicit backward-age draw. Boundary
handling in that vector is opaque. Consequently, differences between the
exponential and power results are not attributable solely to the decay law.

The release also supplies neither random seeds nor the parameter-search code,
although the paper says it used a fixed seed during search and a new seed for
the final run. `expFit` reuses decay/revival templates across batches, whereas
`powerFit` redraws them, flattens revival processes across item boundaries, and
contains an unexplained rule forbidding every 2,999th candidate revival. Exact
simulation reproduction is therefore impossible from the release.

The reported output multiplier \(A\) is not passed into either simulator.
Instead, the scripts rescale raw conditional prediction values directly so
their mean *log* value (equivalently, geometric mean) matches the observations.
This is different from scaling odds and from the paper's statement about
matching mean probability.

### 4. AMPE range is inclusive in code

The released AMPE uses

\[
G=(\text{oldest age}-\text{newest age})+1.
\]

This resolves the examples in Figure 4 but contradicts prose saying range is
the sum of consecutive lags and equals spacing for \(N=2\). The Python code
uses the released inclusive definition and exposes an explicit `range_value`
argument for alternative behavioral encodings.

### 5. Prototype histories are part of the fitted model

For \(N>2\), the environmental fitting functions construct earlier occurrences
at approximately even intervals between the second-most-recent occurrence and
the old edge of the 1,000-text window. Exact rounding differs by model (`ceil`
in P&A; noninteger ages in PPE and AMPE). Thus the reported fits are not merely
fits of the printed equations to observed histories; they are fits of equations
plus model-specific synthetic histories.

ACT-R uses an integral approximation for unobserved earlier ages rather than
the same prototype. This contradicts Appendix B's claim that its grid contains
the full information needed for ACT-R.

### 6. P&A's prototype implementation is not the printed recursion

The core Python implementation evaluates Equation A2 exactly for an arbitrary
history. The released `PAFit.m` builds prototype schedules, reverses component
decays, and calculates the newest component using only one stored earlier
decay where Equation A2 requires the sum over every earlier component. This is
not a general implementation of the printed P&A model.

### 7. MCM has an undocumented truncation

Equation A6 permits a negative increment whenever a cumulative weighted trace
mean exceeds one. The released functions instead use

```matlab
max(0, 1 - strength)
```

and compute all increments simultaneously from the pre-update state. The
Python implementation follows the released behavior and documents it in the
function.

### 8. The readme points to stale A&M display arrays

The release contains both `display46.mat` and `display461.mat`. The readme
directs users to the former, but its saved A&M predictions do not yield the
Table 1 results. The undocumented latter file is the display state consistent
with Table 1. A reproduction that follows the readme therefore shows stale
results.

Recomputing the paper's statistics from the cached surfaces gives:

| Model/cache | RMSE | \(r^2\) | compared cells |
|---|---:|---:|---:|
| GPE | .581342 | .886138 | 513 |
| ACT-R | .824587 | .784265 | 513 |
| P&A | .683574 | .847171 | 513 |
| PPE | .548636 | .898373 | 513 |
| MCM | .578254 | .888356 | 513 |
| Exponential A&M, `display461.mat` | .573769 | .905927 | 513 |
| Power A&M, `display461.mat` | .409184 | .943892 | 481 |
| AMPE | .397781 | .946969 | 513 |

These round to Table 1. Power A&M is evaluated on only 481 cells because its
simulation masks cells with fewer than 1,000 simulated cases; the other fits
use 513. The stale `display46.mat` instead gives .654085/.901801 for
exponential A&M and .505429/.941513 for power A&M.

### 9. Human-fitting scripts retain abandoned alternatives

Exploratory scripts and saved artifacts contain alternative AMPE currency
rules, including most-recent age, arithmetic mean, and hybrids; for example,
the active rule in `fit4s.m` is not Equation 10. These artifacts should not be
mistaken for the final behavioral pipeline.

The final `Appendix.m` implementation does use the paper's harmonic-mean
currency and orders parameters as `[b,tP,gP,threshold,s]`. Its fitting bounds
include undocumented upper limits of 1,000 on both priors. `Appendix(0)` uses
saved fits deterministically; the readme instead recommends unseeded random
starts. It also hard-codes Rumelhart's first condition to chance (1/3) and
silently averages two 64-condition prediction blocks when a schedule has 128
elements. Many released schedule patterns contain age zero; MATLAB's harmonic
mean then returns zero and AMPE currency becomes one after Equation 10's
offset. This is a different time convention from environmental ages, whose
minimum is one.

### 10. Hick-law helpers contain harmless-but-real errors

The GPE and ACT-R helper function names are reversed, although their caller
routes the labels to the intended formulas. The GPE helper also applies the
odds multiplier twice. Because the latency fit has a free multiplicative
scale, that constant is absorbed and the attainable fit is unchanged, but the
helper is not the published equation.

## Missing pipeline pieces and portability problems

- `redditData.mat` lacks the precomputed summary arrays needed by the readme's
  `displayPair` call, and the referenced `summarizeReddit.mat` is absent.
- The first documented `getNs` call has an inconsistent transpose; the readme
  also repeats argument `n` where the simulators take `n,m`.
- `summarizeTwitter` and `summarizeReddit` call Wavelet Toolbox's vector
  reversal function `wrev` on matrices. The intended row reversal is
  `flipud`; the released call is version-dependent and adds an unnecessary
  toolbox dependency.
- File references differ in case (`TwitterData`/`twitterData.mat` and
  `figureC1`/`FigureC1.mat`), failing on case-sensitive systems.
- `getNs` and `freqRange` omit the final valid 2,000-text window for each
  source.
- Reproduction requires MATLAB plus Statistics and Machine Learning and
  Parallel Computing toolboxes, and as written the Wavelet Toolbox. No MATLAB
  or toolbox versions are recorded.

## Paper-level issues not repaired silently

- SAM is not specified or released.
- AMPE has no \(N=0\) rule.
- The 3-1-2 behavioral pseudo-exposure schedule is incomplete.
- Mixed-design currency and range use different clocks.
- Parameter bounds and optimizer settings are not fully reported.
- Table 2's mean P&A \(r^2=.4\) is evidently \(.94\).
- Five occurrences create four consecutive lags, not five.
- A 3,000-event series ordinarily supplies 2,000 history-plus-target windows,
  not the reported 1,999.
- Table 1 assigns \(g_P\) to Equation 12; it is defined in Equation 11.
- The text cross-references Equation 14 as AMPE odds; Equation 14 is the recall
  mapping, while AMPE odds follow from Equations 9-13.
- Equation A8 uses the wrong denominator dummy index, and the following MCM
  constraint names \(\omega\) where decreasing weights require \(\xi<1\).

## False monotonic-practice claim

The paper says adding practice within a fixed range can never lower AMPE odds,
regardless of parameterization. Holding range fixed keeps \(M\) and \(d\)
fixed, but adding an older internal occurrence can raise harmonic-mean
currency enough to dominate the multiplier \(N\).

For \(t_P=1\), ages \(\{1,100\}\), and \(d=10\):

\[
T_{\mathrm{old}}=1+\operatorname{HM}(1,100,1)=2.4925.
\]

Add an internal occurrence at age 99 without changing the endpoints:

\[
T_{\mathrm{new}}=1+\operatorname{HM}(1,99,100,1)=2.9801.
\]

Then

\[
\frac{O_{\mathrm{new}}}{O_{\mathrm{old}}}
=\frac{3}{2}\left(\frac{2.9801}{2.4925}\right)^{-10}
\approx .251.
\]

Added practice lowers predicted odds by about 75%. The claim may hold in a
restricted fitted region, but not for all positive parameterizations.

## Implementation policy used here

The Python module follows the paper's intended mathematical model where that
intent is coherent, uses released code to resolve genuine omissions (inclusive
range, MCM update order and truncation, positive ACT-R/P&A scales), and makes
repairs explicit where neither source is coherent (odds conversion and the
shared discrete convention that the first sampled event after revival has age
one). The A&M output scale is kept out of history generation and applied only
after conditional averaging. The code exposes an explicit
`released_probability` mode for comparison, but does not silently reproduce
stale arrays, model-specific binning approximations, or fitting bugs.

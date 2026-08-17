"""Numerical verification of the results in the manuscript section
"Generalization: which learning and forgetting functions produce the empirical
spacing profile?".

Model:  S(a, b) = g0 * f(a+b) + g(g0 * f(a)) * f(b),   g0 = g(0)
  g : retrieved strength -> new learning (encoding function)
  f : forgetting function, f(0) = 1, strictly decreasing
  a : inter-study interval,  b : retention interval (test delay after study 2)
Derived quantities:
  h(t)   = -f'(t)/f(t)      hazard of forgetting
  m(t)   = -f'(t)           momentary forgetting rate
  lam(t) = -m'(t)/m(t)      proportional decay rate of the forgetting rate
  c(s)   = -g'(s)           encoding-suppression slope

Checks (numbering follows the theorems in the manuscript section):
  V1  Theorem 1: massed-optimality threshold b_c = h^{-1}(c(g0) h(0))
  V2  Theorem 2: exponential f  =>  a*(b) constant in b (any g); delta rule => massed
  V3  Theorem 3/4: sign rule; power law: a* strictly increasing, trapped below a=(1+b)/gamma;
      closed form of the delta-rule solution matches brute-force argmax
  V4  Theorem 3 (converse direction): a g for which a*(b) DEcreases while the
      optimum sits above the line a=(1+b)/gamma, then increases after crossing
  V5  Theorem 5: universal log-log slope 1/(1+gamma p) and constant kappa p g0^(p-1)
  V6  Theorem 6: exponential-tail f (two-exponential mixture) => a*(b) saturates
      at the root of c(g0 f(a)) m(a) = r1 exp(-r1 a)

Figures:
  figures/f3_general_universality.png     (Theorem 5)
  figures/f4_general_forgetting_tails.png (Theorems 2, 3, 6)

Run from the repository root:  python3 src/spacing-general-verification.py
Requires: numpy, scipy, matplotlib.
"""

import numpy as np
from scipy.optimize import brentq

# ----------------------------------------------------------------------------
# brute-force global argmax over a for fixed b
# ----------------------------------------------------------------------------

def argmax_a(S, b, a_max=None, n=4001):
    """Global argmax of a -> S(a, b) on [0, a_max]: dense grid + golden-section
    refinement, with an explicit comparison against the massed boundary a=0."""
    if a_max is None:
        a_max = max(50.0, 400.0 * (1 + b))
    grid = np.concatenate([[0.0], np.logspace(-6, np.log10(a_max), n)])
    vals = S(grid, b)
    i = int(np.argmax(vals))
    if i == 0:
        return 0.0
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    phi = (np.sqrt(5) - 1) / 2
    x1, x2 = hi - phi * (hi - lo), lo + phi * (hi - lo)
    f1, f2 = S(np.array([x1]), b)[0], S(np.array([x2]), b)[0]
    for _ in range(200):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + phi * (hi - lo)
            f2 = S(np.array([x2]), b)[0]
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - phi * (hi - lo)
            f1 = S(np.array([x1]), b)[0]
        if hi - lo < 1e-12 * (1 + hi):
            break
    am = 0.5 * (lo + hi)
    if S(np.array([0.0]), b)[0] >= S(np.array([am]), b)[0]:
        return 0.0
    return am


def make_S(f, g, g0):
    return lambda a, b: g0 * f(a + b) + g(g0 * f(a)) * f(b)


f_pow = lambda gam: (lambda t: (1.0 + t) ** (-gam))
f_exp = lambda gam: (lambda t: np.exp(-gam * t))

PASS = {True: "PASS", False: "FAIL"}
results = []

def report(label, ok):
    results.append(ok)
    print(f"[{PASS[ok]}] {label}")


# ----------------------------------------------------------------------------
# V1  Theorem 1: onset threshold b_c (power law + delta rule: b_c = 1/delta - 1)
# ----------------------------------------------------------------------------
gam, delta = 0.5, 0.4
S = make_S(f_pow(gam), lambda s: delta * (1 - s), delta)
b_c = 1 / delta - 1
below = all(argmax_a(S, b) == 0.0 for b in [0.5 * b_c, 0.9 * b_c, 0.99 * b_c])
above = all(argmax_a(S, b) > 0.0 for b in [1.01 * b_c, 1.5 * b_c, 3.0 * b_c])
report("V1 massed iff b <= b_c = 1/delta - 1 (power+delta, gam=.5, delta=.4)", below and above)

# ----------------------------------------------------------------------------
# V2  Theorem 2: exponential forgetting
# ----------------------------------------------------------------------------
# (a) concave g with suppression slope crossing 1: c(s)=0.5+0.4s = 1 at s=1.25
#     => a* = ln(1/0.625)/gam, the same for every b
g_conc, g0_conc = (lambda s: 2.0 - 0.5 * s - 0.2 * s ** 2), 2.0
ok = True
for gam_e in [0.3, 1.0]:
    S = make_S(f_exp(gam_e), g_conc, g0_conc)
    a_pred = np.log(1 / 0.625) / gam_e
    stars = [argmax_a(S, b, a_max=200.0) for b in [0.1, 1.0, 10.0, 100.0]]
    ok &= max(abs(x - a_pred) for x in stars) < 1e-4 * (1 + a_pred)
report("V2a exponential f: a*(b) constant and equal to f^{-1}(s*/g0), any b", ok)
# (b) delta rule (c = delta < 1): massed for every b
S = make_S(f_exp(0.5), lambda s: 0.7 * (1 - s), 0.7)
ok = all(argmax_a(S, b, a_max=100.0) == 0.0 for b in [0.1, 1.0, 10.0])
report("V2b exponential f + delta rule: massed for every b", ok)

# ----------------------------------------------------------------------------
# V3  Theorems 3/4: power law closed form, monotonicity, trapping region
# ----------------------------------------------------------------------------
gam, delta = 0.5, 0.4
S = make_S(f_pow(gam), lambda s: delta * (1 - s), delta)

def a_star_closed(b, delta, gam):
    xt = 1.0 + b
    alpha = delta ** (-1 / (gam + 1))
    return max((xt - 1) / (alpha * xt ** (gam / (1 + gam)) - 1) - 1, 0.0)

ok_match = ok_inc = ok_trap = True
prev = -1.0
for b in np.logspace(np.log10(1.6), 6, 25):
    a_num = argmax_a(S, b)
    ok_match &= abs(a_num - a_star_closed(b, delta, gam)) < 1e-3 * (1 + a_num)
    ok_inc &= a_num > prev
    ok_trap &= a_num < (1 + b) / gam
    prev = a_num
report("V3a power+delta: brute-force argmax matches the closed form", ok_match)
report("V3b power+delta: a*(b) strictly increasing for b > b_c", ok_inc)
report("V3c power+delta: a*(b) < (1+b)/gamma (invariant region)", ok_trap)

# ----------------------------------------------------------------------------
# V4  Theorem 3 converse: a decreasing stretch of a*(b) above the line
# ----------------------------------------------------------------------------
gam = 1.0
c0_, rho = 1.9, 1.5     # g(s) = 1 - (c0_/rho)(1 - e^{-rho s}):  c(s) = c0_ e^{-rho s}
g_cex = lambda s: 1.0 - (c0_ / rho) * (1.0 - np.exp(-rho * s))
S = make_S(f_pow(gam), g_cex, 1.0)
bs_cex = np.linspace(0.2, 8.0, 40)
a_cex = np.array([argmax_a(S, b, a_max=1e4) for b in bs_cex])
interior = a_cex > 0
dec_above = inc_below = False
for i in range(len(bs_cex) - 1):
    if interior[i] and interior[i + 1]:
        above = a_cex[i] > (1 + bs_cex[i]) / gam
        if above and a_cex[i + 1] < a_cex[i] - 1e-6:
            dec_above = True
        if (not above) and a_cex[i + 1] > a_cex[i] + 1e-6:
            inc_below = True
report("V4  power f, plateau g: a* decreases above a=(1+b)/gam, increases below", dec_above and inc_below)

# ----------------------------------------------------------------------------
# V5  Theorem 5: universality  (1+a*)^(1+gam p)/(1+b) -> kappa p g0^(p-1)
# ----------------------------------------------------------------------------
gam = 0.4
fam = {
    # label:                      (g,                                        g0,  p,   K = kappa p g0^(p-1))
    "delta 0.5      (p=1)":       (lambda s: 0.5 * (1 - s),                  0.5, 1.0, 0.5),
    "delta 0.9      (p=1)":       (lambda s: 0.9 * (1 - s),                  0.9, 1.0, 0.9),
    "convex exp     (p=1)":       (lambda s: 0.8 * np.exp(-0.625 * s),       0.8, 1.0, 0.5),
    "sqrt, kappa=.8 (p=1/2)":     (lambda s: 1.0 - 0.8 * np.sqrt(np.maximum(s, 0)), 1.0, 0.5, 0.4),
    "quadr, kappa=.9(p=2)":       (lambda s: 0.9 * (1 - s ** 2),             0.9, 2.0, 1.62),
}
ok_slope = ok_const = True
slopes = {}
for name, (g, g0, p, K) in fam.items():
    S = make_S(f_pow(gam), g, g0)
    bs = np.logspace(3, 7, 9)
    la = np.array([np.log(1 + argmax_a(S, b, a_max=50 * (1 + b))) for b in bs])
    slope = np.polyfit(np.log(1 + bs), la, 1)[0]
    slopes[name] = slope
    ok_slope &= abs(slope - 1 / (1 + gam * p)) < 0.02
    a8 = argmax_a(S, 1e8, a_max=50 * (1 + 1e8))
    ok_const &= abs((1 + a8) ** (1 + gam * p) / (1 + 1e8) - K) < 0.15 * K
report("V5a fitted log-log slopes within .02 of 1/(1+gamma p) for all five g", ok_slope)
report("V5b (1+a*)^(1+gam p)/(1+b) -> kappa p g0^(p-1) at b = 1e8", ok_const)

# ----------------------------------------------------------------------------
# V6  Theorem 6: exponential-tail mixture => a*(b) saturates at the tail root
# ----------------------------------------------------------------------------
w, r1, r2 = 0.5, 0.1, 2.0
f_mix = lambda t: w * np.exp(-r1 * t) + (1 - w) * np.exp(-r2 * t)
m_mix = lambda t: w * r1 * np.exp(-r1 * t) + (1 - w) * r2 * np.exp(-r2 * t)
c_conc = lambda s: 0.5 + 0.4 * s
S = make_S(f_mix, g_conc, g0_conc)
a_inf = brentq(lambda a: c_conc(g0_conc * f_mix(a)) * m_mix(a) - r1 * np.exp(-r1 * a), 1e-9, 3000)
stars = [argmax_a(S, b, a_max=3000.0) for b in [50.0, 300.0, 1000.0]]
ok = max(abs(x - a_inf) for x in stars) < 1e-3 * (1 + a_inf)
report(f"V6  two-exponential mixture: a*(b) -> {a_inf:.4f} solving c(s_a) m(a) = r1 e^(-r1 a)", ok)

print()
print(f"{sum(results)}/{len(results)} checks passed")

# ============================================================================
# Figures
# ============================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#1f1f1e", "#6e6d66", "#e7e6e2"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11,
})

# ---- Figure 3: universality of the log-log law (Theorem 5) ------------------
gam = 0.4
bs = np.logspace(0, 8, 33)
lb = np.log10(1 + bs)

fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.1), constrained_layout=True)

panelA = [("delta rule, $\\delta=0.5$  ($c_0=0.5$)",            fam["delta 0.5      (p=1)"], BLUE,   "-"),
          ("delta rule, $\\delta=0.9$  ($c_0=0.9$)",            fam["delta 0.9      (p=1)"], ORANGE, "-"),
          ("convex $g=0.8e^{-0.625s}$  ($c_0=0.5$)",            fam["convex exp     (p=1)"], AQUA,   (0, (4, 2)))]
ax = axes[0]
for label, (g, g0, p, K), col, ls in panelA:
    la = np.array([np.log10(1 + argmax_a(make_S(f_pow(gam), g, g0), b, a_max=50 * (1 + b))) for b in bs])
    ax.plot(lb, la, color=col, lw=2, ls=ls, label=label)
x0 = np.array([4.8, 7.6])
ax.plot(x0, -0.35 + x0 / (1 + gam), color=INK2, lw=1.2, ls=(0, (4, 3)))
ax.annotate(f"slope $1/(1+\\gamma)={1/(1+gam):.3f}$", (x0.mean() + 0.35, -0.42 + x0.mean() / (1 + gam)),
            ha="center", va="top", color=INK2, fontsize=9, rotation=33)
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.set_title("Smooth $g$: slope universal, intercept set by $c_0=-g'(0)$")
ax.set_xlabel("$\\log_{10}(1+b)$   (retention interval)")
ax.set_ylabel("$\\log_{10}(1+a^*)$   (optimal lag)")

panelB = [("$p=1/2$: $g=1-0.8\\sqrt{s}$,  slope $\\to 5/6$",  fam["sqrt, kappa=.8 (p=1/2)"], YELLOW),
          ("$p=1$: $g=0.9(1-s)$,  slope $\\to 5/7$",           fam["delta 0.9      (p=1)"],   ORANGE),
          ("$p=2$: $g=0.9(1-s^2)$,  slope $\\to 5/9$",         fam["quadr, kappa=.9(p=2)"],   MAGENTA)]
ax = axes[1]
for label, (g, g0, p, K), col in panelB:
    la = np.array([np.log10(1 + argmax_a(make_S(f_pow(gam), g, g0), b, a_max=50 * (1 + b))) for b in bs])
    ax.plot(lb, la, color=col, lw=2, label=label)
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.set_title("$g(s)=g_0-\\kappa s^p+o(s^p)$: slope $1/(1+\\gamma p)$")
ax.set_xlabel("$\\log_{10}(1+b)$   (retention interval)")
ax.set_ylabel("$\\log_{10}(1+a^*)$")
fig.suptitle(f"Power-law forgetting ($\\gamma={gam}$): the slope of the optimal-lag law is set by $f$"
             " and by the local exponent of $g$ at zero strength", fontsize=11)
fig.savefig("figures/f3_general_universality.png", dpi=200)
plt.close(fig)

# ---- Figure 4: forgetting tails decide the fate of a*(b) (Thms 2, 3, 6) -----
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.1), constrained_layout=True)

# Panel A: same concave g, three forgetting functions
ax = axes[0]
bs4 = np.logspace(-1, 3, 41)
curves = [
    ("power $f=(1+t)^{-1/2}$:   $h\\downarrow 0$,  $a^*\\to\\infty$", make_S(f_pow(0.5), g_conc, g0_conc), AQUA),
    ("mixture $0.5e^{-t/10}\\!+0.5e^{-2t}$:   $h\\to r_1$,  $a^*$ saturates", make_S(f_mix, g_conc, g0_conc), ORANGE),
    ("exponential $f=e^{-0.15t}$:   $h$ constant,  $a^*$ constant", make_S(f_exp(0.15), g_conc, g0_conc), BLUE),
]
for label, Sf, col in curves:
    a_v = np.array([argmax_a(Sf, b, a_max=3000.0) for b in bs4])
    ax.plot(np.log10(bs4), np.log10(np.maximum(a_v, 1e-3)), color=col, lw=2, label=label)
ax.legend(loc="upper left", frameon=False, fontsize=8.5)
ax.set_ylim(-0.35, 2.6)
ax.set_title("Same learning rule $g$, three forgetting tails")
ax.set_xlabel("$\\log_{10} b$   (retention interval)")
ax.set_ylabel("$\\log_{10} a^*(b)$   (optimal lag)")

# Panel B: the sign rule at work (counterexample trajectory vs the line)
ax = axes[1]
bline = np.linspace(0.0, 8.0, 100)
ax.fill_between(bline, (1 + bline) / 1.0, 12.0, color=GRID, alpha=0.55, lw=0)
ax.plot(bline, (1 + bline) / 1.0, color=INK2, lw=1.2, ls=(0, (4, 3)))
ax.annotate("$a=(1+b)/\\gamma$\nabove it $da^*/db<0$\nbelow it $da^*/db>0$",
            (2.0, 8.6), color=INK2, fontsize=9, ha="left")
mask0 = a_cex == 0
ax.plot(bs_cex[~mask0], a_cex[~mask0], color=BLUE, lw=2)
ax.plot(bs_cex[mask0], a_cex[mask0], color=BLUE, lw=2)
ax.plot([bs_cex[mask0][-1], bs_cex[~mask0][0]], [0, a_cex[~mask0][0]],
        color=BLUE, lw=1.2, ls=(0, (2, 2)))
ax.annotate("massed", (bs_cex[mask0].mean(), 0.15), color=INK, fontsize=9, ha="center")
ax.annotate("onset jump", (bs_cex[~mask0][0] + 0.12, a_cex[~mask0][0] * 0.5),
            color=INK, fontsize=9, ha="left")
ax.set_ylim(-0.4, 12)
ax.set_title("Power $f$ ($\\gamma=1$), plateau $g$: non-monotone $a^*(b)$")
ax.set_xlabel("$b$   (retention interval)")
ax.set_ylabel("$a^*(b)$")
fig.suptitle("The retention-interval dependence of the optimal lag is decided by the forgetting function alone",
             fontsize=11)
fig.savefig("figures/f4_general_forgetting_tails.png", dpi=200)
plt.close(fig)

print("wrote figures/f3_general_universality.png, figures/f4_general_forgetting_tails.png")

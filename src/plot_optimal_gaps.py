"""Reproduce the meta-analysis figure of Mozer et al. (2009): optimal gap by
test delay on log-log axes, with equal decade spacing on both axes (square
plot box). Data points only, plus the least-squares power-law fit.

Usage: python src/plot_optimal_gaps.py
Reads  data/cepeda2006_optimal_gaps.csv, writes figures/optimal_gap_by_test_delay.png
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent

delay, gap = [], []
with open(ROOT / "data" / "cepeda2006_optimal_gaps.csv") as f:
    for row in csv.DictReader(f):
        delay.append(float(row["test_delay_days"]))
        gap.append(float(row["optimal_gap_days"]))
delay, gap = np.array(delay), np.array(gap)

# power-law fit (OLS in log10 space)
b, a = np.polyfit(np.log10(delay), np.log10(gap), 1)
xfit = np.array([delay.min(), delay.max()])
yfit = 10 ** (a + b * np.log10(xfit))

fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.plot(xfit, yfit, "--", color="#333333", lw=1.5, zorder=1)
ax.plot(delay, gap, "o", mfc="none", mec="#D62728", mew=1.4, ms=7, zorder=2)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-5, 1e3)
ax.set_ylim(1e-6, 1e2)
ax.set_aspect(1)  # one decade spans the same distance on both axes -> square box
ax.set_xlabel("Test Delay (days)")
ax.set_ylabel("Optimal Gap (days)")
ax.grid(True, which="major", color="#DDDDDD", lw=0.6)
ax.tick_params(which="minor", length=0)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)

fig.tight_layout()
out = ROOT / "figures" / "optimal_gap_by_test_delay.png"
fig.savefig(out, dpi=300)
print(f"wrote {out}  (fit: gap = {10**a:.3f} * delay^{b:.3f})")

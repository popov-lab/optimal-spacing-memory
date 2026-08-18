"""Plot the recall data of Cepeda et al. (2008, 2009). Panels a-c: forgetting
and spacing functions of Cepeda et al. (2009) Experiments 1, 2a, and 2b, read
from Figures 3 and 4 of that paper. Panel d: spacing functions of
Cepeda et al. (2008) for the four retention intervals.

Each panel is produced twice: with a linear x axis (as in the source figures)
and with logarithmic x spacing. Because the massed condition (ISI = 0) has no
place on a log axis, the log versions use a symlog axis: linear below 1 day,
logarithmic above.

Usage: python src/plot_spacing_recall.py
Reads  data/cepeda_spacing_recall.csv,
writes figures/spacing_recall_{a,b,c,d}.png and figures/spacing_recall_{a,b,c,d}_log.png
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent

rows = list(csv.DictReader(open(ROOT / "data" / "cepeda_spacing_recall.csv")))
for r in rows:
    r["isi_days"] = float(r["isi_days"])
    r["recall_pct"] = float(r["recall_pct"])

panels = {
    "a": dict(title="(a) RI = 10 days", xticks=[0, 1, 2, 4, 7, 14],
              xticks_linear=[0, 1, 2, 4, 7, 14]),
    "b": dict(title="(b) RI = 168 days", xticks=[0, 1, 7, 28, 84, 168],
              xticks_linear=[0, 7, 28, 84, 168]),  # "1" would overlap "0"
    "c": dict(title="(c) RI = 168 days", xticks=[0, 1, 7, 28, 84, 168],
              xticks_linear=[0, 7, 28, 84, 168]),
    "d": dict(title="(d) RIs = 7, 35, 70, 350 days",
              xticks=[0, 1, 7, 14, 21, 35, 70, 105],
              xticks_linear=[0, 7, 14, 21, 35, 70, 105],
              # unlabeled marks at 14, 35, 70: their labels would collide in log spacing
              xlabels_log=[0, 1, 7, 21, 105]),
}
OLIVE = "#666600"
RI_COLORS = {7: "#FF0000", 35: "#A95400", 70: "#54A900", 350: "#00FF00"}

for p, cfg in panels.items():
    for scale in ("linear", "log"):
        fig, ax = plt.subplots(figsize=(4.2, 4.0))
        prows = [r for r in rows if r["panel"] == p]
        if p == "d":
            by_ri = defaultdict(list)
            for r in prows:
                by_ri[int(r["ri_days"])].append((r["isi_days"], r["recall_pct"]))
            for ri, pts in sorted(by_ri.items()):
                x, y = zip(*sorted(pts))
                ax.plot(x, y, "o:", color=RI_COLORS[ri], mfc=RI_COLORS[ri],
                        ms=6, lw=1, label=f"RI = {ri} d")
            ax.legend(frameon=False, fontsize=9, ncol=2, loc="lower center",
                      bbox_to_anchor=(0.55, 0.02), columnspacing=1.2)
        else:
            for func, style in [("forgetting", dict(marker="s", ls=":", color="#0000FF",
                                                    mfc="none", label="forgetting function")),
                                ("spacing", dict(marker="o", ls=":", color=OLIVE,
                                                 mfc=OLIVE, label="spacing function"))]:
                pts = sorted((r["isi_days"], r["recall_pct"])
                             for r in prows if r["function"] == func)
                x, y = zip(*pts)
                ax.plot(x, y, ms=6, lw=1, **style)
            ax.legend(frameon=False, fontsize=9)
        xmax = cfg["xticks"][-1]
        if scale == "log":
            ax.set_xscale("symlog", linthresh=1)
            ax.set_xlim(-0.12, xmax * 1.4)
            labeled = cfg.get("xlabels_log", cfg["xticks"])
            ax.set_xticks(cfg["xticks"])
            ax.set_xticklabels([f"{t:g}" if t in labeled else "" for t in cfg["xticks"]])
            ax.tick_params(which="minor", length=0)
        else:
            ax.set_xlim(-0.025 * xmax, 1.025 * xmax)
            ax.set_xticks(cfg["xticks_linear"])
        ax.set_ylim(0, 100)
        ax.set_xlabel("ISI (days)")
        ax.set_ylabel("% recall")
        ax.set_title(cfg["title"], fontsize=11)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        fig.tight_layout()
        suffix = "_log" if scale == "log" else ""
        out = ROOT / "figures" / f"spacing_recall_{p}{suffix}.png"
        fig.savefig(out, dpi=300)
        plt.close(fig)
        print("wrote", out)

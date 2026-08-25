"""Plot the four Cepeda et al. (2008) MCM history-ablation fits.

This script uses the post-hoc parameter fits saved by ``ablate_mcm_histories.py``
and overlays each fitted model on the 26 published spacing observations.

Usage:
    python src/plot_mcm_history_ablation_curves.py

Outputs:
    figures/mcm_2008_full_stochastic.svg
    figures/mcm_2008_no_retrieval_branching.svg
    figures/mcm_2008_no_encoding_branching.svg
    figures/mcm_2008_no_branching.svg
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ablate_mcm_histories import VARIANTS, ZERO_GAP_2008, two_session_variant
from mcm import MCMParams

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "cepeda_spacing_recall.csv"
PARAMS = ROOT / "results" / "mcm_history_ablation_2008_params.csv"
FIGURES = ROOT / "figures"

SLUGS = {
    "full stochastic": "full_stochastic",
    "no retrieval branching": "no_retrieval_branching",
    "no encoding branching": "no_encoding_branching",
    "no branching": "no_branching",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def load_params() -> dict[str, tuple[MCMParams, float]]:
    out: dict[str, tuple[MCMParams, float]] = {}
    for row in read_csv(PARAMS):
        out[row["variant"]] = (
            MCMParams(
                mu=float(row["mu"]),
                nu=float(row["nu"]),
                omega=float(row["omega"]),
                xi=float(row["xi"]),
                epsilon_r=float(row["epsilon_r"]),
                n=100,
            ),
            float(row["rmse_pp"]),
        )
    return out


def main() -> None:
    rows = [
        row
        for row in read_csv(DATA)
        if row["panel"] == "d" and row["function"] == "spacing"
    ]
    fitted = load_params()
    FIGURES.mkdir(parents=True, exist_ok=True)

    for variant, (retrieval_mode, encoding_mode) in VARIANTS.items():
        params, rmse = fitted[variant]
        fig, ax = plt.subplots(figsize=(7.3, 5.1))

        for ri in (7, 35, 70, 350):
            points = sorted(
                [
                    (float(row["isi_days"]), float(row["recall_pct"]))
                    for row in rows
                    if float(row["ri_days"]) == ri
                ]
            )
            x_obs = np.array([x for x, _ in points])
            y_obs = np.array([y for _, y in points])

            x_grid = np.r_[np.linspace(0, 1, 120), np.geomspace(1.01, 105, 400)]
            y_grid = np.array(
                [
                    100.0
                    * two_session_variant(
                        ZERO_GAP_2008 if x == 0 else x,
                        ri,
                        params,
                        retrieval_mode,
                        encoding_mode,
                    )
                    for x in x_grid
                ]
            )

            line, = ax.plot(x_grid, y_grid, label=f"RI = {ri} d")
            ax.plot(x_obs, y_obs, "o", color=line.get_color())

        ax.set_xscale("symlog", linthresh=1)
        ax.set_xticks([0, 1, 7, 21, 105])
        ax.set_xticklabels(["0", "1", "7", "21", "105"])
        ax.set_xlabel("ISI (days)")
        ax.set_ylabel("Final-test recall (%)")
        ax.set_ylim(0, 100)
        ax.set_title(f"Cepeda et al. (2008): {variant.title()} (RMSE {rmse:.2f} pp)")
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout()

        output = FIGURES / f"mcm_2008_{SLUGS[variant]}.svg"
        fig.savefig(output)
        plt.close(fig)
        print("wrote", output)


if __name__ == "__main__":
    main()

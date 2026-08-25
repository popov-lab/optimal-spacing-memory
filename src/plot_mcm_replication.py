"""Plot the MCM replication analyses for Cepeda et al. (2008, 2009).

Requires the result files produced by ``fit_mcm_cepeda2009.py`` and
``fit_mcm_cepeda2008.py``. Writes individual observed-vs-model figures for the
three 2009 forgetting and spacing functions and two views of the 2008 joint
post-hoc spacing fit.

Usage
-----
python src/plot_mcm_replication.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mcm import MCMParams, single_study_recall, two_session_spacing

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "cepeda_spacing_recall.csv"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
ZERO_GAP_2008 = 0.00256


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def params_from_row(row: dict[str, str]) -> MCMParams:
    return MCMParams(
        mu=float(row["mu"]),
        nu=float(row["nu"]),
        omega=float(row["omega"]),
        xi=float(row["xi"]),
        epsilon_r=float(row["epsilon_r"]),
        n=int(row["n"]),
    )


def plot_2009(data: list[dict[str, str]]) -> None:
    fit_rows = {
        r["panel"]: r for r in read_csv(RESULTS / "mcm_cepeda2009_fits.csv")
    }
    titles = {
        "a": "Cepeda et al. (2009), Experiment 1",
        "b": "Cepeda et al. (2009), Experiment 2a (facts)",
        "c": "Cepeda et al. (2009), Experiment 2b (object names)",
    }

    for panel in ("a", "b", "c"):
        params = params_from_row(fit_rows[panel])
        prows = [r for r in data if r["panel"] == panel]
        forgetting = sorted(
            (r for r in prows if r["function"] == "forgetting"),
            key=lambda r: float(r["isi_days"]),
        )
        spacing = sorted(
            (r for r in prows if r["function"] == "spacing"),
            key=lambda r: float(r["isi_days"]),
        )
        xf = np.array([float(r["isi_days"]) for r in forgetting])
        yf = np.array([float(r["recall_pct"]) for r in forgetting])
        xs = np.array([float(r["isi_days"]) for r in spacing])
        ys = np.array([float(r["recall_pct"]) for r in spacing])
        ri = float(spacing[0]["ri_days"])
        grid = np.linspace(0, max(xf.max(), xs.max()), 500)

        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        ax.plot(grid, 100 * single_study_recall(grid, params), label="MCM fit")
        ax.plot(xf, yf, "o", label="Observed forgetting")
        ax.set(
            xlabel="ISI (days)",
            ylabel="Recall (%)",
            title=f"{titles[panel]}: forgetting function",
            ylim=(0, 100),
            xlim=(0, None),
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIGURES / f"mcm_forgetting_panel_{panel}.png", dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        pred = np.array([two_session_spacing(x, ri, params) for x in grid])
        ax.plot(grid, 100 * pred, label="MCM prediction")
        ax.plot(xs, ys, "o", label="Observed spacing")
        ax.set(
            xlabel="ISI (days)",
            ylabel="Final-test recall (%)",
            title=f"{titles[panel]}: spacing function (RI = {ri:g} days)",
            ylim=(0, 100),
            xlim=(0, None),
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIGURES / f"mcm_spacing_panel_{panel}.png", dpi=200)
        plt.close(fig)


def plot_2008(data: list[dict[str, str]]) -> None:
    fit_row = read_csv(RESULTS / "mcm_cepeda2008_fit.csv")[0]
    params = params_from_row(fit_row)
    rows = [
        r for r in data if r["panel"] == "d" and r["function"] == "spacing"
    ]
    ris = [7, 35, 70, 350]

    for suffix, symlog in (("linear", False), ("symlog", True)):
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        for ri in ris:
            rrows = sorted(
                (r for r in rows if float(r["ri_days"]) == ri),
                key=lambda r: float(r["isi_days"]),
            )
            xobs = np.array([float(r["isi_days"]) for r in rrows])
            yobs = np.array([float(r["recall_pct"]) for r in rrows])
            if symlog:
                grid = np.r_[np.linspace(0, 1, 120), np.geomspace(1.01, 105, 380)]
            else:
                grid = np.linspace(0, 105, 500)
            pred = np.array(
                [
                    two_session_spacing(ZERO_GAP_2008 if x == 0 else x, ri, params)
                    for x in grid
                ]
            )
            line, = ax.plot(grid, 100 * pred, label=f"MCM fit, RI = {ri} d")
            ax.plot(
                xobs,
                yobs,
                "o",
                color=line.get_color(),
                label=f"Data, RI = {ri} d",
            )

        if symlog:
            ax.set_xscale("symlog", linthresh=1)
            ax.set_xticks([0, 1, 2, 4, 7, 11, 14, 21, 35, 70, 105])
        else:
            ax.set_xlim(0, 105)
        ax.set(
            xlabel="ISI (days)",
            ylabel="Final-test recall (%)",
            title="Cepeda et al. (2008): direct joint fit of MCM to spacing data",
            ylim=(0, 100),
        )
        ax.legend(ncol=2, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIGURES / f"mcm_2008_direct_fit_{suffix}.png", dpi=220)
        plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    data = read_csv(DATA)
    plot_2009(data)
    plot_2008(data)


if __name__ == "__main__":
    main()

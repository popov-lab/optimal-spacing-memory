"""Post-hoc joint fit of MCM to Cepeda et al. (2008) spacing curves.

Mozer et al. (2009) fit MCM's four primitive parameters {mu, nu, omega, xi}
using the single-session forgetting function and then predicted the four spacing
curves. The published Cepeda et al. (2008) article and recovered dataset do not
contain that forgetting function, so the original parameter-free procedure
cannot be reconstructed from the available data.

This script therefore fits ONE common MCM parameterization jointly to all 26
published spacing observations across RI = 7, 35, 70, and 350 days. It is a
post-hoc fit and must not be interpreted as the same predictive test as the
Cepeda et al. (2009) analyses in ``fit_mcm_cepeda2009.py``.

Usage
-----
python src/fit_mcm_cepeda2008.py
python src/fit_mcm_cepeda2008.py --starts 128 --seed 20260825
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit

from mcm import MCMParams, two_session_spacing

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA = ROOT / "data" / "cepeda_spacing_recall.csv"
DEFAULT_RESULTS = ROOT / "results"

# Cepeda et al. (2008) report that the nominal zero-day gap was about 3 min.
# Their text gives an actual value around 0.00256 d. We use that value here.
# This choice has negligible impact on day-to-year predictions but avoids
# treating the two sessions as simultaneous.
ZERO_GAP_DAYS = 0.00256


def unpack(z: np.ndarray) -> MCMParams:
    return MCMParams(
        mu=float(np.exp(z[0])),
        nu=float(1.0 + np.exp(z[1])),
        omega=float(expit(z[2])),
        xi=float(expit(z[3])),
        epsilon_r=9.0,
        n=100,
    )


def pack(params: MCMParams) -> np.ndarray:
    return np.array(
        [
            np.log(params.mu),
            np.log(params.nu - 1.0),
            logit(params.omega),
            logit(params.xi),
        ],
        dtype=float,
    )


def predict(rows: list[dict[str, str]], params: MCMParams) -> np.ndarray:
    values = []
    for row in rows:
        isi = float(row["isi_days"])
        ri = float(row["ri_days"])
        model_isi = ZERO_GAP_DAYS if isi == 0 else isi
        values.append(two_session_spacing(model_isi, ri, params))
    return np.asarray(values)


def fit_spacing(
    rows: list[dict[str, str]],
    recall: np.ndarray,
    starts: int,
    seed: int,
) -> tuple[MCMParams, float]:
    """Multistart least-squares fit to all four RI curves jointly."""

    rng = np.random.default_rng(seed)
    best_sse = np.inf
    best_params: MCMParams | None = None

    def residual(z: np.ndarray) -> np.ndarray:
        return predict(rows, unpack(z)) - recall

    initials = [MCMParams(mu=10.0, nu=1.2, omega=0.8, xi=0.9)]
    for _ in range(max(0, starts - 1)):
        initials.append(
            MCMParams(
                mu=10.0 ** rng.uniform(-4.0, 3.0),
                nu=1.0 + 10.0 ** rng.uniform(-3.0, 1.5),
                omega=rng.uniform(0.45, 0.999),
                xi=rng.uniform(0.05, 0.999),
            )
        )

    for initial in initials:
        result = least_squares(
            residual,
            pack(initial),
            max_nfev=5000,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        sse = float(np.sum(result.fun**2))
        if np.isfinite(sse) and sse < best_sse:
            best_sse = sse
            best_params = unpack(result.x)

    if best_params is None:
        raise RuntimeError("all optimization starts failed")
    return best_params, best_sse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--starts", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()

    with args.data.open(newline="") as f:
        all_rows = list(csv.DictReader(f))

    rows = [r for r in all_rows if r["panel"] == "d" and r["function"] == "spacing"]
    rows.sort(key=lambda r: (float(r["ri_days"]), float(r["isi_days"])))
    observed = np.array([float(r["recall_pct"]) / 100.0 for r in rows])

    params, spacing_sse = fit_spacing(rows, observed, args.starts, args.seed)
    predicted = predict(rows, params)
    residuals = predicted - observed
    rmse = float(np.sqrt(np.mean(residuals**2)))

    args.output_dir.mkdir(parents=True, exist_ok=True)

    fit_path = args.output_dir / "mcm_cepeda2008_fit.csv"
    with fit_path.open("w", newline="") as f:
        fields = [
            "experiment", "fit_target", "mu", "nu", "omega", "xi",
            "epsilon_r", "n", "spacing_sse", "spacing_rmse"
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "experiment": "Cepeda et al. (2008)",
                "fit_target": "26 spacing observations jointly across RI = 7, 35, 70, 350 days",
                "mu": f"{params.mu:.12g}",
                "nu": f"{params.nu:.12g}",
                "omega": f"{params.omega:.12g}",
                "xi": f"{params.xi:.12g}",
                "epsilon_r": f"{params.epsilon_r:g}",
                "n": params.n,
                "spacing_sse": f"{spacing_sse:.12g}",
                "spacing_rmse": f"{rmse:.12g}",
            }
        )

    pred_path = args.output_dir / "mcm_cepeda2008_predictions.csv"
    with pred_path.open("w", newline="") as f:
        fields = ["isi_days", "ri_days", "observed_recall", "predicted_recall", "residual"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row, obs, pred in zip(rows, observed, predicted):
            writer.writerow(
                {
                    "isi_days": row["isi_days"],
                    "ri_days": row["ri_days"],
                    "observed_recall": f"{obs:.6f}",
                    "predicted_recall": f"{pred:.6f}",
                    "residual": f"{pred - obs:.6f}",
                }
            )

    print(
        f"mu={params.mu:.8g}, nu={params.nu:.8g}, omega={params.omega:.8g}, "
        f"xi={params.xi:.8g}; spacing RMSE={100*rmse:.3f} pp"
    )
    print("wrote", fit_path)
    print("wrote", pred_path)


if __name__ == "__main__":
    main()

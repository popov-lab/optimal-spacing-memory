"""Fit MCM to the Cepeda et al. (2009) forgetting curves and predict spacing.

The four primitive parameters {mu, nu, omega, xi} are estimated *only* from
single-session forgetting data. The retrieval-success learning rate is fixed to
9 and N to 100, as in Mozer et al. (2009). The fitted parameters are then frozen
and used to predict the final-test spacing function.

Usage
-----
python src/fit_mcm_cepeda2009.py
python src/fit_mcm_cepeda2009.py --starts 256 --seed 20260825
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit

from mcm import MCMParams, single_study_recall, two_session_spacing

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATA = ROOT / "data" / "cepeda_spacing_recall.csv"
DEFAULT_RESULTS = ROOT / "results"


def unpack(z: np.ndarray) -> MCMParams:
    """Transform unconstrained optimizer coordinates into MCM parameters."""

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


def fit_forgetting(
    lags: np.ndarray,
    recall: np.ndarray,
    starts: int,
    seed: int,
) -> tuple[MCMParams, float]:
    """Multistart least-squares fit to a single forgetting function."""

    rng = np.random.default_rng(seed)
    best_sse = np.inf
    best_params: MCMParams | None = None

    def residual(z: np.ndarray) -> np.ndarray:
        return single_study_recall(lags, unpack(z)) - recall

    for _ in range(starts):
        initial = MCMParams(
            mu=10.0 ** rng.uniform(-6.0, 3.0),
            nu=1.0 + 10.0 ** rng.uniform(-3.0, 2.0),
            omega=rng.uniform(0.70, 0.999),
            xi=rng.uniform(0.02, 0.999),
        )
        result = least_squares(
            residual,
            pack(initial),
            max_nfev=4000,
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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--starts", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()

    rows = read_rows(args.data)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fit_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []

    for offset, panel in enumerate(("a", "b", "c")):
        forgetting = [r for r in rows if r["panel"] == panel and r["function"] == "forgetting"]
        spacing = [r for r in rows if r["panel"] == panel and r["function"] == "spacing"]
        forgetting.sort(key=lambda r: float(r["isi_days"]))
        spacing.sort(key=lambda r: float(r["isi_days"]))

        lags = np.array([float(r["isi_days"]) for r in forgetting])
        y = np.array([float(r["recall_pct"]) / 100.0 for r in forgetting])
        params, forgetting_sse = fit_forgetting(
            lags,
            y,
            starts=args.starts,
            seed=args.seed + offset,
        )

        spacing_errors = []
        for row in spacing:
            isi = float(row["isi_days"])
            ri = float(row["ri_days"])
            observed = float(row["recall_pct"]) / 100.0
            predicted = two_session_spacing(isi, ri, params)
            spacing_errors.append(predicted - observed)
            prediction_rows.append(
                {
                    "panel": panel,
                    "experiment": row["experiment"],
                    "material": row["material"],
                    "isi_days": f"{isi:g}",
                    "ri_days": f"{ri:g}",
                    "observed_recall": f"{observed:.6f}",
                    "predicted_recall": f"{predicted:.6f}",
                    "residual": f"{predicted - observed:.6f}",
                }
            )

        rmse = float(np.sqrt(np.mean(np.square(spacing_errors))))
        fit_rows.append(
            {
                "panel": panel,
                "experiment": forgetting[0]["experiment"],
                "material": forgetting[0]["material"],
                "mu": f"{params.mu:.12g}",
                "nu": f"{params.nu:.12g}",
                "omega": f"{params.omega:.12g}",
                "xi": f"{params.xi:.12g}",
                "epsilon_r": f"{params.epsilon_r:g}",
                "n": params.n,
                "forgetting_sse": f"{forgetting_sse:.12g}",
                "spacing_rmse": f"{rmse:.12g}",
            }
        )

        print(
            f"panel {panel}: mu={params.mu:.6g}, nu={params.nu:.6g}, "
            f"omega={params.omega:.6g}, xi={params.xi:.6g}; "
            f"forgetting SSE={forgetting_sse:.6g}; spacing RMSE={100*rmse:.2f} pp"
        )

    fits_path = args.output_dir / "mcm_cepeda2009_fits.csv"
    preds_path = args.output_dir / "mcm_cepeda2009_predictions.csv"

    with fits_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fit_rows[0]))
        writer.writeheader()
        writer.writerows(fit_rows)

    with preds_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(prediction_rows[0]))
        writer.writeheader()
        writer.writerows(prediction_rows)

    print("wrote", fits_path)
    print("wrote", preds_path)


if __name__ == "__main__":
    main()

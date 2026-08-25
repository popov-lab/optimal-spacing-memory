"""Ablate stochastic retrieval/encoding histories in the two-study MCM analyses.

The full MCM branches on retrieval success and encoding success. This script
compares it with mean-field replacements that remove either branch while
preserving the corresponding one-step expectation:

- retrieval mean field: epsilon_bar = 1 + p_recall * (epsilon_r - 1)
- encoding mean field: Delta x_i = omega * epsilon * (1 - s_i)

For Cepeda et al. (2009), the variants reuse the same forgetting-constrained
parameters from results/mcm_cepeda2009_fits.csv, because both mean-field
replacements preserve the single-study expectation exactly.

For Cepeda et al. (2008), where the forgetting data are unavailable, each
variant is fit directly to all 26 spacing observations under common broad
bounds. These are post-hoc flexibility checks, not predictive tests.

Usage:
    python src/ablate_mcm_histories.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit

from mcm import MCMParams, components

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "cepeda_spacing_recall.csv"
FITS_2009 = ROOT / "results" / "mcm_cepeda2009_fits.csv"
OUTPUT = ROOT / "results" / "mcm_history_ablation.csv"

EPSILON_R = 9.0
ZERO_GAP_2008 = 0.00256

VARIANTS = {
    "full stochastic": ("stochastic", "stochastic"),
    "no retrieval branching": ("mean", "stochastic"),
    "no encoding branching": ("stochastic", "mean"),
    "no branching": ("mean", "mean"),
}


def partial_strengths(x: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    return np.cumsum(gamma * x) / np.cumsum(gamma)


def two_session_variant(
    isi: float,
    ri: float,
    params: MCMParams,
    retrieval_mode: str,
    encoding_mode: str,
) -> float:
    """Exact two-study prediction under stochastic or mean-field branching."""

    tau, gamma = components(params)
    states: list[tuple[float, np.ndarray]] = [(1.0, np.zeros(params.n))]
    previous = 0.0

    for episode, study_time in enumerate((0.0, float(isi))):
        dt = 0.0 if episode == 0 else study_time - previous
        if dt:
            d = np.exp(-dt / tau)
            states = [(p, x * d) for p, x in states]

        next_states: list[tuple[float, np.ndarray]] = []
        for history_prob, x in states:
            s = partial_strengths(x, gamma)
            p_recall = float(np.clip(np.dot(gamma, x), 0.0, 1.0))

            if retrieval_mode == "stochastic":
                retrieval = [
                    (1.0 - p_recall, 1.0),
                    (p_recall, params.epsilon_r),
                ]
            elif retrieval_mode == "mean":
                epsilon = 1.0 + p_recall * (params.epsilon_r - 1.0)
                retrieval = [(1.0, epsilon)]
            else:
                raise ValueError(retrieval_mode)

            for retrieval_prob, epsilon in retrieval:
                if retrieval_prob <= 0:
                    continue
                branch_prob = history_prob * retrieval_prob

                if encoding_mode == "stochastic":
                    if 1.0 - params.omega > 0:
                        next_states.append(
                            (branch_prob * (1.0 - params.omega), x.copy())
                        )
                    next_states.append(
                        (
                            branch_prob * params.omega,
                            x + epsilon * (1.0 - s),
                        )
                    )
                elif encoding_mode == "mean":
                    next_states.append(
                        (
                            branch_prob,
                            x + params.omega * epsilon * (1.0 - s),
                        )
                    )
                else:
                    raise ValueError(encoding_mode)

        states = next_states
        previous = study_time

    d = np.exp(-float(ri) / tau)
    return sum(
        p * float(np.clip(np.dot(gamma, x * d), 0.0, 1.0))
        for p, x in states
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def load_2009_params() -> dict[str, MCMParams]:
    out: dict[str, MCMParams] = {}
    for row in read_csv(FITS_2009):
        out[row["panel"]] = MCMParams(
            mu=float(row["mu"]),
            nu=float(row["nu"]),
            omega=float(row["omega"]),
            xi=float(row["xi"]),
            epsilon_r=EPSILON_R,
            n=100,
        )
    return out


def predict_rows(
    rows: list[dict[str, str]],
    params: MCMParams,
    retrieval_mode: str,
    encoding_mode: str,
    zero_gap: float | None = None,
) -> np.ndarray:
    pred = []
    for row in rows:
        isi = float(row["isi_days"])
        if zero_gap is not None and isi == 0:
            isi = zero_gap
        pred.append(
            two_session_variant(
                isi,
                float(row["ri_days"]),
                params,
                retrieval_mode,
                encoding_mode,
            )
        )
    return np.asarray(pred)


# Broad bounds used only for the post-hoc 2008 direct-fit sensitivity check.
LOWER = np.array([
    np.log(1e-4),
    np.log(1e-3),
    logit(0.45),
    logit(0.05),
])
UPPER = np.array([
    np.log(1e3),
    np.log(10**1.5),
    logit(0.999),
    logit(0.999),
])


def unpack_2008(z: np.ndarray) -> MCMParams:
    return MCMParams(
        mu=float(np.exp(z[0])),
        nu=float(1.0 + np.exp(z[1])),
        omega=float(expit(z[2])),
        xi=float(expit(z[3])),
        epsilon_r=EPSILON_R,
        n=100,
    )


def fit_2008_variant(
    rows: list[dict[str, str]],
    observed: np.ndarray,
    retrieval_mode: str,
    encoding_mode: str,
    starts: int,
    seed: int,
) -> tuple[MCMParams, float]:
    rng = np.random.default_rng(seed)
    best_sse = np.inf
    best_params: MCMParams | None = None

    def residual(z: np.ndarray) -> np.ndarray:
        params = unpack_2008(z)
        return (
            predict_rows(
                rows,
                params,
                retrieval_mode,
                encoding_mode,
                zero_gap=ZERO_GAP_2008,
            )
            - observed
        )

    for _ in range(starts):
        initial = LOWER + rng.random(4) * (UPPER - LOWER)
        result = least_squares(
            residual,
            initial,
            bounds=(LOWER, UPPER),
            max_nfev=3500,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        sse = float(np.sum(result.fun**2))
        if sse < best_sse:
            best_sse = sse
            best_params = unpack_2008(result.x)

    if best_params is None:
        raise RuntimeError("all starts failed")
    return best_params, best_sse


def main() -> None:
    all_rows = read_csv(DATA)
    params_2009 = load_2009_params()
    output: list[dict[str, object]] = []

    pooled_errors: dict[str, list[float]] = {name: [] for name in VARIANTS}

    # 2009: genuine forgetting-constrained prediction.
    for panel in ("a", "b", "c"):
        rows = [
            r for r in all_rows
            if r["panel"] == panel and r["function"] == "spacing"
        ]
        rows.sort(key=lambda r: float(r["isi_days"]))
        observed = np.array([float(r["recall_pct"]) / 100.0 for r in rows])

        for name, (retrieval_mode, encoding_mode) in VARIANTS.items():
            predicted = predict_rows(
                rows,
                params_2009[panel],
                retrieval_mode,
                encoding_mode,
            )
            errors = predicted - observed
            pooled_errors[name].extend(errors.tolist())
            output.append(
                {
                    "dataset": "Cepeda 2009 prediction",
                    "panel": panel,
                    "variant": name,
                    "rmse_pp": 100.0 * np.sqrt(np.mean(errors**2)),
                }
            )

    for name, errors in pooled_errors.items():
        errors = np.asarray(errors)
        output.append(
            {
                "dataset": "Cepeda 2009 prediction",
                "panel": "pooled",
                "variant": name,
                "rmse_pp": 100.0 * np.sqrt(np.mean(errors**2)),
            }
        )

    # 2008: post-hoc direct fit to all 26 spacing observations.
    rows_2008 = [
        r for r in all_rows
        if r["panel"] == "d" and r["function"] == "spacing"
    ]
    rows_2008.sort(key=lambda r: (float(r["ri_days"]), float(r["isi_days"])))
    observed_2008 = np.array([float(r["recall_pct"]) / 100.0 for r in rows_2008])

    for offset, (name, (retrieval_mode, encoding_mode)) in enumerate(VARIANTS.items()):
        _, sse = fit_2008_variant(
            rows_2008,
            observed_2008,
            retrieval_mode,
            encoding_mode,
            starts=64,
            seed=20260825 + offset,
        )
        output.append(
            {
                "dataset": "Cepeda 2008 direct fit",
                "panel": "all 26",
                "variant": name,
                "rmse_pp": 100.0 * np.sqrt(sse / len(rows_2008)),
            }
        )

    with OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["dataset", "panel", "variant", "rmse_pp"],
        )
        writer.writeheader()
        for row in output:
            row = dict(row)
            row["rmse_pp"] = f"{float(row['rmse_pp']):.6f}"
            writer.writerow(row)

    print("wrote", OUTPUT)


if __name__ == "__main__":
    main()

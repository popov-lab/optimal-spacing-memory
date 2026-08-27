"""Fit two-event SAC with session-specific learning to Cepeda et al. (2008).

The original one-study-per-session abstraction is retained:

    Session 1 study = 0
    Session 2 study = ISI
    final test       = ISI + RI

The first study increment is fixed to one.  At Session 2, learning is a free
bounded gain applied to the remaining capacity at the current strength:

    u_1 = 1
    B(ISI) = f(ISI)
    u_2 = delta_2 * (1 - B(ISI))

with f(t) = (1 + t / tau)**(-d).  Thus the strengths underlying the observed
forgetting and spacing functions are

    B_F(a) = f(a)
    B_S(a, b) = f(a + b) + delta_2 * (1 - f(a)) * f(b).

A common logistic mapping converts strength to recall probability.  Two fits
are reported: a joint fit to all 11 forgetting and 26 spacing observations,
and a spacing-only fit to the 26 final-test observations.  A forgetting-only
fit cannot identify delta_2 because its observations precede Session 2 study.

Run from the repository root with:

    python3 src/fit_sac_cepeda2008_two_event.py
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit

mpl.rcParams["svg.hashsalt"] = "sac-cepeda2008-two-event"


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "cepeda_spacing_recall.csv"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

ZERO_ISI_DAYS = 0.00256
RI_LEVELS = (7.0, 35.0, 70.0, 350.0)
PROTOCOLS = ("joint", "spacing_only")
PROTOCOL_LABELS = {
    "joint": "Joint forgetting + spacing fit",
    "spacing_only": "Spacing-only fit",
}
COLORS = {
    7.0: "#0072B2",
    35.0: "#009E73",
    70.0: "#D55E00",
    350.0: "#CC79A7",
}


@dataclass(frozen=True)
class Cepeda2008Data:
    forgetting_rows: tuple[dict[str, str], ...]
    spacing_rows: tuple[dict[str, str], ...]
    observed_lag: np.ndarray
    model_lag: np.ndarray
    forgetting: np.ndarray
    observed_isi: np.ndarray
    model_isi: np.ndarray
    ri: np.ndarray
    spacing: np.ndarray


def _replace_nominal_zero(values: np.ndarray) -> np.ndarray:
    modeled = np.asarray(values, dtype=float).copy()
    modeled[modeled == 0.0] = ZERO_ISI_DAYS
    return modeled


def load_data(path: Path = DATA_PATH) -> Cepeda2008Data:
    """Load the 11 forgetting and 26 spacing observations in panel d."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = tuple(
            row for row in csv.DictReader(handle) if row["panel"] == "d"
        )
    forgetting_rows = tuple(
        sorted(
            (row for row in rows if row["function"] == "forgetting"),
            key=lambda row: float(row["isi_days"]),
        )
    )
    spacing_rows = tuple(
        sorted(
            (row for row in rows if row["function"] == "spacing"),
            key=lambda row: (float(row["ri_days"]), float(row["isi_days"])),
        )
    )
    observed_lag = np.array(
        [float(row["isi_days"]) for row in forgetting_rows]
    )
    observed_isi = np.array(
        [float(row["isi_days"]) for row in spacing_rows]
    )
    return Cepeda2008Data(
        forgetting_rows=forgetting_rows,
        spacing_rows=spacing_rows,
        observed_lag=observed_lag,
        model_lag=_replace_nominal_zero(observed_lag),
        forgetting=np.array(
            [float(row["recall_pct"]) / 100.0 for row in forgetting_rows]
        ),
        observed_isi=observed_isi,
        model_isi=_replace_nominal_zero(observed_isi),
        ri=np.array([float(row["ri_days"]) for row in spacing_rows]),
        spacing=np.array(
            [float(row["recall_pct"]) / 100.0 for row in spacing_rows]
        ),
    )


def _forgetting_from_log_tau(
    lag: np.ndarray | float, d: float, log_tau: float
) -> np.ndarray:
    """Evaluate shifted-power forgetting stably for very small tau."""
    lag_array = np.asarray(lag, dtype=float)
    if np.any(lag_array < 0.0):
        raise ValueError("Forgetting lags must be nonnegative")
    log_ratio = np.full_like(lag_array, -np.inf, dtype=float)
    positive = lag_array > 0.0
    log_ratio[positive] = np.log(lag_array[positive]) - log_tau
    return np.exp(-d * np.logaddexp(0.0, log_ratio))


def forgetting_function(
    lag: np.ndarray | float, d: float, tau: float
) -> np.ndarray:
    if d < 0.0:
        raise ValueError("d must be nonnegative")
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    return _forgetting_from_log_tau(lag, d, float(np.log(tau)))


def forgetting_strength_from_log_tau(
    lag: np.ndarray | float, d: float, log_tau: float
) -> np.ndarray:
    """Strength at the Session-2 pre-study test after u_1=1."""
    return _forgetting_from_log_tau(lag, d, log_tau)


def spacing_strength_from_log_tau(
    isi: np.ndarray | float,
    ri: np.ndarray | float,
    delta_2: float,
    d: float,
    log_tau: float,
) -> np.ndarray:
    """Final strength after u_1=1 and u_2=delta_2*(1-B(ISI))."""
    if not 0.0 < delta_2 <= 1.0:
        raise ValueError("delta_2 must be in (0, 1]")
    isi_array, ri_array = np.broadcast_arrays(
        np.asarray(isi, dtype=float), np.asarray(ri, dtype=float)
    )
    if np.any(isi_array < 0.0) or np.any(ri_array < 0.0):
        raise ValueError("ISI and RI must be nonnegative")
    f_isi = _forgetting_from_log_tau(isi_array, d, log_tau)
    f_ri = _forgetting_from_log_tau(ri_array, d, log_tau)
    f_total = _forgetting_from_log_tau(isi_array + ri_array, d, log_tau)
    return f_total + delta_2 * (1.0 - f_isi) * f_ri


def spacing_strength(
    isi: np.ndarray | float,
    ri: np.ndarray | float,
    delta_2: float,
    d: float,
    tau: float,
) -> np.ndarray:
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    return spacing_strength_from_log_tau(
        isi, ri, delta_2, d, float(np.log(tau))
    )


def response_probability(
    strength: np.ndarray | float, theta: float, sigma: float
) -> np.ndarray:
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    return expit((np.asarray(strength, dtype=float) - theta) / sigma)


def _parameter_bounds() -> tuple[np.ndarray, np.ndarray]:
    # The log(tau) coordinate is intentionally unbounded, matching the recent
    # two- and four-event Cepeda 2008 analyses.
    lower = np.array(
        [logit(1e-4), np.log(1e-4), -np.inf, -2.0, np.log(1e-3)]
    )
    upper = np.array(
        [logit(0.9999), np.log(10.0), np.inf, 3.0, np.log(5.0)]
    )
    return lower, upper


def _unpack(z: np.ndarray) -> dict[str, float]:
    logit_delta_2, log_d, log_tau, theta, log_sigma = z
    with np.errstate(over="ignore", under="ignore"):
        tau = float(np.exp(log_tau))
    return {
        "delta_1": 1.0,
        "delta_2": float(expit(logit_delta_2)),
        "d": float(np.exp(log_d)),
        "log_tau": float(log_tau),
        "tau": tau,
        "theta": float(theta),
        "sigma": float(np.exp(log_sigma)),
    }


def _initial_values(starts: int, seed: int) -> list[np.ndarray]:
    lower, upper = _parameter_bounds()
    initials: list[np.ndarray] = []
    for delta_2 in (0.05, 0.15, 0.30, 0.50, 0.75, 0.95):
        for d, tau in (
            (0.05, 1e-12),
            (0.08, 1e-6),
            (0.12, 0.01),
            (0.20, 0.10),
            (0.40, 1.0),
            (0.80, 30.0),
        ):
            initials.append(
                np.array(
                    [logit(delta_2), np.log(d), np.log(tau), 0.25, np.log(0.05)]
                )
            )
    rng = np.random.default_rng(seed)
    practical_lower = np.array(
        [logit(0.005), np.log(0.005), np.log(1e-16), -0.5, np.log(0.003)]
    )
    practical_upper = np.array(
        [logit(0.995), np.log(2.0), np.log(1e3), 1.5, np.log(0.8)]
    )
    while len(initials) < starts:
        z = practical_lower + rng.random(5) * (
            practical_upper - practical_lower
        )
        initials.append(np.minimum(np.maximum(z, lower), upper))
    return initials[:starts]


def predict(
    data: Cepeda2008Data, parameters: dict[str, float]
) -> tuple[np.ndarray, np.ndarray]:
    forgetting_prediction = response_probability(
        forgetting_strength_from_log_tau(
            data.model_lag, parameters["d"], parameters["log_tau"]
        ),
        parameters["theta"],
        parameters["sigma"],
    )
    spacing_prediction = response_probability(
        spacing_strength_from_log_tau(
            data.model_isi,
            data.ri,
            parameters["delta_2"],
            parameters["d"],
            parameters["log_tau"],
        ),
        parameters["theta"],
        parameters["sigma"],
    )
    return forgetting_prediction, spacing_prediction


def residual_vector(
    z: np.ndarray, data: Cepeda2008Data, protocol: str
) -> np.ndarray:
    forgetting_prediction, spacing_prediction = predict(data, _unpack(z))
    spacing_residual = spacing_prediction - data.spacing
    if protocol == "spacing_only":
        return spacing_residual
    if protocol == "joint":
        return np.concatenate(
            [forgetting_prediction - data.forgetting, spacing_residual]
        )
    raise ValueError(f"Unknown protocol: {protocol}")


def fit_model(
    data: Cepeda2008Data,
    protocol: str,
    starts: int = 128,
    seed: int = 20260827,
) -> dict[str, object]:
    """Fit one protocol by multistart least squares."""
    if protocol not in PROTOCOLS:
        raise ValueError(f"Unknown protocol: {protocol}")
    lower, upper = _parameter_bounds()
    solutions: list[tuple[float, object]] = []
    for z0 in _initial_values(starts, seed):
        fit = least_squares(
            residual_vector,
            z0,
            args=(data, protocol),
            bounds=(lower, upper),
            max_nfev=10000,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            x_scale="jac",
        )
        sse = float(fit.fun @ fit.fun)
        if np.isfinite(sse):
            solutions.append((sse, fit))
    if not solutions:
        raise RuntimeError(f"No finite fit for protocol={protocol}")
    solutions.sort(key=lambda item: item[0])
    objective_sse, best = solutions[0]
    parameters = _unpack(best.x)
    forgetting_prediction, spacing_prediction = predict(data, parameters)
    forgetting_residual = forgetting_prediction - data.forgetting
    spacing_residual = spacing_prediction - data.spacing
    singular_values = np.linalg.svd(best.jac, compute_uv=False)
    condition = (
        np.inf
        if singular_values[-1] == 0.0
        else float(singular_values[0] / singular_values[-1])
    )
    return {
        "protocol": protocol,
        "parameter_count": len(best.x),
        "objective_n": len(best.fun),
        "parameters": parameters,
        "forgetting_prediction": forgetting_prediction,
        "spacing_prediction": spacing_prediction,
        "forgetting_rmse_pp": 100.0
        * float(np.sqrt(np.mean(forgetting_residual**2))),
        "spacing_rmse_pp": 100.0
        * float(np.sqrt(np.mean(spacing_residual**2))),
        "joint_rmse_pp": 100.0
        * float(
            np.sqrt(
                np.mean(
                    np.concatenate([forgetting_residual, spacing_residual]) ** 2
                )
            )
        ),
        "objective_sse": objective_sse,
        "jacobian_condition": condition,
        "optimality": float(best.optimality),
        "nfev": int(best.nfev),
        "best_sse_values": [float(item[0]) for item in solutions[:10]],
    }


def continuous_optimum_isi(
    d: float, tau: float, ri: float, delta_2: float
) -> float:
    """Unconstrained nonnegative optimum for the two-event strength curve."""
    f_ri = float(forgetting_function(ri, d, tau))
    q = (delta_2 * f_ri) ** (-1.0 / (d + 1.0))
    return max(0.0, ri / (q - 1.0) - tau)


def _fit_row(fit: dict[str, object]) -> dict[str, object]:
    parameters = fit["parameters"]
    row: dict[str, object] = {
        "protocol": fit["protocol"],
        "parameters": fit["parameter_count"],
        "objective_n": fit["objective_n"],
        "delta_1": parameters["delta_1"],
        "delta_2": parameters["delta_2"],
        "d": parameters["d"],
        "tau_days": parameters["tau"],
        "log_tau_days": parameters["log_tau"],
        "theta": parameters["theta"],
        "sigma": parameters["sigma"],
        "forgetting_rmse_pp": fit["forgetting_rmse_pp"],
        "spacing_rmse_pp": fit["spacing_rmse_pp"],
        "joint_rmse_pp": fit["joint_rmse_pp"],
        "objective_sse": fit["objective_sse"],
        "jacobian_condition": fit["jacobian_condition"],
        "optimality": fit["optimality"],
        "nfev": fit["nfev"],
    }
    for ri in RI_LEVELS:
        row[f"optimum_isi_ri_{int(ri)}_days"] = continuous_optimum_isi(
            parameters["d"], parameters["tau"], ri, parameters["delta_2"]
        )
    return row


def write_fit_table(fits: list[dict[str, object]]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2008_two_event_fits.csv"
    rows = [_fit_row(fit) for fit in fits]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_prediction_table(
    fits: list[dict[str, object]], data: Cepeda2008Data
) -> Path:
    path = RESULTS_DIR / "sac_cepeda2008_two_event_predictions.csv"
    fieldnames = [
        "protocol",
        "function",
        "isi_days",
        "model_isi_days",
        "ri_days",
        "in_objective",
        "observed_pct",
        "predicted_pct",
        "residual_pp",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        for fit in fits:
            for function, observed_x, model_x, ri, observed, predicted in (
                (
                    "forgetting",
                    data.observed_lag,
                    data.model_lag,
                    np.full(data.observed_lag.shape, np.nan),
                    data.forgetting,
                    fit["forgetting_prediction"],
                ),
                (
                    "spacing",
                    data.observed_isi,
                    data.model_isi,
                    data.ri,
                    data.spacing,
                    fit["spacing_prediction"],
                ),
            ):
                in_objective = function == "spacing" or fit["protocol"] == "joint"
                for x_i, model_x_i, ri_i, observed_i, predicted_i in zip(
                    observed_x, model_x, ri, observed, predicted
                ):
                    writer.writerow(
                        {
                            "protocol": fit["protocol"],
                            "function": function,
                            "isi_days": x_i,
                            "model_isi_days": model_x_i,
                            "ri_days": "" if np.isnan(ri_i) else ri_i,
                            "in_objective": in_objective,
                            "observed_pct": 100.0 * observed_i,
                            "predicted_pct": 100.0 * predicted_i,
                            "residual_pp": 100.0 * (predicted_i - observed_i),
                        }
                    )
    return path


def _style_axis(axis: plt.Axes) -> None:
    axis.set_ylim(-2.0, 102.0)
    axis.set_yticks([0, 20, 40, 60, 80, 100])
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.8, alpha=0.75)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=8.7)


def _normalize_generated_svg(path: Path) -> None:
    content = path.read_text(encoding="utf-8")
    path.write_text(
        "\n".join(line.rstrip() for line in content.splitlines()) + "\n",
        encoding="utf-8",
    )


def plot_fit(fit: dict[str, object], data: Cepeda2008Data) -> Path:
    parameters = fit["parameters"]
    figure = plt.figure(figsize=(14.2, 5.5), layout="constrained")
    grid = figure.add_gridspec(1, 3, width_ratios=(1.15, 2.1, 0.82))
    forgetting_axis = figure.add_subplot(grid[0, 0])
    spacing_axis = figure.add_subplot(grid[0, 1])
    info_axis = figure.add_subplot(grid[0, 2])
    info_axis.set_axis_off()

    figure.suptitle(
        "Cepeda et al. (2008): two-event SAC with "
        "$\\delta_1=1$ and free $\\delta_2$\n"
        + PROTOCOL_LABELS[fit["protocol"]],
        fontsize=14.5,
        fontweight="semibold",
    )

    lag_grid = np.linspace(ZERO_ISI_DAYS, 105.0, 800)
    forgetting_curve = response_probability(
        forgetting_strength_from_log_tau(
            lag_grid, parameters["d"], parameters["log_tau"]
        ),
        parameters["theta"],
        parameters["sigma"],
    )
    forgetting_axis.plot(
        lag_grid, 100.0 * forgetting_curve, color="#3A3A3A", linewidth=2.3
    )
    forgetting_axis.plot(
        data.observed_lag,
        100.0 * data.forgetting,
        linestyle="none",
        marker="s",
        markersize=6.2,
        markerfacecolor="white",
        markeredgecolor="#202020",
        markeredgewidth=1.4,
        zorder=3,
    )
    forgetting_axis.set_xlim(-4.0, 109.0)
    forgetting_axis.set_xticks([0, 21, 35, 70, 105])
    forgetting_axis.set_xlabel("Study-test lag (days)")
    forgetting_axis.set_ylabel("Recall (%)")
    forgetting_axis.set_title("Session-2 pre-study test")
    _style_axis(forgetting_axis)
    forgetting_axis.text(
        0.97,
        0.95,
        f"RMSE = {fit['forgetting_rmse_pp']:.2f} pp",
        transform=forgetting_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.7,
    )

    isi_grid = np.linspace(0.0, 105.0, 900)
    model_isi_grid = _replace_nominal_zero(isi_grid)
    legend_handles = []
    legend_labels = []
    for ri in RI_LEVELS:
        color = COLORS[ri]
        ri_grid = np.full_like(model_isi_grid, ri)
        spacing_curve = response_probability(
            spacing_strength_from_log_tau(
                model_isi_grid,
                ri_grid,
                parameters["delta_2"],
                parameters["d"],
                parameters["log_tau"],
            ),
            parameters["theta"],
            parameters["sigma"],
        )
        line, = spacing_axis.plot(
            isi_grid, 100.0 * spacing_curve, color=color, linewidth=2.2
        )
        condition = data.ri == ri
        points = spacing_axis.plot(
            data.observed_isi[condition],
            100.0 * data.spacing[condition],
            linestyle="none",
            marker="o",
            markersize=6.0,
            markerfacecolor=color,
            markeredgecolor="#202020",
            markeredgewidth=0.6,
            zorder=3,
        )[0]
        optimum = continuous_optimum_isi(
            parameters["d"], parameters["tau"], ri, parameters["delta_2"]
        )
        if 0.0 < optimum < 105.0:
            spacing_axis.axvline(
                optimum,
                color=color,
                linewidth=1.0,
                alpha=0.55,
                linestyle=(0, (1.5, 2.5)),
            )
        legend_handles.extend([line, points])
        legend_labels.extend([f"Fit, RI={ri:g} d", f"Data, RI={ri:g} d"])

    spacing_axis.set_xlim(-5.25, 110.25)
    spacing_axis.set_xticks([0, 7, 14, 21, 35, 70, 105])
    spacing_axis.set_xlabel("ISI (days)")
    spacing_axis.set_ylabel("Final-test recall (%)")
    spacing_axis.set_title("Final spacing test")
    _style_axis(spacing_axis)
    spacing_axis.text(
        0.98,
        0.95,
        f"RMSE = {fit['spacing_rmse_pp']:.2f} pp",
        transform=spacing_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.7,
    )

    parameter_lines = [
        "Parameters",
        "delta_1 = 1 (fixed)",
        f"delta_2 = {parameters['delta_2']:.4f}",
        f"d = {parameters['d']:.4f}",
        f"tau = {parameters['tau']:.5g} d",
        f"theta = {parameters['theta']:.4f}",
        f"sigma = {parameters['sigma']:.4f}",
        "",
        "RMSE",
        f"forgetting = {fit['forgetting_rmse_pp']:.2f} pp",
        f"spacing = {fit['spacing_rmse_pp']:.2f} pp",
        f"all 37 = {fit['joint_rmse_pp']:.2f} pp",
        "",
        "Dotted lines: continuous optima",
    ]
    info_axis.text(
        0.0,
        0.99,
        "\n".join(parameter_lines),
        transform=info_axis.transAxes,
        ha="left",
        va="top",
        fontsize=9.2,
        linespacing=1.32,
    )
    info_axis.legend(
        legend_handles,
        legend_labels,
        loc="lower left",
        frameon=False,
        fontsize=8.0,
        borderaxespad=0.0,
    )

    stem = FIGURES_DIR / f"sac_cepeda2008_two_event_{fit['protocol']}"
    svg_path = stem.with_suffix(".svg")
    figure.savefig(svg_path, facecolor="white", metadata={"Date": None})
    _normalize_generated_svg(svg_path)
    plt.close(figure)
    return svg_path


def _format_tau(value: float) -> str:
    if value < 1e-4 or value >= 1e4:
        return f"{value:.3e}"
    return f"{value:.5f}"


def write_report(fits: list[dict[str, object]]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2008_two_event_report.md"
    by_protocol = {fit["protocol"]: fit for fit in fits}
    lines = [
        "# Two-event SAC fits to Cepeda et al. (2008)",
        "",
        "## Model",
        "",
        "This analysis returns to the original abstraction of one study event per session. The first-session increment is fixed at one, while the second-session learning rate is estimated and acts on the unlearned proportion at the current strength:",
        "",
        "$$",
        "u_1=1,\\qquad B(a)=f(a),\\qquad u_2=\\delta_2[1-B(a)].",
        "$$",
        "",
        "Consequently, the strengths underlying the forgetting and spacing observations are",
        "",
        "$$",
        "B_F(a)=f(a),",
        "$$",
        "",
        "$$",
        "B_S(a,b)=f(a+b)+\\delta_2[1-f(a)]f(b),",
        "$$",
        "",
        "where $a$ is the ISI, $b$ is the RI, and $f(t)=(1+t/\\tau)^{-d}$. A shared logistic response function maps strength to recall probability:",
        "",
        "$$",
        "p(B)=\\operatorname{logistic}\\!\\left(\\frac{B-\\theta}{\\sigma}\\right).",
        "$$",
        "",
        "The update itself is $u_2=\\delta_2[1-B(a)]$; its contribution at the final test is then multiplied by $f(b)$. The pre-update interval $a$ does not decay a trace that is only created at Session 2.",
        "",
        "The nominal zero-day condition is represented as the reported approximately 3-minute interval, $0.00256$ days. All residuals are unweighted point-level recall-probability residuals; no trial denominators or standard errors are available.",
        "",
        "## Fit targets",
        "",
        "Two requested fits are reported:",
        "",
        "- **Joint:** estimates $\\{\\delta_2,d,\\tau,\\theta,\\sigma\\}$ from all 11 forgetting and 26 spacing observations.",
        "- **Spacing only:** estimates the same five parameters from the 26 spacing observations; the forgetting curve is an out-of-sample diagnostic.",
        "",
        "A forgetting-only fit cannot estimate $\\delta_2$, because every forgetting observation is the Session-2 pre-study test and therefore occurs before the update governed by $\\delta_2$. A purported forgetting-only prediction with a freely estimated second-session learning rate would be unidentified; it would require fixing or calibrating $\\delta_2$ from some spacing data.",
        "",
        "## Results",
        "",
        "| Fit target | Parameters | $\\delta_2$ | $d$ | $\\tau$ (days) | $\\theta$ | $\\sigma$ | Forgetting RMSE | Spacing RMSE |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for protocol in PROTOCOLS:
        fit = by_protocol[protocol]
        parameters = fit["parameters"]
        lines.append(
            f"| {PROTOCOL_LABELS[protocol]} | {fit['parameter_count']} | "
            f"{parameters['delta_2']:.5f} | {parameters['d']:.5f} | "
            f"{_format_tau(parameters['tau'])} | {parameters['theta']:.5f} | "
            f"{parameters['sigma']:.5f} | {fit['forgetting_rmse_pp']:.2f} pp | "
            f"{fit['spacing_rmse_pp']:.2f} pp |"
        )
    joint = by_protocol["joint"]
    spacing_only = by_protocol["spacing_only"]
    lines.extend(
        [
            "",
            "The joint fit is the direct test of whether one parameterization can reconcile the single-study forgetting function with the four spacing curves. The spacing-only fit shows the best descriptive fit available to this session-specific-learning version and reveals what it sacrifices on the forgetting curve.",
            "",
            f"For comparison, the preceding two-event spacing-only model used one common $\\delta$ at both studies and obtained 4.02 pp spacing RMSE. Fixing the first increment to one and estimating only $\\delta_2$ changes that value to {spacing_only['spacing_rmse_pp']:.2f} pp.",
            "",
            "## Interpretation",
            "",
            f"The session-specific learning rates improve the descriptive spacing-only fit by {4.0199395450585085 - spacing_only['spacing_rmse_pp']:.2f} percentage points relative to the preceding common-$\\delta$ two-event fit. They do not reconcile the two datasets: imposing the forgetting observations raises spacing RMSE from {spacing_only['spacing_rmse_pp']:.2f} to {joint['spacing_rmse_pp']:.2f} pp, while the spacing-only solution misses the forgetting curve by {spacing_only['forgetting_rmse_pp']:.2f} pp.",
            "",
            f"The estimated time scale also changes sharply with the target. The joint fit places $\\tau$ at {_format_tau(joint['parameters']['tau'])} days (about {24.0 * 60.0 * joint['parameters']['tau']:.1f} minutes), below the shortest modeled lag of about 3.7 minutes. The spacing-only fit places it at {_format_tau(spacing_only['parameters']['tau'])} days (about {24.0 * spacing_only['parameters']['tau']:.1f} hours). Thus fixing the first increment exposes, rather than removes, the tension between the forgetting curve and the spacing surface.",
            "",
            "## Continuous optimal ISIs",
            "",
            "For this model, an interior optimum satisfies",
            "",
            "$$",
            "a^*=\\frac{b}{[\\delta_2 f(b)]^{-1/(d+1)}-1}-\\tau,",
            "$$",
            "",
            "with negative values replaced by the boundary $a^*=0$.",
            "",
            "| Fit target | RI = 7 d | RI = 35 d | RI = 70 d | RI = 350 d |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for protocol in PROTOCOLS:
        row = _fit_row(by_protocol[protocol])
        lines.append(
            f"| {PROTOCOL_LABELS[protocol]} | "
            f"{row['optimum_isi_ri_7_days']:.2f} d | "
            f"{row['optimum_isi_ri_35_days']:.2f} d | "
            f"{row['optimum_isi_ri_70_days']:.2f} d | "
            f"{row['optimum_isi_ri_350_days']:.2f} d |"
        )
    lines.extend(
        [
            "",
            "## Diagnostic figures",
            "",
            "### Joint fit",
            "",
            "![Joint two-event SAC fit](../figures/sac_cepeda2008_two_event_joint.svg)",
            "",
            "### Spacing-only fit",
            "",
            "![Spacing-only two-event SAC fit](../figures/sac_cepeda2008_two_event_spacing_only.svg)",
            "",
            "Point predictions are in [`sac_cepeda2008_two_event_predictions.csv`](sac_cepeda2008_two_event_predictions.csv), and full-precision parameters and diagnostics are in [`sac_cepeda2008_two_event_fits.csv`](sac_cepeda2008_two_event_fits.csv).",
            "",
            "## Reproduction",
            "",
            "```bash",
            "python3 src/fit_sac_cepeda2008_two_event.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260827)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    data = load_data()
    fits: list[dict[str, object]] = []
    for index, protocol in enumerate(PROTOCOLS):
        print(f"Fitting {protocol} ({args.starts} starts)", flush=True)
        fit = fit_model(
            data, protocol, starts=args.starts, seed=args.seed + 1000 * index
        )
        fits.append(fit)
        parameters = fit["parameters"]
        print(
            f"  delta_2={parameters['delta_2']:.8g}, "
            f"d={parameters['d']:.8g}, tau={parameters['tau']:.8g}, "
            f"forgetting RMSE={fit['forgetting_rmse_pp']:.3f} pp, "
            f"spacing RMSE={fit['spacing_rmse_pp']:.3f} pp",
            flush=True,
        )

    fit_path = write_fit_table(fits)
    prediction_path = write_prediction_table(fits, data)
    figure_paths = [plot_fit(fit, data) for fit in fits]
    report_path = write_report(fits)
    print("Wrote", fit_path)
    print("Wrote", prediction_path)
    for figure_path in figure_paths:
        print("Wrote", figure_path)
    print("Wrote", report_path)


if __name__ == "__main__":
    main()

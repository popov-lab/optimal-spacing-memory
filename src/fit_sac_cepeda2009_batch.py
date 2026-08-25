"""Fit session-batch SAC to the Cepeda et al. (2009) data.

The model retains the nonlinear consequences of repeated learning events while
collapsing all events within a session to the same time.  If a session contains
``m`` events with a common SAC learning rate ``delta``, its immediate strength
from an initially absent trace is

    A_m = 1 - (1 - delta)**m.

Session 1 contains ``m1`` events and Session 2 contains two events.  The first
Session-2 test measures the pre-feedback strength.  The final spacing-test
strength is the surviving Session-1 contribution plus the Session-2 increment.

Three protocols are fitted for each assumed Experiment-1 event count m1=3..7:

* forgetting_calibrated: fit the 18 forgetting observations conditional on
  delta, and select delta using only the Experiment-1 spacing curve;
* joint: fit all 36 forgetting and spacing observations;
* spacing_only: fit the 18 spacing observations without forgetting data.

Experiment 2 has four Session-1 events.  All experiments have two Session-2
events.  The nominal zero-day gaps are corrected to 5 minutes in Experiment 1
and 20 minutes in Experiment 2.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import expit, logit


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "cepeda_spacing_recall.csv"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

PANEL_KEYS = ("a", "b", "c")
PANEL_LABELS = {
    "a": "Experiment 1\nSwahili-English",
    "b": "Experiment 2a\nObscure facts",
    "c": "Experiment 2b\nObject names",
}
PROTOCOL_LABELS = {
    "forgetting_calibrated": "Forgetting fit + Exp. 1 delta calibration",
    "joint": "Joint fit",
    "spacing_only": "Spacing-only fit",
}
SHORTEST_ISI_DAYS = {
    "a": 5.0 / (24.0 * 60.0),
    "b": 20.0 / (24.0 * 60.0),
    "c": 20.0 / (24.0 * 60.0),
}
M1_EXPERIMENT_2 = 4
M2_ALL = 2


@dataclass(frozen=True)
class Panel:
    key: str
    lag: np.ndarray
    forgetting: np.ndarray
    isi: np.ndarray
    ri: float
    spacing: np.ndarray


def load_panels(path: Path = DATA_PATH) -> dict[str, Panel]:
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["panel"] in PANEL_KEYS:
                rows.append(row)

    panels: dict[str, Panel] = {}
    for key in PANEL_KEYS:
        forgetting_rows = [
            row for row in rows
            if row["panel"] == key and row["function"] == "forgetting"
        ]
        spacing_rows = [
            row for row in rows
            if row["panel"] == key and row["function"] == "spacing"
        ]
        correction = SHORTEST_ISI_DAYS[key]
        lag = np.array([float(row["isi_days"]) for row in forgetting_rows])
        isi = np.array([float(row["isi_days"]) for row in spacing_rows])
        lag[lag == 0.0] = correction
        isi[isi == 0.0] = correction
        panels[key] = Panel(
            key=key,
            lag=lag,
            forgetting=np.array(
                [float(row["recall_pct"]) / 100.0 for row in forgetting_rows]
            ),
            isi=isi,
            ri=float(spacing_rows[0]["ri_days"]),
            spacing=np.array(
                [float(row["recall_pct"]) / 100.0 for row in spacing_rows]
            ),
        )
    return panels


def forgetting_function(t: np.ndarray | float, d: float, tau: float) -> np.ndarray:
    return (1.0 + np.asarray(t, dtype=float) / tau) ** (-d)


def session_gain(m: int, delta: float) -> float:
    return 1.0 - (1.0 - delta) ** m


def response_probability(strength: np.ndarray, theta: float, sigma: float) -> np.ndarray:
    return expit((np.asarray(strength) - theta) / sigma)


def unpack(z: np.ndarray) -> dict[str, object]:
    ds = np.exp(z[:3])
    taus = np.exp(z[3:6])
    theta = float(z[6])
    sigma = float(np.exp(z[7]))
    delta = float(expit(z[8]))
    return {
        "panel": {
            key: {"d": float(ds[j]), "tau": float(taus[j])}
            for j, key in enumerate(PANEL_KEYS)
        },
        "theta": theta,
        "sigma": sigma,
        "delta": delta,
    }


def panel_m1(key: str, experiment1_m1: int) -> int:
    return experiment1_m1 if key == "a" else M1_EXPERIMENT_2


def forgetting_strength(
    t: np.ndarray, d: float, tau: float, delta: float, m1: int
) -> np.ndarray:
    return session_gain(m1, delta) * forgetting_function(t, d, tau)


def spacing_strength(
    a: np.ndarray, b: float, d: float, tau: float, delta: float, m1: int
) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    gain1 = session_gain(m1, delta)
    gain2 = session_gain(M2_ALL, delta)
    f_a = forgetting_function(a, d, tau)
    f_b = forgetting_function(b, d, tau)
    f_ab = forgetting_function(a + b, d, tau)
    return gain1 * f_ab + gain2 * (1.0 - gain1 * f_a) * f_b


def predict_forgetting(panel: Panel, p: dict[str, float], pars: dict[str, object], m1: int) -> np.ndarray:
    strength = forgetting_strength(
        panel.lag, p["d"], p["tau"], pars["delta"], m1
    )
    return response_probability(strength, pars["theta"], pars["sigma"])


def predict_spacing(panel: Panel, p: dict[str, float], pars: dict[str, object], m1: int) -> np.ndarray:
    strength = spacing_strength(
        panel.isi, panel.ri, p["d"], p["tau"], pars["delta"], m1
    )
    return response_probability(strength, pars["theta"], pars["sigma"])


def residual_vector(
    z: np.ndarray,
    panels: dict[str, Panel],
    experiment1_m1: int,
    protocol: str,
) -> np.ndarray:
    pars = unpack(z)
    out: list[float] = []
    if protocol in {"forgetting_only", "joint"}:
        for key in PANEL_KEYS:
            panel = panels[key]
            prediction = predict_forgetting(
                panel, pars["panel"][key], pars,
                panel_m1(key, experiment1_m1),
            )
            out.extend(prediction - panel.forgetting)
    if protocol in {"spacing_only", "joint"}:
        for key in PANEL_KEYS:
            panel = panels[key]
            prediction = predict_spacing(
                panel, pars["panel"][key], pars,
                panel_m1(key, experiment1_m1),
            )
            out.extend(prediction - panel.spacing)
    return np.asarray(out)


def parameter_bounds() -> tuple[np.ndarray, np.ndarray]:
    lower = np.array(
        [np.log(1e-3)] * 3
        + [np.log(1e-5)] * 3
        + [-3.0, np.log(1e-3), logit(1e-4)]
    )
    upper = np.array(
        [np.log(3.0)] * 3
        + [np.log(1e4)] * 3
        + [3.0, np.log(2.0), logit(0.9999)]
    )
    return lower, upper


def initial_values(starts: int, seed: int) -> list[np.ndarray]:
    lower, upper = parameter_bounds()
    initials: list[np.ndarray] = []

    previous = (
        [0.1932, 0.1736, 0.2298],
        [0.7237, 3.1758, 1.9105],
        0.6719,
        0.1209,
        0.626,
    )
    d_sets = (
        previous[0],
        [0.08, 0.08, 0.08],
        [0.20, 0.20, 0.20],
        [0.45, 0.30, 0.40],
        [0.80, 0.80, 0.80],
    )
    tau_sets = (
        previous[1],
        [0.01, 0.01, 0.01],
        [0.20, 1.0, 0.5],
        [1.0, 10.0, 3.0],
        [10.0, 100.0, 30.0],
    )
    for delta in (0.03, 0.10, 0.30, 0.626, 0.85, 0.98):
        for ds, taus in zip(d_sets, tau_sets):
            initials.append(np.array(
                [*(np.log(ds)), *(np.log(taus)), previous[2],
                 np.log(previous[3]), logit(delta)]
            ))

    rng = np.random.default_rng(seed)
    practical_lower = np.array(
        [np.log(0.01)] * 3
        + [np.log(5e-4)] * 3
        + [-0.5, np.log(0.008), logit(0.005)]
    )
    practical_upper = np.array(
        [np.log(1.5)] * 3
        + [np.log(500.0)] * 3
        + [1.5, np.log(0.8), logit(0.995)]
    )
    while len(initials) < starts:
        z = practical_lower + rng.random(9) * (practical_upper - practical_lower)
        initials.append(np.minimum(np.maximum(z, lower), upper))
    return initials[:starts]


def fit_model(
    panels: dict[str, Panel],
    experiment1_m1: int,
    protocol: str,
    starts: int,
    seed: int,
) -> dict[str, object]:
    lower, upper = parameter_bounds()
    solutions: list[tuple[float, np.ndarray, np.ndarray, float, int]] = []
    for z0 in initial_values(starts, seed):
        fit = least_squares(
            residual_vector,
            z0,
            args=(panels, experiment1_m1, protocol),
            bounds=(lower, upper),
            max_nfev=7000,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        sse = float(fit.fun @ fit.fun)
        if np.isfinite(sse):
            solutions.append((sse, fit.x, fit.jac, fit.optimality, fit.nfev))
    if not solutions:
        raise RuntimeError(f"No finite fit for {protocol}, m1={experiment1_m1}")
    solutions.sort(key=lambda item: item[0])
    sse, z, jacobian, optimality, nfev = solutions[0]
    singular = np.linalg.svd(jacobian, compute_uv=False)
    condition = np.inf if singular[-1] == 0 else float(singular[0] / singular[-1])
    return {
        "protocol": protocol,
        "experiment1_m1": experiment1_m1,
        "parameters": unpack(z),
        "objective_sse": sse,
        "jacobian_condition": condition,
        "jacobian_singular_values": singular,
        "optimality": float(optimality),
        "nfev": int(nfev),
        "next_sse": [float(item[0]) for item in solutions[:10]],
    }


def fit_forgetting_given_delta(
    panels: dict[str, Panel],
    experiment1_m1: int,
    delta: float,
    starts: int,
    seed: int,
    warm: np.ndarray | None = None,
) -> dict[str, object]:
    """Optimize all forgetting parameters while holding delta fixed."""
    lower, upper = parameter_bounds()
    lower8, upper8 = lower[:8], upper[:8]
    candidates: list[np.ndarray] = []
    if warm is not None:
        candidates.append(np.asarray(warm, dtype=float))

    # Structured full-model initials repeat the first eight coordinates across
    # delta values. Deduplicate them and retain enough random starts.
    for z9 in initial_values(30 + starts, seed):
        z8 = z9[:8]
        if not any(np.allclose(z8, old, rtol=0.0, atol=1e-12) for old in candidates):
            candidates.append(z8)
        if len(candidates) >= starts:
            break

    fixed_delta_z = logit(delta)

    def conditional_residual(z8: np.ndarray) -> np.ndarray:
        z9 = np.concatenate([z8, [fixed_delta_z]])
        return residual_vector(z9, panels, experiment1_m1, "forgetting_only")

    solutions: list[tuple[float, np.ndarray, np.ndarray, float, int]] = []
    for z0 in candidates:
        fit = least_squares(
            conditional_residual,
            z0,
            bounds=(lower8, upper8),
            max_nfev=6000,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        sse = float(fit.fun @ fit.fun)
        if np.isfinite(sse):
            solutions.append((sse, fit.x, fit.jac, fit.optimality, fit.nfev))
    if not solutions:
        raise RuntimeError(
            f"No conditional forgetting fit for delta={delta}, m1={experiment1_m1}"
        )
    solutions.sort(key=lambda item: item[0])
    sse, z8, jacobian, optimality, nfev = solutions[0]
    singular = np.linalg.svd(jacobian, compute_uv=False)
    condition = np.inf if singular[-1] == 0 else float(singular[0] / singular[-1])
    z9 = np.concatenate([z8, [fixed_delta_z]])
    return {
        "z8": z8,
        "parameters": unpack(z9),
        "forgetting_sse": sse,
        "jacobian_condition": condition,
        "jacobian_singular_values": singular,
        "optimality": float(optimality),
        "nfev": int(nfev),
    }


def experiment1_spacing_sse(
    conditional_fit: dict[str, object],
    panels: dict[str, Panel],
    experiment1_m1: int,
) -> float:
    pars = conditional_fit["parameters"]
    panel = panels["a"]
    prediction = predict_spacing(
        panel, pars["panel"]["a"], pars, experiment1_m1
    )
    residual = prediction - panel.spacing
    return float(residual @ residual)


def fit_forgetting_calibrated(
    panels: dict[str, Panel],
    experiment1_m1: int,
    starts: int,
    seed: int,
    grid_size: int,
) -> dict[str, object]:
    """Bilevel fit: forgetting determines all parameters except delta.

    For every candidate delta, the remaining parameters minimize forgetting
    SSE. Delta is then selected solely by Experiment 1's spacing SSE.
    """
    conditional_starts = max(10, starts // 8)
    delta_grid = np.linspace(0.02, 0.98, grid_size)
    grid_fits: list[tuple[float, dict[str, object]]] = []
    warm: np.ndarray | None = None
    for index, delta in enumerate(delta_grid):
        candidate = fit_forgetting_given_delta(
            panels,
            experiment1_m1,
            float(delta),
            conditional_starts,
            seed + index,
            warm=warm,
        )
        warm = candidate["z8"]
        grid_fits.append((float(delta), candidate))

    grid_objectives = np.array([
        experiment1_spacing_sse(candidate, panels, experiment1_m1)
        for _, candidate in grid_fits
    ])
    best_index = int(np.argmin(grid_objectives))
    left = float(delta_grid[max(0, best_index - 1)])
    right = float(delta_grid[min(len(delta_grid) - 1, best_index + 1)])

    cache: dict[float, dict[str, object]] = {
        round(delta, 12): candidate for delta, candidate in grid_fits
    }
    best_warm = grid_fits[best_index][1]["z8"]

    def outer_objective(delta: float) -> float:
        key = round(float(delta), 12)
        if key not in cache:
            cache[key] = fit_forgetting_given_delta(
                panels,
                experiment1_m1,
                float(delta),
                conditional_starts,
                seed + 10000 + len(cache),
                warm=best_warm,
            )
        return experiment1_spacing_sse(cache[key], panels, experiment1_m1)

    if right > left:
        refinement = minimize_scalar(
            outer_objective,
            bounds=(left, right),
            method="bounded",
            options={"xatol": 1e-7, "maxiter": 80},
        )
        calibrated_delta = float(refinement.x)
        outer_objective(calibrated_delta)
        conditional = cache[round(calibrated_delta, 12)]
    else:
        calibrated_delta, conditional = grid_fits[best_index]

    calibration_sse = experiment1_spacing_sse(
        conditional, panels, experiment1_m1
    )
    return {
        "protocol": "forgetting_calibrated",
        "experiment1_m1": experiment1_m1,
        "parameters": conditional["parameters"],
        "objective_sse": conditional["forgetting_sse"],
        "calibration_sse": calibration_sse,
        "jacobian_condition": conditional["jacobian_condition"],
        "jacobian_singular_values": conditional["jacobian_singular_values"],
        "optimality": conditional["optimality"],
        "nfev": conditional["nfev"],
        "delta_grid": delta_grid,
        "delta_grid_calibration_sse": grid_objectives,
    }


def continuous_optimum(d: float, tau: float, ri: float, delta: float) -> float:
    gain2 = session_gain(M2_ALL, delta)
    f_b = float(forgetting_function(ri, d, tau))
    q = (gain2 * f_b) ** (-1.0 / (d + 1.0))
    return max(0.0, ri / (q - 1.0) - tau)


def add_metrics(fit: dict[str, object], panels: dict[str, Panel]) -> None:
    pars = fit["parameters"]
    experiment1_m1 = fit["experiment1_m1"]
    forgetting_residuals: list[float] = []
    spacing_residuals: list[float] = []
    metrics: dict[str, object] = {"panel": {}}
    for key in PANEL_KEYS:
        panel = panels[key]
        p = pars["panel"][key]
        m1 = panel_m1(key, experiment1_m1)
        forgetting_prediction = predict_forgetting(panel, p, pars, m1)
        spacing_prediction = predict_spacing(panel, p, pars, m1)
        forgetting_residual = forgetting_prediction - panel.forgetting
        spacing_residual = spacing_prediction - panel.spacing
        forgetting_residuals.extend(forgetting_residual)
        spacing_residuals.extend(spacing_residual)
        metrics["panel"][key] = {
            "forgetting_rmse_pp": 100.0 * float(np.sqrt(np.mean(forgetting_residual**2))),
            "spacing_rmse_pp": 100.0 * float(np.sqrt(np.mean(spacing_residual**2))),
            "optimum_isi_days": continuous_optimum(
                p["d"], p["tau"], panel.ri, pars["delta"]
            ),
            "forgetting_prediction": forgetting_prediction,
            "spacing_prediction": spacing_prediction,
        }
    forgetting_residuals_array = np.asarray(forgetting_residuals)
    spacing_residuals_array = np.asarray(spacing_residuals)
    metrics["forgetting_rmse_pp"] = 100.0 * float(
        np.sqrt(np.mean(forgetting_residuals_array**2))
    )
    metrics["spacing_rmse_pp"] = 100.0 * float(
        np.sqrt(np.mean(spacing_residuals_array**2))
    )
    metrics["joint_rmse_pp"] = 100.0 * float(np.sqrt(np.mean(np.concatenate(
        [forgetting_residuals_array, spacing_residuals_array]
    ) ** 2)))
    fit["metrics"] = metrics


def fit_row(fit: dict[str, object]) -> dict[str, object]:
    pars = fit["parameters"]
    metrics = fit["metrics"]
    row: dict[str, object] = {
        "protocol": fit["protocol"],
        "experiment1_m1": fit["experiment1_m1"],
        "experiment2_m1": M1_EXPERIMENT_2,
        "session2_m2": M2_ALL,
        "delta": pars["delta"],
        "theta": pars["theta"],
        "sigma": pars["sigma"],
        "forgetting_rmse_pp": metrics["forgetting_rmse_pp"],
        "spacing_rmse_pp": metrics["spacing_rmse_pp"],
        "joint_rmse_pp": metrics["joint_rmse_pp"],
        "objective_sse": fit["objective_sse"],
        "calibration_sse": fit.get("calibration_sse", ""),
        "jacobian_condition": fit["jacobian_condition"],
    }
    for key in PANEL_KEYS:
        row[f"d_{key}"] = pars["panel"][key]["d"]
        row[f"tau_{key}_days"] = pars["panel"][key]["tau"]
        row[f"forgetting_rmse_{key}_pp"] = metrics["panel"][key]["forgetting_rmse_pp"]
        row[f"spacing_rmse_{key}_pp"] = metrics["panel"][key]["spacing_rmse_pp"]
        row[f"optimum_{key}_days"] = metrics["panel"][key]["optimum_isi_days"]
    return row


def write_fit_table(fits: list[dict[str, object]]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2009_batch_fits.csv"
    rows = [fit_row(fit) for fit in fits]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_prediction_table(fits: list[dict[str, object]], panels: dict[str, Panel]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2009_batch_predictions.csv"
    fieldnames = [
        "protocol", "experiment1_m1", "panel", "function", "isi_days",
        "ri_days", "observed_pct", "predicted_pct",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        for fit in fits:
            metrics = fit["metrics"]
            for key in PANEL_KEYS:
                panel = panels[key]
                for function, x, observed, predicted, ri in (
                    ("forgetting", panel.lag, panel.forgetting,
                     metrics["panel"][key]["forgetting_prediction"], ""),
                    ("spacing", panel.isi, panel.spacing,
                     metrics["panel"][key]["spacing_prediction"], panel.ri),
                ):
                    for x_i, observed_i, predicted_i in zip(x, observed, predicted):
                        writer.writerow({
                            "protocol": fit["protocol"],
                            "experiment1_m1": fit["experiment1_m1"],
                            "panel": key,
                            "function": function,
                            "isi_days": x_i,
                            "ri_days": ri,
                            "observed_pct": 100.0 * observed_i,
                            "predicted_pct": 100.0 * predicted_i,
                        })
    return path


def select_central(fits: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    """Return the m1=4 fits used for the full diagnostic figures.

    Event count is an assumed sensitivity variable rather than another fitted
    model parameter, so selecting it by fit would overstate what these aggregate
    observations identify.
    """
    central: dict[str, dict[str, object]] = {}
    for protocol in ("forgetting_calibrated", "joint", "spacing_only"):
        central[protocol] = next(
            fit for fit in fits
            if fit["protocol"] == protocol and fit["experiment1_m1"] == 4
        )
    return central


def style_axis(axis: plt.Axes) -> None:
    axis.set_ylim(-2.0, 102.0)
    axis.set_yticks([0, 20, 40, 60, 80, 100])
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.8, alpha=0.75)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(labelsize=8.5)


def add_x_padding(axis: plt.Axes, xmax: float) -> None:
    axis.set_xlim(-0.035 * xmax, 1.04 * xmax)


def normalize_generated_svg(path: Path) -> None:
    """Remove Matplotlib's path-line trailing spaces for clean Git diffs."""
    content = path.read_text(encoding="utf-8")
    path.write_text(
        "\n".join(line.rstrip() for line in content.splitlines()) + "\n",
        encoding="utf-8",
    )


def plot_full_fit(fit: dict[str, object], panels: dict[str, Panel]) -> Path:
    protocol = fit["protocol"]
    experiment1_m1 = fit["experiment1_m1"]
    pars = fit["parameters"]
    metrics = fit["metrics"]
    colors = {"a": "#0072B2", "b": "#009E73", "c": "#CC79A7"}

    figure = plt.figure(figsize=(15.3, 8.4), layout="constrained")
    grid = figure.add_gridspec(2, 4, width_ratios=(1, 1, 1, 0.60))
    forgetting_axes = [figure.add_subplot(grid[0, j]) for j in range(3)]
    spacing_axes = [figure.add_subplot(grid[1, j]) for j in range(3)]
    info_axis = figure.add_subplot(grid[:, 3])
    info_axis.set_axis_off()

    protocol_label = PROTOCOL_LABELS[protocol]
    figure.suptitle(
        f"Session-batch SAC: {protocol_label}\n"
        f"Experiment 1 has {experiment1_m1} Session-1 events; "
        "Experiment 2 has 4",
        fontsize=15,
        fontweight="semibold",
    )

    for j, key in enumerate(PANEL_KEYS):
        panel = panels[key]
        p = pars["panel"][key]
        m1 = panel_m1(key, experiment1_m1)
        color = colors[key]

        axis = forgetting_axes[j]
        xmax = float(panel.lag.max())
        xgrid = np.linspace(SHORTEST_ISI_DAYS[key], xmax, 600)
        curve = response_probability(
            forgetting_strength(xgrid, p["d"], p["tau"], pars["delta"], m1),
            pars["theta"], pars["sigma"],
        )
        axis.plot(xgrid, 100.0 * curve, color=color, linewidth=2.3)
        axis.plot(
            panel.lag, 100.0 * panel.forgetting, linestyle="none", marker="s",
            markersize=6.3, markerfacecolor="white", markeredgecolor=color,
            markeredgewidth=1.5, zorder=3,
        )
        style_axis(axis)
        add_x_padding(axis, xmax)
        axis.set_title(PANEL_LABELS[key], fontsize=11.3, pad=7)
        axis.set_xlabel("Study-test lag (days)", fontsize=9.5)
        if j == 0:
            axis.set_ylabel("Recall (%)", fontsize=10)
        axis.text(
            0.98, 0.96,
            f"RMSE = {metrics['panel'][key]['forgetting_rmse_pp']:.2f} pp",
            transform=axis.transAxes, ha="right", va="top", fontsize=8.5,
        )

        axis = spacing_axes[j]
        optimum = metrics["panel"][key]["optimum_isi_days"]
        xmax = max(float(panel.isi.max()), 1.08 * optimum)
        xgrid = np.linspace(SHORTEST_ISI_DAYS[key], xmax, 900)
        curve = response_probability(
            spacing_strength(xgrid, panel.ri, p["d"], p["tau"], pars["delta"], m1),
            pars["theta"], pars["sigma"],
        )
        axis.plot(xgrid, 100.0 * curve, color=color, linewidth=2.3)
        axis.plot(
            panel.isi, 100.0 * panel.spacing, linestyle="none", marker="o",
            markersize=6.3, markerfacecolor=color, markeredgecolor="#202020",
            markeredgewidth=0.6, zorder=3,
        )
        axis.axvline(
            optimum, color=color, linewidth=1.25, alpha=0.7,
            linestyle=(0, (1.5, 2.5)),
        )
        style_axis(axis)
        add_x_padding(axis, xmax)
        axis.set_title(f"Spacing curve, RI = {panel.ri:g} days", fontsize=10.3, pad=7)
        axis.set_xlabel("ISI (days)", fontsize=9.5)
        if j == 0:
            axis.set_ylabel("Final-test recall (%)", fontsize=10)
        axis.text(
            0.98, 0.96,
            f"RMSE = {metrics['panel'][key]['spacing_rmse_pp']:.2f} pp\n"
            f"optimum = {optimum:.2f} d",
            transform=axis.transAxes, ha="right", va="top", fontsize=8.5,
        )

    parameter_lines = [
        "Shared parameters",
        f"delta = {pars['delta']:.4f}",
        f"theta = {pars['theta']:.4f}",
        f"sigma = {pars['sigma']:.4f}",
        "",
        "Experiment-specific decay",
    ]
    for key in PANEL_KEYS:
        p = pars["panel"][key]
        parameter_lines.extend([
            f"{key}: d = {p['d']:.4f}",
            f"   tau = {p['tau']:.4g} d",
        ])
    parameter_lines.extend([
        "",
        "Overall RMSE",
        f"forgetting = {metrics['forgetting_rmse_pp']:.2f} pp",
        f"spacing = {metrics['spacing_rmse_pp']:.2f} pp",
        f"joint = {metrics['joint_rmse_pp']:.2f} pp",
        "",
        "Lines: model",
        "Points: observations",
        "Dotted lines: continuous optima",
        "",
        "Zero-day gaps corrected to",
        "5 min (Exp. 1) and 20 min",
        "(Experiments 2a and 2b).",
    ])
    info_axis.text(
        0.0, 0.98, "\n".join(parameter_lines), transform=info_axis.transAxes,
        ha="left", va="top", fontsize=9.2, linespacing=1.32,
    )

    stem = FIGURES_DIR / f"sac_cepeda2009_batch_{protocol}_m1_{experiment1_m1}"
    figure.savefig(stem.with_suffix(".png"), dpi=220, facecolor="white")
    figure.savefig(stem.with_suffix(".svg"), facecolor="white")
    normalize_generated_svg(stem.with_suffix(".svg"))
    plt.close(figure)
    return stem


def plot_sensitivity(fits: list[dict[str, object]]) -> Path:
    labels = {
        "forgetting_calibrated": "Forgetting + Exp. 1 calibration",
        "joint": "Joint",
        "spacing_only": "Spacing only",
    }
    colors = {
        "forgetting_calibrated": "#0072B2",
        "joint": "#D55E00",
        "spacing_only": "#009E73",
    }
    figure, axes = plt.subplots(1, 4, figsize=(15.5, 4.3), layout="constrained")
    for protocol in labels:
        subset = sorted(
            (fit for fit in fits if fit["protocol"] == protocol),
            key=lambda fit: fit["experiment1_m1"],
        )
        x = [fit["experiment1_m1"] for fit in subset]
        axes[0].plot(
            x, [fit["metrics"]["forgetting_rmse_pp"] for fit in subset],
            marker="o", linewidth=2, label=labels[protocol], color=colors[protocol],
        )
        axes[1].plot(
            x, [fit["metrics"]["spacing_rmse_pp"] for fit in subset],
            marker="o", linewidth=2, label=labels[protocol], color=colors[protocol],
        )
        axes[2].plot(
            x, [fit["parameters"]["delta"] for fit in subset],
            marker="o", linewidth=2, label=labels[protocol], color=colors[protocol],
        )
        axes[3].plot(
            x, [fit["metrics"]["panel"]["a"]["optimum_isi_days"] for fit in subset],
            marker="o", linewidth=2, label=labels[protocol], color=colors[protocol],
        )

    titles = (
        "Forgetting fit",
        "Spacing fit",
        "Learning rate",
        "Experiment 1 optimum",
    )
    ylabels = ("RMSE (percentage points)", "RMSE (percentage points)", "delta", "ISI (days)")
    for axis, title, ylabel in zip(axes, titles, ylabels):
        axis.set_title(title, fontsize=11)
        axis.set_xlabel("Assumed Exp. 1 Session-1 events")
        axis.set_ylabel(ylabel)
        axis.set_xticks([3, 4, 5, 6, 7])
        axis.grid(color="#D9D9D9", linewidth=0.8, alpha=0.75)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8.8)
    figure.suptitle(
        "Session-batch SAC sensitivity to the unknown Experiment-1 event count",
        fontsize=14, fontweight="semibold",
    )
    stem = FIGURES_DIR / "sac_cepeda2009_batch_sensitivity"
    figure.savefig(stem.with_suffix(".png"), dpi=220, facecolor="white")
    figure.savefig(stem.with_suffix(".svg"), facecolor="white")
    normalize_generated_svg(stem.with_suffix(".svg"))
    plt.close(figure)
    return stem


def write_report(fits: list[dict[str, object]], central: dict[str, dict[str, object]]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2009_batch_report.md"
    lines = [
        "# Session-batch SAC fits to Cepeda et al. (2009)",
        "",
        "## Model and approximation",
        "",
        "All learning events within a session are collapsed to the same time, but their nonlinear SAC updates are retained. For a session containing $m$ events,",
        "",
        "$$",
        "A_m = 1-(1-\\delta)^m.",
        "$$",
        "",
        "The forgetting and spacing strengths are",
        "",
        "$$",
        "B_F(a)=A_{m_1}f(a),",
        "$$",
        "",
        "$$",
        "B_S(a,b)=A_{m_1}f(a+b)+A_{m_2}[1-A_{m_1}f(a)]f(b),",
        "$$",
        "",
        "with $f(t)=(1+t/\\tau)^{-d}$ and a shared logistic response mapping. Experiment 2 uses $m_1=4$ and all experiments use $m_2=2$. Because the number of Session-1 trials in Experiment 1 is unavailable, $m_1=3,\\ldots,7$ is treated as a sensitivity analysis. The nominal zero-day gaps are represented as 5 minutes in Experiment 1 and 20 minutes in Experiment 2.",
        "",
        "Separate $d$ and $\\tau$ parameters are estimated for the three curves; $\\delta$, $\\theta$, and $\\sigma$ are shared.",
        "",
        "The forgetting-trained protocol uses a bilevel fit: for each candidate $\\delta$, all remaining parameters are optimized using only the forgetting observations; $\\delta$ is then calibrated using only the Experiment-1 spacing curve. The joint and spacing-only protocols estimate every parameter from their named datasets.",
        "",
        "## Central sensitivity assumption: Experiment 1 $m_1=4$",
        "",
        "| Protocol | Exp. 1 $m_1$ | $\\delta$ | Forgetting RMSE | Spacing RMSE | Exp. 1 optimum |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for protocol in ("forgetting_calibrated", "joint", "spacing_only"):
        fit = central[protocol]
        metrics = fit["metrics"]
        lines.append(
            f"| {PROTOCOL_LABELS[protocol]} | {fit['experiment1_m1']} | "
            f"{fit['parameters']['delta']:.4f} | {metrics['forgetting_rmse_pp']:.2f} pp | "
            f"{metrics['spacing_rmse_pp']:.2f} pp | "
            f"{metrics['panel']['a']['optimum_isi_days']:.2f} d |"
        )
    lines.extend([
        "",
        "## Change from the single-event approximation",
        "",
        "The table compares the central $m_1=4$ batch model with the preceding single-event/normalized-interaction analyses. Parameter counts are unchanged within each protocol.",
        "",
        "| Protocol | Spacing RMSE, single event | Spacing RMSE, batch | Exp. 1 optimum, single event | Exp. 1 optimum, batch |",
        "|---|---:|---:|---:|---:|",
        f"| Forgetting fit + Exp. 1 calibration | 7.93 pp | {central['forgetting_calibrated']['metrics']['spacing_rmse_pp']:.2f} pp | 12.55 d | {central['forgetting_calibrated']['metrics']['panel']['a']['optimum_isi_days']:.2f} d |",
        f"| Joint fit | 4.72 pp | {central['joint']['metrics']['spacing_rmse_pp']:.2f} pp | 6.53 d | {central['joint']['metrics']['panel']['a']['optimum_isi_days']:.2f} d |",
        f"| Spacing-only fit | 2.61 pp | {central['spacing_only']['metrics']['spacing_rmse_pp']:.2f} pp | 1.53 d | {central['spacing_only']['metrics']['panel']['a']['optimum_isi_days']:.2f} d |",
        "",
        "The repeated-event correction therefore removes most of the extreme optimum displacement in the forgetting-trained analysis and improves the joint fit moderately. It does not fully reconcile the forgetting and spacing constraints: the spacing-only fit still prefers an Experiment-1 optimum near 1.7 days, whereas the joint fit prefers about 5.2 days. The spacing-only Experiment-1 time scale also remains at the numerical lower bound, so its excellent spacing fit does not supply a stable estimate of the underlying forgetting time scale.",
        "",
        "## Complete sensitivity results",
        "",
        "| Protocol | Exp. 1 $m_1$ | $\\delta$ | Forgetting RMSE | Spacing RMSE | Joint RMSE | Exp. 1 optimum |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for fit in sorted(fits, key=lambda item: (item["protocol"], item["experiment1_m1"])):
        metrics = fit["metrics"]
        lines.append(
            f"| {PROTOCOL_LABELS[fit['protocol']]} | {fit['experiment1_m1']} | "
            f"{fit['parameters']['delta']:.4f} | {metrics['forgetting_rmse_pp']:.2f} pp | "
            f"{metrics['spacing_rmse_pp']:.2f} pp | {metrics['joint_rmse_pp']:.2f} pp | "
            f"{metrics['panel']['a']['optimum_isi_days']:.2f} d |"
        )
    lines.extend([
        "",
        "![Sensitivity summary](../figures/sac_cepeda2009_batch_sensitivity.svg)",
        "",
        "## Full fits for the central $m_1=4$ assumption",
        "",
    ])
    for protocol in ("forgetting_calibrated", "joint", "spacing_only"):
        m1 = central[protocol]["experiment1_m1"]
        label = PROTOCOL_LABELS[protocol]
        lines.extend([
            f"### {label}",
            "",
            f"![{label}](../figures/sac_cepeda2009_batch_{protocol}_m1_{m1}.svg)",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def print_summary(fits: list[dict[str, object]], central: dict[str, dict[str, object]]) -> None:
    for fit in fits:
        pars = fit["parameters"]
        metrics = fit["metrics"]
        print(
            fit["protocol"],
            "m1", fit["experiment1_m1"],
            "delta", f"{pars['delta']:.6f}",
            "forget_rmse", f"{metrics['forgetting_rmse_pp']:.4f}",
            "spacing_rmse", f"{metrics['spacing_rmse_pp']:.4f}",
            "joint_rmse", f"{metrics['joint_rmse_pp']:.4f}",
            "optimum_a", f"{metrics['panel']['a']['optimum_isi_days']:.4f}",
            "condition", f"{fit['jacobian_condition']:.3e}",
            flush=True,
        )
    print("\nCENTRAL M1=4", flush=True)
    for protocol, fit in central.items():
        print(protocol, fit_row(fit), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=96)
    parser.add_argument("--seed", type=int, default=20260825)
    parser.add_argument("--calibration-grid", type=int, default=25)
    args = parser.parse_args()

    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    panels = load_panels()
    fits: list[dict[str, object]] = []
    for protocol_index, protocol in enumerate(
        ("forgetting_calibrated", "joint", "spacing_only")
    ):
        for experiment1_m1 in range(3, 8):
            print(f"Fitting {protocol}, Experiment-1 m1={experiment1_m1}", flush=True)
            if protocol == "forgetting_calibrated":
                fit = fit_forgetting_calibrated(
                    panels,
                    experiment1_m1,
                    starts=args.starts,
                    seed=args.seed + experiment1_m1,
                    grid_size=args.calibration_grid,
                )
            else:
                fit = fit_model(
                    panels,
                    experiment1_m1,
                    protocol,
                    starts=args.starts,
                    seed=args.seed + 1000 * protocol_index + experiment1_m1,
                )
            add_metrics(fit, panels)
            fits.append(fit)

    central = select_central(fits)
    fit_path = write_fit_table(fits)
    prediction_path = write_prediction_table(fits, panels)
    sensitivity_stem = plot_sensitivity(fits)
    for fit in central.values():
        plot_full_fit(fit, panels)
    report_path = write_report(fits, central)
    print_summary(fits, central)
    print("\nWROTE", fit_path, prediction_path, sensitivity_stem, report_path, flush=True)


if __name__ == "__main__":
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.9,
        "savefig.bbox": "tight",
        "svg.fonttype": "none",
    })
    main()

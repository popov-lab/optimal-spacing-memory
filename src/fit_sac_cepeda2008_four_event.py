"""Fit time-resolved four-event SAC to Cepeda et al. (2008).

The primary model follows the approximate event schedule supplied for the
experiment (all times are converted to days):

    event 1 = 0
    event 2 = 320 seconds
    event 3 = 320 seconds + ISI
    event 4 = 640 seconds + ISI
    final test = 960 seconds + ISI + RI

Each event uses the general SAC recursion

    u_n = delta * (1 - B(t_n))
    B(t) = sum_k u_k * (1 + (t - t_k) / tau)**(-d).

The script fits the 26 final-test spacing observations only. It fits both a
free-delta model and a delta=1 model. For a like-for-like assessment of the
effect of representing within-session repetitions, it also refits the prior
two-event schedule using the same parameterization and optimizer.

Run from the repository root with:

    python3 src/fit_sac_cepeda2008_four_event.py
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import expit, logit


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "cepeda_spacing_recall.csv"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

SECONDS_PER_DAY = 24.0 * 60.0 * 60.0
WITHIN_SESSION_GAP_DAYS = 320.0 / SECONDS_PER_DAY
ZERO_ISI_DAYS = 0.00256
RI_LEVELS = (7.0, 35.0, 70.0, 350.0)
SCHEDULES = ("two_event", "four_event")
DELTA_VARIANTS = ("free", "fixed_one")
SCHEDULE_LABELS = {
    "two_event": "Two-event baseline",
    "four_event": "Four-event timing model",
}
DELTA_LABELS = {
    "free": "Free delta",
    "fixed_one": "delta = 1",
}


@dataclass(frozen=True)
class SpacingData:
    rows: tuple[dict[str, str], ...]
    observed_isi: np.ndarray
    model_isi: np.ndarray
    ri: np.ndarray
    recall: np.ndarray


def load_spacing_data(path: Path = DATA_PATH) -> SpacingData:
    """Load only the 26 Cepeda (2008) final spacing-test observations."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = tuple(
            row
            for row in csv.DictReader(handle)
            if row["panel"] == "d" and row["function"] == "spacing"
        )
    rows = tuple(sorted(
        rows,
        key=lambda row: (float(row["ri_days"]), float(row["isi_days"])),
    ))
    observed_isi = np.array([float(row["isi_days"]) for row in rows])
    model_isi = observed_isi.copy()
    model_isi[model_isi == 0.0] = ZERO_ISI_DAYS
    return SpacingData(
        rows=rows,
        observed_isi=observed_isi,
        model_isi=model_isi,
        ri=np.array([float(row["ri_days"]) for row in rows]),
        recall=np.array([float(row["recall_pct"]) / 100.0 for row in rows]),
    )


def _forgetting_from_log_tau(
    lag: np.ndarray | float, d: float, log_tau: float
) -> np.ndarray:
    """Evaluate shifted-power forgetting stably for arbitrarily small tau."""
    lag_array = np.asarray(lag, dtype=float)
    if np.any(lag_array < 0.0):
        raise ValueError("Forgetting lags must be nonnegative")
    log_ratio = np.full_like(lag_array, -np.inf, dtype=float)
    positive = lag_array > 0.0
    log_ratio[positive] = np.log(lag_array[positive]) - log_tau
    log_one_plus_ratio = np.logaddexp(0.0, log_ratio)
    return np.exp(-d * log_one_plus_ratio)


def sac_forgetting(
    lag: np.ndarray | float, d: float, tau: float
) -> np.ndarray:
    """Return ``(1 + lag / tau)**(-d)`` with stable small-tau arithmetic."""
    if d < 0.0:
        raise ValueError("d must be nonnegative")
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    return _forgetting_from_log_tau(lag, d, float(np.log(tau)))


def _sac_strength_from_log_tau(
    event_times: np.ndarray,
    test_time: float,
    delta: float,
    d: float,
    log_tau: float,
) -> float:
    """Apply the general SAC recursion to one ordered event schedule."""
    times = np.asarray(event_times, dtype=float)
    if times.ndim != 1 or times.size == 0:
        raise ValueError("event_times must be a nonempty one-dimensional array")
    if np.any(np.diff(times) < 0.0):
        raise ValueError("event_times must be ordered")
    if test_time < times[-1]:
        raise ValueError("test_time must not precede the final event")
    if not 0.0 < delta <= 1.0:
        raise ValueError("delta must be in (0, 1]")

    increments = np.empty(times.size, dtype=float)
    for event_index, event_time in enumerate(times):
        if event_index == 0:
            strength_before_event = 0.0
        else:
            surviving = _forgetting_from_log_tau(
                event_time - times[:event_index], d, log_tau
            )
            strength_before_event = float(
                increments[:event_index] @ surviving
            )
        increments[event_index] = delta * (1.0 - strength_before_event)

    test_survival = _forgetting_from_log_tau(test_time - times, d, log_tau)
    return float(increments @ test_survival)


def sac_strength_from_times(
    event_times: np.ndarray,
    test_time: float,
    delta: float,
    d: float,
    tau: float,
) -> float:
    """Public tau-parameterized wrapper for the general SAC recursion."""
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    return _sac_strength_from_log_tau(
        event_times, test_time, delta, d, float(np.log(tau))
    )


def two_event_times(isi: float, ri: float) -> tuple[np.ndarray, float]:
    """Return the prior abstraction: studies at 0 and ISI, then test after RI."""
    if isi < 0.0 or ri < 0.0:
        raise ValueError("ISI and RI must be nonnegative")
    return np.array([0.0, isi]), isi + ri


def four_event_times(isi: float, ri: float) -> tuple[np.ndarray, float]:
    """Return the supplied four-event Cepeda (2008) approximate schedule."""
    if isi < 0.0 or ri < 0.0:
        raise ValueError("ISI and RI must be nonnegative")
    within = WITHIN_SESSION_GAP_DAYS
    events = np.array([0.0, within, within + isi, 2.0 * within + isi])
    test_time = 3.0 * within + isi + ri
    return events, test_time


def _schedule_function(
    schedule: str,
) -> Callable[[float, float], tuple[np.ndarray, float]]:
    if schedule == "two_event":
        return two_event_times
    if schedule == "four_event":
        return four_event_times
    raise ValueError(f"Unknown schedule: {schedule}")


def _sac_strength_for_schedule_log_tau(
    isi: np.ndarray | float,
    ri: np.ndarray | float,
    delta: float,
    d: float,
    log_tau: float,
    schedule: str,
) -> np.ndarray:
    isi_array, ri_array = np.broadcast_arrays(
        np.asarray(isi, dtype=float), np.asarray(ri, dtype=float)
    )
    build_times = _schedule_function(schedule)
    output = np.empty(isi_array.size, dtype=float)
    for index, (isi_value, ri_value) in enumerate(
        zip(isi_array.ravel(), ri_array.ravel())
    ):
        events, test_time = build_times(float(isi_value), float(ri_value))
        output[index] = _sac_strength_from_log_tau(
            events, test_time, delta, d, log_tau
        )
    return output.reshape(isi_array.shape)


def sac_strength_for_schedule(
    isi: np.ndarray | float,
    ri: np.ndarray | float,
    delta: float,
    d: float,
    tau: float,
    schedule: str = "four_event",
) -> np.ndarray:
    """Return SAC final-test strength for one of the fitted schedules."""
    if tau <= 0.0:
        raise ValueError("tau must be positive")
    return _sac_strength_for_schedule_log_tau(
        isi, ri, delta, d, float(np.log(tau)), schedule
    )


def response_probability(
    strength: np.ndarray | float, theta: float, sigma: float
) -> np.ndarray:
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    return expit((np.asarray(strength, dtype=float) - theta) / sigma)


def _parameter_bounds(delta_variant: str) -> tuple[np.ndarray, np.ndarray]:
    # log(tau) is intentionally unbounded. The bounds on other parameters
    # exclude only numerically degenerate regions, as in the earlier fit.
    lower = [np.log(1e-4), -np.inf, -2.0, np.log(1e-3)]
    upper = [np.log(10.0), np.inf, 3.0, np.log(5.0)]
    if delta_variant == "free":
        lower.insert(0, logit(1e-4))
        upper.insert(0, logit(0.9999))
    elif delta_variant != "fixed_one":
        raise ValueError(f"Unknown delta variant: {delta_variant}")
    return np.asarray(lower), np.asarray(upper)


def _unpack_parameters(z: np.ndarray, delta_variant: str) -> dict[str, float]:
    if delta_variant == "free":
        logit_delta, log_d, log_tau, theta, log_sigma = z
        delta = float(expit(logit_delta))
    else:
        log_d, log_tau, theta, log_sigma = z
        delta = 1.0
    with np.errstate(over="ignore", under="ignore"):
        tau = float(np.exp(log_tau))
    return {
        "delta": delta,
        "d": float(np.exp(log_d)),
        "log_tau": float(log_tau),
        "tau": tau,
        "theta": float(theta),
        "sigma": float(np.exp(log_sigma)),
    }


def _initial_values(
    delta_variant: str, starts: int, seed: int
) -> list[np.ndarray]:
    lower, upper = _parameter_bounds(delta_variant)
    initials: list[np.ndarray] = []
    decay_starts = (
        (0.05, 1e-12),
        (0.08, 1e-8),
        (0.12, 1e-4),
        (0.14, 0.03),
        (0.30, 1.0),
        (0.80, 30.0),
    )
    delta_starts = (0.10, 0.25, 0.40, 0.70, 0.95)
    for delta in delta_starts if delta_variant == "free" else (1.0,):
        for d, tau in decay_starts:
            values = [np.log(d), np.log(tau), 0.25, np.log(0.05)]
            if delta_variant == "free":
                values.insert(0, logit(delta))
            initials.append(np.asarray(values))

    rng = np.random.default_rng(seed)
    practical_lower = [np.log(0.005), np.log(1e-16), -0.5, np.log(0.003)]
    practical_upper = [np.log(2.0), np.log(1e3), 1.5, np.log(0.8)]
    if delta_variant == "free":
        practical_lower.insert(0, logit(0.01))
        practical_upper.insert(0, logit(0.99))
    practical_lower_array = np.asarray(practical_lower)
    practical_upper_array = np.asarray(practical_upper)
    while len(initials) < starts:
        z = practical_lower_array + rng.random(len(lower)) * (
            practical_upper_array - practical_lower_array
        )
        initials.append(np.minimum(np.maximum(z, lower), upper))
    return initials[:starts]


def fit_model(
    data: SpacingData,
    schedule: str,
    delta_variant: str,
    starts: int = 128,
    seed: int = 20260825,
) -> dict[str, object]:
    """Fit one schedule and delta variant by multistart least squares."""
    lower, upper = _parameter_bounds(delta_variant)

    def residual(z: np.ndarray) -> np.ndarray:
        parameters = _unpack_parameters(z, delta_variant)
        strength = _sac_strength_for_schedule_log_tau(
            data.model_isi,
            data.ri,
            parameters["delta"],
            parameters["d"],
            parameters["log_tau"],
            schedule,
        )
        return response_probability(
            strength, parameters["theta"], parameters["sigma"]
        ) - data.recall

    solutions: list[tuple[float, object]] = []
    for z0 in _initial_values(delta_variant, starts, seed):
        fit = least_squares(
            residual,
            z0,
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
        raise RuntimeError(
            f"No finite fit for schedule={schedule}, delta={delta_variant}"
        )
    solutions.sort(key=lambda item: item[0])
    sse, best = solutions[0]
    parameters = _unpack_parameters(best.x, delta_variant)
    prediction = data.recall + best.fun
    singular_values = np.linalg.svd(best.jac, compute_uv=False)
    condition = (
        np.inf
        if singular_values[-1] == 0.0
        else float(singular_values[0] / singular_values[-1])
    )
    return {
        "schedule": schedule,
        "delta_variant": delta_variant,
        "parameter_count": len(best.x),
        "parameters": parameters,
        "prediction": prediction,
        "sse": sse,
        "rmse_pp": 100.0 * float(np.sqrt(sse / data.recall.size)),
        "jacobian_condition": condition,
        "optimality": float(best.optimality),
        "nfev": int(best.nfev),
        "best_sse_values": [float(item[0]) for item in solutions[:10]],
    }


def optimum_isi(
    fit: dict[str, object], ri: float, upper: float = 105.0
) -> float:
    """Return the maximum-strength ISI within the observed 0--105 day range."""
    parameters = fit["parameters"]

    def negative_strength(isi: float) -> float:
        model_isi = ZERO_ISI_DAYS if isi == 0.0 else isi
        strength = _sac_strength_for_schedule_log_tau(
            model_isi,
            ri,
            parameters["delta"],
            parameters["d"],
            parameters["log_tau"],
            fit["schedule"],
        )
        return -float(strength)

    optimized = minimize_scalar(
        negative_strength,
        bounds=(ZERO_ISI_DAYS, upper),
        method="bounded",
        options={"xatol": 1e-8},
    )
    candidates = (ZERO_ISI_DAYS, float(optimized.x), upper)
    return min(candidates, key=negative_strength)


def _fit_row(fit: dict[str, object]) -> dict[str, object]:
    parameters = fit["parameters"]
    row: dict[str, object] = {
        "schedule": fit["schedule"],
        "delta_variant": fit["delta_variant"],
        "parameters": fit["parameter_count"],
        "delta": parameters["delta"],
        "d": parameters["d"],
        "tau_days": parameters["tau"],
        "log_tau_days": parameters["log_tau"],
        "theta": parameters["theta"],
        "sigma": parameters["sigma"],
        "rmse_pp": fit["rmse_pp"],
        "sse": fit["sse"],
        "jacobian_condition": fit["jacobian_condition"],
        "optimality": fit["optimality"],
        "nfev": fit["nfev"],
    }
    for ri in RI_LEVELS:
        row[f"optimum_isi_ri_{int(ri)}_days"] = optimum_isi(fit, ri)
    return row


def write_fit_table(fits: list[dict[str, object]]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2008_four_event_fits.csv"
    rows = [_fit_row(fit) for fit in fits]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_prediction_table(
    fits: list[dict[str, object]], data: SpacingData
) -> Path:
    path = RESULTS_DIR / "sac_cepeda2008_four_event_predictions.csv"
    fieldnames = [
        "schedule",
        "delta_variant",
        "isi_days",
        "model_isi_days",
        "ri_days",
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
            for isi, model_isi, ri, observed, predicted in zip(
                data.observed_isi,
                data.model_isi,
                data.ri,
                data.recall,
                fit["prediction"],
            ):
                writer.writerow({
                    "schedule": fit["schedule"],
                    "delta_variant": fit["delta_variant"],
                    "isi_days": isi,
                    "model_isi_days": model_isi,
                    "ri_days": ri,
                    "observed_pct": 100.0 * observed,
                    "predicted_pct": 100.0 * predicted,
                    "residual_pp": 100.0 * (predicted - observed),
                })
    return path


def normalize_generated_svg(path: Path) -> None:
    content = path.read_text(encoding="utf-8")
    path.write_text(
        "\n".join(line.rstrip() for line in content.splitlines()) + "\n",
        encoding="utf-8",
    )


def plot_four_event_fits(
    fits: list[dict[str, object]], data: SpacingData
) -> Path:
    selected = {
        fit["delta_variant"]: fit
        for fit in fits
        if fit["schedule"] == "four_event"
    }
    colors = {7.0: "#0072B2", 35.0: "#009E73", 70.0: "#D55E00", 350.0: "#CC79A7"}
    figure = plt.figure(figsize=(13.3, 5.4), layout="constrained")
    grid = figure.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 0.48))
    axes = [figure.add_subplot(grid[0, index]) for index in range(2)]
    info_axis = figure.add_subplot(grid[0, 2])
    info_axis.set_axis_off()

    for axis, delta_variant in zip(axes, DELTA_VARIANTS):
        fit = selected[delta_variant]
        parameters = fit["parameters"]
        for ri in RI_LEVELS:
            subset = data.ri == ri
            xgrid = np.linspace(ZERO_ISI_DAYS, 105.0, 800)
            strength = _sac_strength_for_schedule_log_tau(
                xgrid,
                np.full_like(xgrid, ri),
                parameters["delta"],
                parameters["d"],
                parameters["log_tau"],
                "four_event",
            )
            curve = response_probability(
                strength, parameters["theta"], parameters["sigma"]
            )
            color = colors[ri]
            axis.plot(
                xgrid, 100.0 * curve, color=color, linewidth=2.0,
                label=f"RI = {ri:g} d",
            )
            axis.plot(
                data.observed_isi[subset],
                100.0 * data.recall[subset],
                linestyle="none",
                marker="o",
                markersize=5.7,
                markerfacecolor="white",
                markeredgecolor=color,
                markeredgewidth=1.35,
            )
        axis.set_xlim(-4.0, 109.0)
        axis.set_ylim(-2.0, 102.0)
        axis.set_xticks([0, 7, 14, 21, 35, 70, 105])
        axis.set_yticks([0, 20, 40, 60, 80, 100])
        axis.set_xlabel("ISI between events 2 and 3 (days)")
        axis.set_title(
            f"{DELTA_LABELS[delta_variant]}\nRMSE = {fit['rmse_pp']:.2f} pp",
            fontsize=11.5,
        )
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.8, alpha=0.75)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Final-test recall (%)")
    axes[0].legend(frameon=False, fontsize=8.8, loc="lower right")

    info_lines = [
        "Four-event schedule",
        "1: 0 s",
        "2: 320 s",
        "3: 320 s + ISI",
        "4: 640 s + ISI",
        "test: 960 s + ISI + RI",
        "",
        "Lines: fitted model",
        "Points: observations",
        "",
        "Spacing observations only",
        "Nominal zero ISI: 3.69 min",
    ]
    info_axis.text(
        0.0,
        0.98,
        "\n".join(info_lines),
        transform=info_axis.transAxes,
        ha="left",
        va="top",
        fontsize=9.2,
        linespacing=1.35,
    )
    figure.suptitle(
        "Cepeda et al. (2008): time-resolved four-event SAC",
        fontsize=14.5,
        fontweight="semibold",
    )

    path = FIGURES_DIR / "sac_cepeda2008_four_event_spacing.svg"
    figure.savefig(path, facecolor="white", metadata={"Date": None})
    normalize_generated_svg(path)
    plt.close(figure)
    return path


def _format_tau(tau: float) -> str:
    if tau < 0.001 or tau >= 1000.0:
        return f"{tau:.3e}"
    return f"{tau:.5f}"


def write_report(fits: list[dict[str, object]]) -> Path:
    path = RESULTS_DIR / "sac_cepeda2008_four_event_report.md"
    by_key = {
        (fit["schedule"], fit["delta_variant"]): fit for fit in fits
    }
    lines = [
        "# Time-resolved SAC fits to Cepeda et al. (2008)",
        "",
        "## Scope",
        "",
        "These are post-hoc fits to the 26 final spacing-test observations. The newly available Session-2 forgetting observations are deliberately excluded. The nominal zero-day ISI is represented as the reported approximately 3-minute interval, $0.00256$ days.",
        "",
        "## Four-event schedule and model",
        "",
        "The timing approximation uses $h=320/86400=0.0037037$ days and",
        "",
        "| Point | Absolute time |",
        "|---|---:|",
        "| Event 1 (Session 1) | $0$ |",
        "| Event 2 (Session 1) | $h$ |",
        "| Event 3 (Session 2) | $h+\\mathrm{ISI}$ |",
        "| Event 4 (Session 2) | $2h+\\mathrm{ISI}$ |",
        "| Final test | $3h+\\mathrm{ISI}+\\mathrm{RI}$ |",
        "",
        "Thus the manipulated RI is the variable part of the post-event-4 delay; the supplied absolute-time schedule also places a fixed 320 seconds between event 4 and the RI anchor.",
        "",
        "At every event, the implementation evaluates the general SAC recursion directly:",
        "",
        "$$",
        "u_n=\\delta[1-B(t_n)],",
        "\\qquad",
        "B(t)=\\sum_{k:t_k<t}u_k f(t-t_k),",
        "$$",
        "",
        "with",
        "",
        "$$",
        "f(t)=\\left(1+\\frac{t}{\\tau}\\right)^{-d}.",
        "$$",
        "",
        "Latent strength is mapped to recall probability with the same logistic rule as the preceding two-event analysis:",
        "",
        "$$",
        "P(\\mathrm{recall})=\\frac{1}{1+\\exp[-(B-\\theta)/\\sigma]}.",
        "$$",
        "",
        "The free-$\\delta$ variant estimates $\\{\\delta,d,\\tau,\\theta,\\sigma\\}$; the fixed variant sets $\\delta=1$ and estimates $\\{d,\\tau,\\theta,\\sigma\\}$. The optimization is performed in $\\log\\tau$ with no lower or upper bound on that coordinate.",
        "",
        "## Results",
        "",
        "The two-event baseline was refit with the same code and unconstrained $\\log\\tau$, so the change due to the event schedule is not confounded with the earlier optimizer bounds.",
        "",
        "| Event representation | Learning-rate variant | Parameters | $\\delta$ | $d$ | $\\tau$ (days) | $\\theta$ | $\\sigma$ | RMSE |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for schedule in SCHEDULES:
        for delta_variant in DELTA_VARIANTS:
            fit = by_key[(schedule, delta_variant)]
            parameters = fit["parameters"]
            lines.append(
                f"| {SCHEDULE_LABELS[schedule]} | {DELTA_LABELS[delta_variant]} | "
                f"{fit['parameter_count']} | {parameters['delta']:.5f} | "
                f"{parameters['d']:.5f} | {_format_tau(parameters['tau'])} | "
                f"{parameters['theta']:.5f} | {parameters['sigma']:.5f} | "
                f"{fit['rmse_pp']:.2f} pp |"
            )
    lines.extend([
        "",
        "### Change produced by the four-event schedule",
        "",
        "| Learning-rate variant | Two-event RMSE | Four-event RMSE | Change |",
        "|---|---:|---:|---:|",
    ])
    for delta_variant in DELTA_VARIANTS:
        baseline = by_key[("two_event", delta_variant)]["rmse_pp"]
        four_event = by_key[("four_event", delta_variant)]["rmse_pp"]
        lines.append(
            f"| {DELTA_LABELS[delta_variant]} | {baseline:.2f} pp | "
            f"{four_event:.2f} pp | {four_event - baseline:+.2f} pp |"
        )
    free_four = by_key[("four_event", "free")]
    fixed_four = by_key[("four_event", "fixed_one")]
    free_parameters = free_four["parameters"]
    fixed_parameters = fixed_four["parameters"]
    free_prediction_change = 100.0 * float(np.max(np.abs(
        by_key[("four_event", "free")]["prediction"]
        - by_key[("two_event", "free")]["prediction"]
    )))
    fixed_prediction_change = 100.0 * float(np.max(np.abs(
        by_key[("four_event", "fixed_one")]["prediction"]
        - by_key[("two_event", "fixed_one")]["prediction"]
    )))
    lines.extend([
        "",
        f"With free $\\delta$, the four-event timing model is substantively indistinguishable from the two-event abstraction: RMSE improves by only {by_key[('two_event', 'free')]['rmse_pp'] - free_four['rmse_pp']:.4f} percentage points, and the largest change among the 26 fitted predictions is {free_prediction_change:.4f} percentage points. The extra repetitions are absorbed mainly by a change in the learning-rate estimate, from $\\delta={by_key[('two_event', 'free')]['parameters']['delta']:.4f}$ to $\\delta={free_parameters['delta']:.4f}$.",
        "",
        f"With $\\delta=1$, four events improve RMSE by {by_key[('two_event', 'fixed_one')]['rmse_pp'] - fixed_four['rmse_pp']:.2f} percentage points and change the fitted observations by at most {fixed_prediction_change:.2f} percentage points. Within the four-event representation, freeing $\\delta$ improves RMSE by only {fixed_four['rmse_pp'] - free_four['rmse_pp']:.2f} percentage points. The additional learning-rate parameter therefore has little effect on spacing-only fit quality.",
        "",
        f"The free four-event fit estimates $\\tau={_format_tau(free_parameters['tau'])}$ days (about {24.0 * 60.0 * free_parameters['tau']:.1f} minutes). By contrast, the $\\delta=1$ estimate is $\\tau={_format_tau(fixed_parameters['tau'])}$ days, effectively zero on the experimental time scale. This is an asymptotic, weakly identified scale rather than evidence for a meaningful microsecond forgetting constant; over every positive observed lag, the shifted power law is behaving like its unshifted power-law limit. The large Jacobian condition number in the fit table is consistent with that interpretation.",
        "",
        "The fitted optima below are constrained to the observed ISI range from the corrected zero gap through 105 days.",
        "",
        "| Four-event variant | RI = 7 d | RI = 35 d | RI = 70 d | RI = 350 d |",
        "|---|---:|---:|---:|---:|",
    ])
    for delta_variant in DELTA_VARIANTS:
        fit = by_key[("four_event", delta_variant)]
        optima = [optimum_isi(fit, ri) for ri in RI_LEVELS]
        lines.append(
            f"| {DELTA_LABELS[delta_variant]} | "
            + " | ".join(f"{value:.2f} d" for value in optima)
            + " |"
        )
    lines.extend([
        "",
        "![Four-event spacing fits](../figures/sac_cepeda2008_four_event_spacing.svg)",
        "",
        "Point predictions are in [`sac_cepeda2008_four_event_predictions.csv`](sac_cepeda2008_four_event_predictions.csv), and full-precision parameters and diagnostics are in [`sac_cepeda2008_four_event_fits.csv`](sac_cepeda2008_four_event_fits.csv).",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def print_summary(fits: list[dict[str, object]]) -> None:
    for fit in fits:
        parameters = fit["parameters"]
        print(
            fit["schedule"],
            fit["delta_variant"],
            f"delta={parameters['delta']:.8g}",
            f"d={parameters['d']:.8g}",
            f"tau={parameters['tau']:.8g}",
            f"theta={parameters['theta']:.8g}",
            f"sigma={parameters['sigma']:.8g}",
            f"rmse={fit['rmse_pp']:.8g}",
            f"condition={fit['jacobian_condition']:.4e}",
            flush=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    if args.starts < 1:
        parser.error("--starts must be at least 1")

    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    data = load_spacing_data()
    fits: list[dict[str, object]] = []
    for schedule_index, schedule in enumerate(SCHEDULES):
        for delta_index, delta_variant in enumerate(DELTA_VARIANTS):
            print(
                f"Fitting {schedule}, {delta_variant} delta",
                flush=True,
            )
            fits.append(fit_model(
                data,
                schedule,
                delta_variant,
                starts=args.starts,
                seed=args.seed + 1000 * schedule_index + delta_index,
            ))

    fit_path = write_fit_table(fits)
    prediction_path = write_prediction_table(fits, data)
    figure_path = plot_four_event_fits(fits, data)
    report_path = write_report(fits)
    print_summary(fits)
    print(
        "WROTE",
        fit_path,
        prediction_path,
        figure_path,
        report_path,
        flush=True,
    )


if __name__ == "__main__":
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.9,
        "savefig.bbox": "tight",
        "svg.fonttype": "none",
        "svg.hashsalt": "sac-cepeda2008-four-event",
    })
    main()

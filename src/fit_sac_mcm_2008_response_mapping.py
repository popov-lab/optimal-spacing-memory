"""Fit deterministic MCM and SAC variants to Cepeda et al. (2008).

Models retained here:
1. Deterministic MCM with one bounded gain delta in (0,1] at every study event
   plus a logistic response mapping.
2. SAC with f(t)=(1+t)^(-d) plus the same logistic response mapping.
3. SAC with f(t)=(1+t/tau)^(-d) plus the same logistic response mapping.

All models are fit jointly to the 26 spacing observations in
``data/cepeda_spacing_recall.csv``. Figures use a linear ISI axis.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "cepeda_spacing_recall.csv"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

N = 100
ZERO_GAP_DAYS = 0.00256


def read_data():
    with DATA.open(newline="") as f:
        rows = list(csv.DictReader(f))
    rows = [r for r in rows if r["panel"] == "d" and r["function"] == "spacing"]
    rows.sort(key=lambda r: (float(r["ri_days"]), float(r["isi_days"])))
    isi = np.array([float(r["isi_days"]) for r in rows])
    ri = np.array([float(r["ri_days"]) for r in rows])
    recall = np.array([float(r["recall_pct"]) / 100.0 for r in rows])
    model_isi = isi.copy()
    model_isi[model_isi == 0] = ZERO_GAP_DAYS
    return rows, isi, model_isi, ri, recall


def logistic_response(strength, theta, sigma):
    return expit((np.asarray(strength) - theta) / sigma)


def mcm_components(mu, nu, xi):
    i = np.arange(1, N + 1, dtype=float)
    log_tau = np.log(mu) + i * np.log(nu)
    limit = np.log(np.finfo(float).max)
    tau = np.empty_like(log_tau)
    finite = log_tau < limit
    tau[finite] = np.exp(log_tau[finite])
    tau[~finite] = np.inf

    log_gamma = i * np.log(xi)
    log_gamma -= log_gamma.max()
    gamma = np.exp(log_gamma)
    gamma /= gamma.sum()
    return tau, gamma


def mcm_strength(isi, ri, mu, nu, xi, delta):
    """Two-study deterministic MCM with the same delta at both study events."""
    tau, gamma = mcm_components(mu, nu, xi)
    isi = np.atleast_1d(np.asarray(isi, dtype=float))
    ri = np.atleast_1d(np.asarray(ri, dtype=float))
    out = np.empty_like(isi)
    cumulative_gamma = np.cumsum(gamma)

    for j, (a, b) in enumerate(zip(isi, ri)):
        # First study: x_i = delta because s_i=0.
        x_pre = delta * np.exp(-a / tau)
        s = np.cumsum(gamma * x_pre) / cumulative_gamma
        x_post = x_pre + delta * (1.0 - s)
        out[j] = float(np.dot(gamma, x_post * np.exp(-b / tau)))
    return out


def sac_forgetting(t, d, tau=None):
    t = np.asarray(t, dtype=float)
    if tau is None:
        return (1.0 + t) ** (-d)
    return (1.0 + t / tau) ** (-d)


def sac_strength(isi, ri, delta, d, tau=None):
    fa = sac_forgetting(isi, d, tau)
    fb = sac_forgetting(ri, d, tau)
    fab = sac_forgetting(isi + ri, d, tau)
    return delta * fb + delta * fab - delta**2 * fa * fb


def fit_mcm(model_isi, ri, observed, starts=96, seed=20260825):
    rng = np.random.default_rng(seed)
    lo = np.array([
        np.log(1e-5), np.log(1e-4), logit(0.005), logit(0.005),
        -2.0, np.log(1e-3),
    ])
    hi = np.array([
        np.log(1e4), np.log(1e2), logit(0.9999), logit(0.9999),
        3.0, np.log(5.0),
    ])

    def unpack(z):
        return {
            "mu": float(np.exp(z[0])),
            "nu": float(1.0 + np.exp(z[1])),
            "xi": float(expit(z[2])),
            "delta": float(expit(z[3])),
            "theta": float(z[4]),
            "sigma": float(np.exp(z[5])),
        }

    def residual(z):
        p = unpack(z)
        strength = mcm_strength(model_isi, ri, p["mu"], p["nu"], p["xi"], p["delta"])
        return logistic_response(strength, p["theta"], p["sigma"]) - observed

    initials = [
        np.array([np.log(10), np.log(.2), logit(.9), logit(.6), .5, np.log(.2)]),
        np.array([np.log(1), np.log(.1), logit(.95), logit(.5), .3, np.log(.15)]),
    ]
    while len(initials) < starts:
        initials.append(lo + rng.random(len(lo)) * (hi - lo))

    best = None
    for z0 in initials:
        result = least_squares(
            residual, z0, bounds=(lo, hi), max_nfev=5000,
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
        )
        sse = float(np.sum(result.fun**2))
        if best is None or sse < best[0]:
            best = (sse, result.x)

    sse, z = best
    p = unpack(z)
    p["rmse_pp"] = 100 * np.sqrt(sse / len(observed))
    return p


def fit_sac(model_isi, ri, observed, with_tau, starts=96, seed=20260826):
    rng = np.random.default_rng(seed)
    lo = [logit(.005), np.log(1e-4)]
    hi = [logit(.9999), np.log(10.0)]
    if with_tau:
        lo += [np.log(1e-4)]
        hi += [np.log(1e5)]
    lo += [-2.0, np.log(1e-3)]
    hi += [3.0, np.log(5.0)]
    lo = np.asarray(lo)
    hi = np.asarray(hi)

    def unpack(z):
        k = 0
        p = {"delta": float(expit(z[k]))}
        k += 1
        p["d"] = float(np.exp(z[k]))
        k += 1
        p["tau"] = float(np.exp(z[k])) if with_tau else None
        if with_tau:
            k += 1
        p["theta"] = float(z[k])
        k += 1
        p["sigma"] = float(np.exp(z[k]))
        return p

    def residual(z):
        p = unpack(z)
        strength = sac_strength(model_isi, ri, p["delta"], p["d"], p["tau"])
        return logistic_response(strength, p["theta"], p["sigma"]) - observed

    initials = []
    for delta0, d0 in [(0.25, .1), (0.5, .3), (0.75, .5), (0.95, 1.0)]:
        z = [logit(delta0), np.log(d0)]
        if with_tau:
            z += [np.log(1.0)]
        z += [.5, np.log(.2)]
        initials.append(np.asarray(z))
    while len(initials) < starts:
        initials.append(lo + rng.random(len(lo)) * (hi - lo))

    best = None
    for z0 in initials:
        result = least_squares(
            residual, z0, bounds=(lo, hi), max_nfev=5000,
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
        )
        sse = float(np.sum(result.fun**2))
        if best is None or sse < best[0]:
            best = (sse, result.x)

    sse, z = best
    p = unpack(z)
    p["rmse_pp"] = 100 * np.sqrt(sse / len(observed))
    return p


def plot_fit(rows, model_name, strength_fn, params, outfile):
    plt.figure(figsize=(7.4, 5.2))
    for RI in [7, 35, 70, 350]:
        condition = [r for r in rows if float(r["ri_days"]) == RI]
        condition.sort(key=lambda r: float(r["isi_days"]))
        xobs = np.array([float(r["isi_days"]) for r in condition])
        yobs = np.array([float(r["recall_pct"]) for r in condition])

        xgrid = np.linspace(0, 105, 500)
        xmodel = xgrid.copy()
        xmodel[xmodel == 0] = ZERO_GAP_DAYS
        rig = np.full_like(xmodel, RI)
        strength = strength_fn(xmodel, rig, params)
        prediction = logistic_response(strength, params["theta"], params["sigma"])

        line, = plt.plot(xgrid, 100 * prediction, label=f"Fit, RI={RI} d")
        plt.plot(xobs, yobs, "o", color=line.get_color(), label=f"Data, RI={RI} d")

    plt.xlim(0, 105)
    plt.ylim(0, 100)
    plt.xticks([0, 7, 14, 21, 35, 70, 105])
    plt.xlabel("ISI (days)")
    plt.ylabel("Final-test recall (%)")
    plt.title(model_name)
    plt.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def main():
    rows, _, model_isi, ri, observed = read_data()
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)

    mcm = fit_mcm(model_isi, ri, observed)
    sac_fixed = fit_sac(model_isi, ri, observed, with_tau=False)
    sac_tau = fit_sac(model_isi, ri, observed, with_tau=True)

    with (RESULTS / "sac_mcm_2008_response_mapping_fits.csv").open("w", newline="") as f:
        fields = [
            "model", "parameters", "delta", "mu", "nu", "xi", "d",
            "tau_days", "theta", "sigma", "rmse_pp"
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow({
            "model": "Deterministic MCM + logistic", "parameters": 6,
            "delta": mcm["delta"], "mu": mcm["mu"], "nu": mcm["nu"],
            "xi": mcm["xi"], "d": "", "tau_days": "",
            "theta": mcm["theta"], "sigma": mcm["sigma"],
            "rmse_pp": mcm["rmse_pp"],
        })
        for label, p, npar in [
            ("SAC fixed scale + logistic", sac_fixed, 4),
            ("SAC with tau + logistic", sac_tau, 5),
        ]:
            writer.writerow({
                "model": label, "parameters": npar, "delta": p["delta"],
                "mu": "", "nu": "", "xi": "", "d": p["d"],
                "tau_days": "" if p["tau"] is None else p["tau"],
                "theta": p["theta"], "sigma": p["sigma"],
                "rmse_pp": p["rmse_pp"],
            })

    plot_fit(
        rows,
        "Cepeda et al. (2008): deterministic MCM + logistic response\n"
        "single bounded gain at all study events",
        lambda a, b, p: mcm_strength(a, b, p["mu"], p["nu"], p["xi"], p["delta"]),
        mcm,
        FIGURES / "mcm_2008_logistic_linear.svg",
    )
    plot_fit(
        rows,
        "Cepeda et al. (2008): SAC + logistic response\n" + r"$f(t)=(1+t)^{-d}$",
        lambda a, b, p: sac_strength(a, b, p["delta"], p["d"], None),
        sac_fixed,
        FIGURES / "sac_2008_logistic_fixed_scale_linear.svg",
    )
    plot_fit(
        rows,
        "Cepeda et al. (2008): SAC + logistic response\n" + r"$f(t)=(1+t/\tau)^{-d}$",
        lambda a, b, p: sac_strength(a, b, p["delta"], p["d"], p["tau"]),
        sac_tau,
        FIGURES / "sac_2008_logistic_tau_linear.svg",
    )


if __name__ == "__main__":
    main()

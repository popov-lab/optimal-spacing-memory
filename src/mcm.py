"""Multiscale Context Model (MCM) from Mozer et al. (2009).

This module implements the leaky-integrator formulation in Eqs. 3, 5--7 and
the stochastic encoding/retrieval marginalization described in Section 4.

Times are expressed in one consistent unit. For the Cepeda replications we use
days, so ``mu`` and all study/test times are in days.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class MCMParams:
    """Primitive MCM parameters.

    Parameters follow the paper's notation: tau_i = mu * nu**i and
    gamma_i proportional to xi**i. ``omega`` is the probability that a study
    episode is encoded. ``epsilon_r`` is the learning rate following successful
    retrieval; following retrieval failure the learning rate is 1.
    """

    mu: float
    nu: float
    omega: float
    xi: float
    epsilon_r: float = 9.0
    n: int = 100

    def __post_init__(self) -> None:
        if self.mu <= 0:
            raise ValueError("mu must be positive")
        if self.nu <= 1:
            raise ValueError("nu must exceed 1")
        if not 0 < self.omega <= 1:
            raise ValueError("omega must lie in (0, 1]")
        if not 0 < self.xi < 1:
            raise ValueError("xi must lie in (0, 1)")
        if self.epsilon_r <= 0:
            raise ValueError("epsilon_r must be positive")
        if self.n < 1:
            raise ValueError("n must be at least 1")


def components(params: MCMParams) -> tuple[np.ndarray, np.ndarray]:
    """Return time constants ``tau`` and normalized weights ``gamma``."""

    i = np.arange(1, params.n + 1, dtype=float)

    # tau can overflow for very large nu at high i. Those components are
    # effectively nondecaying on any finite experimental time scale, so inf is
    # the numerically correct limiting representation.
    log_tau = np.log(params.mu) + i * np.log(params.nu)
    max_log = np.log(np.finfo(float).max)
    tau = np.empty_like(log_tau)
    finite = log_tau < max_log
    tau[finite] = np.exp(log_tau[finite])
    tau[~finite] = np.inf

    # Normalize the geometric weights in log space.
    log_gamma = i * np.log(params.xi)
    log_gamma -= log_gamma.max()
    gamma = np.exp(log_gamma)
    gamma /= gamma.sum()
    return tau, gamma


def decay_factors(dt: float, params: MCMParams) -> np.ndarray:
    """Exponential decay factors over lag ``dt`` for all integrators."""

    if dt < 0:
        raise ValueError("time lags must be nonnegative")
    tau, _ = components(params)
    return np.exp(-dt / tau)


def partial_strengths(x: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    """Return s_i, the gamma-weighted average strength through scale i."""

    return np.cumsum(gamma * x) / np.cumsum(gamma)


def single_study_recall(lags: Iterable[float], params: MCMParams) -> np.ndarray:
    """Expected recall after one study episode as a function of lag.

    Section 4 clarifies that omega is an encoding probability: conditional on
    encoding, all integrators start at 1; otherwise they start at 0. Therefore
    the expectation after one study is omega times the weighted exponential
    mixture.
    """

    lags = np.asarray(list(lags), dtype=float)
    if np.any(lags < 0):
        raise ValueError("lags must be nonnegative")
    tau, gamma = components(params)
    mixture = np.exp(-lags[:, None] / tau[None, :]) @ gamma
    return params.omega * mixture


def expected_recall(
    study_times: Iterable[float],
    test_time: float,
    params: MCMParams,
) -> float:
    """Exact expected recall for a sequence of study episodes.

    The function explicitly marginalizes over retrieval success/failure and
    encoding success/failure at each study episode, exactly as described in
    Section 4 of Mozer et al. (2009). With ``s`` study episodes there are at
    most 2**(2*s - 1) nonzero outcome histories because recall necessarily
    fails before the first study.

    Each element of ``study_times`` is treated as one MCM study episode. This
    matches the paper's abstraction of the Cepeda experiments at the session
    level; within-session passes are not modeled as separate episodes.
    """

    times = np.asarray(list(study_times), dtype=float)
    if times.ndim != 1 or len(times) == 0:
        raise ValueError("study_times must be a nonempty one-dimensional sequence")
    if np.any(np.diff(times) < 0):
        raise ValueError("study_times must be nondecreasing")
    if test_time < times[-1]:
        raise ValueError("test_time must not precede the last study")

    tau, gamma = components(params)
    zero = np.zeros(params.n, dtype=float)
    states: list[tuple[float, np.ndarray]] = [(1.0, zero)]
    previous_time = times[0]

    for episode, study_time in enumerate(times):
        dt = 0.0 if episode == 0 else study_time - previous_time
        if dt:
            d = np.exp(-dt / tau)
            states = [(prob, x * d) for prob, x in states]

        new_states: list[tuple[float, np.ndarray]] = []
        for history_prob, x in states:
            s = partial_strengths(x, gamma)
            p_recall = float(np.clip(s[-1], 0.0, 1.0))

            # Retrieval is assessed before study and selects the learning rate.
            for recall_prob, epsilon in (
                (1.0 - p_recall, 1.0),
                (p_recall, params.epsilon_r),
            ):
                if recall_prob <= 0:
                    continue

                branch_prob = history_prob * recall_prob

                # Encoding failure: study has no effect.
                p_no_encode = branch_prob * (1.0 - params.omega)
                if p_no_encode > 0:
                    new_states.append((p_no_encode, x.copy()))

                # Encoding success: Eq. 7, Delta x_i = epsilon (1 - s_i).
                p_encode = branch_prob * params.omega
                if p_encode > 0:
                    updated = x + epsilon * (1.0 - s)
                    new_states.append((p_encode, updated))

        states = new_states
        previous_time = study_time

    final_decay = np.exp(-(test_time - times[-1]) / tau)
    expected = 0.0
    for history_prob, x in states:
        strength = float(np.dot(gamma, x * final_decay))
        # Paper specifies min(1, s_N). Numerical lower clipping only protects
        # probability semantics if a pathological parameter/history combination
        # drives strength below zero after an overshooting update.
        p_recall = float(np.clip(strength, 0.0, 1.0))
        expected += history_prob * p_recall

    return expected


def two_session_spacing(
    isi: float,
    ri: float,
    params: MCMParams,
) -> float:
    """Convenience wrapper for study at 0 and ``isi``, tested after ``ri``."""

    if isi < 0 or ri < 0:
        raise ValueError("isi and ri must be nonnegative")
    return expected_recall([0.0, isi], isi + ri, params)

"""Models compared in Anderson et al. (2023), *Psychological Review*.

Environmental ``ages`` are positive times from each occurrence to prediction;
an occurrence in the immediately preceding text has age 1. The authors'
behavioral AMPE schedules additionally use age 0, which the recall function
accepts explicitly. Functions whose paper equation produces odds are named
``*_odds``; use :func:`odds_to_probability` for environmental probabilities.

The implementation follows the equations rather than the paper's data-binning
and fitting code.  Documented corrections and code/paper discrepancies are in
``docs/anderson_2023_model_reference.md``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np


ArrayLike = Iterable[float] | np.ndarray
DecayKind = Literal["exponential", "power"]
OccurrenceMapping = Literal["odds", "released_probability"]
OutputScaleSemantics = Literal["odds", "released_probability"]


def _ages_oldest_to_newest(ages: ArrayLike) -> np.ndarray:
    """Validate ages and return them in chronological occurrence order."""

    values = np.asarray(list(ages), dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("ages must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("all ages must be finite and strictly positive")
    return np.sort(values)[::-1]


def _behavioral_ages_oldest_to_newest(ages: ArrayLike) -> np.ndarray:
    """Validate behavioral AMPE ages, for which the release permits zero."""

    values = np.asarray(list(ages), dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("ages must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("behavioral ages must be finite and nonnegative")
    return np.sort(values)[::-1]


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return value


def _nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def odds_to_probability(odds: float | np.ndarray) -> float | np.ndarray:
    """Convert nonnegative odds to probability without overflowing."""

    values = np.asarray(odds, dtype=float)
    if np.any(values < 0) or np.any(np.isnan(values)):
        raise ValueError("odds must be nonnegative and not NaN")
    probabilities = np.divide(
        values,
        1.0 + values,
        out=np.ones_like(values),
        where=~np.isinf(values),
    )
    return float(probabilities) if probabilities.ndim == 0 else probabilities


def logistic(value: float | np.ndarray) -> float | np.ndarray:
    """Numerically stable logistic function."""

    values = np.asarray(value, dtype=float)
    positive = values >= 0
    out = np.empty_like(values)
    out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    exp_values = np.exp(values[~positive])
    out[~positive] = exp_values / (1.0 + exp_values)
    return float(out) if out.ndim == 0 else out


def harmonic_mean(values: ArrayLike) -> float:
    """Harmonic mean of finite, strictly positive values."""

    values = np.asarray(list(values), dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("values must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("harmonic-mean inputs must be finite and positive")
    return float(values.size / np.sum(1.0 / values))


def inclusive_range(ages: ArrayLike) -> float:
    """Inclusive first-to-last span used by the released AMPE code.

    A singleton has range 1 and occurrences in adjacent texts have range 2.
    """

    values = _ages_oldest_to_newest(ages)
    return float(values[0] - values[-1] + 1.0)


@dataclass(frozen=True)
class GPEParameters:
    odds_scale: float
    frequency_exponent: float
    decay: float


def gpe_odds(ages: ArrayLike, params: GPEParameters) -> float:
    """General Performance Equation (paper Equation 3)."""

    values = _ages_oldest_to_newest(ages)
    scale = _positive(params.odds_scale, "odds_scale")
    decay = _nonnegative(params.decay, "decay")
    frequency_exponent = _finite(params.frequency_exponent, "frequency_exponent")
    return float(scale * values.size**frequency_exponent * values[-1] ** (-decay))


@dataclass(frozen=True)
class ACTRParameters:
    odds_scale: float
    decay: float


def actr_odds(ages: ArrayLike, params: ACTRParameters) -> float:
    """ACT-R base-level odds (paper Equation 4)."""

    values = _ages_oldest_to_newest(ages)
    scale = _positive(params.odds_scale, "odds_scale")
    decay = _nonnegative(params.decay, "decay")
    return float(scale * np.sum(values ** (-decay)))


@dataclass(frozen=True)
class PavlikAndersonParameters:
    odds_scale: float
    minimum_decay: float
    activation_sensitivity: float


def pavlik_anderson_component_decays(
    ages: ArrayLike, params: PavlikAndersonParameters
) -> np.ndarray:
    """Assign one decay exponent per occurrence (paper Equations A1-A2)."""

    values = _ages_oldest_to_newest(ages)
    if np.any(values[:-1] <= values[1:]):
        raise ValueError("P&A requires distinct occurrence times")
    minimum = _nonnegative(params.minimum_decay, "minimum_decay")
    sensitivity = _nonnegative(
        params.activation_sensitivity, "activation_sensitivity"
    )

    decays = np.empty(values.size, dtype=float)
    decays[0] = minimum
    for occurrence in range(1, values.size):
        elapsed = values[:occurrence] - values[occurrence]
        prior_strength = np.sum(elapsed ** (-decays[:occurrence]))
        decays[occurrence] = minimum + sensitivity * prior_strength
    return decays


def pavlik_anderson_odds(
    ages: ArrayLike, params: PavlikAndersonParameters
) -> float:
    """Pavlik-Anderson odds (paper Equation 5 with Equations A1-A2)."""

    values = _ages_oldest_to_newest(ages)
    scale = _positive(params.odds_scale, "odds_scale")
    decays = pavlik_anderson_component_decays(values, params)
    return float(scale * np.sum(values ** (-decays)))


@dataclass(frozen=True)
class PPEParameters:
    odds_scale: float
    frequency_exponent: float
    recency_weight: float
    minimum_decay: float
    spacing_sensitivity: float


@dataclass(frozen=True)
class PPEComponents:
    effective_age: float
    decay: float


def ppe_components(ages: ArrayLike, params: PPEParameters) -> PPEComponents:
    """Return PPE effective age and spacing-dependent decay (A3-A4)."""

    values = _ages_oldest_to_newest(ages)
    recency_weight = float(params.recency_weight)
    if not np.isfinite(recency_weight):
        raise ValueError("recency_weight must be finite")
    minimum = _nonnegative(params.minimum_decay, "minimum_decay")
    sensitivity = _nonnegative(params.spacing_sensitivity, "spacing_sensitivity")

    log_weights = -recency_weight * np.log(values)
    weights = np.exp(log_weights - np.max(log_weights))
    weights /= np.sum(weights)
    effective_age = float(np.sum(weights * values))

    if values.size == 1:
        decay = minimum
    else:
        lags = values[:-1] - values[1:]
        if np.any(lags < 0):
            raise ValueError("occurrence ages imply negative lags")
        spacing_term = np.mean(1.0 / np.log(lags + np.e))
        decay = float(minimum + sensitivity * spacing_term)
    return PPEComponents(effective_age=effective_age, decay=decay)


def ppe_odds(ages: ArrayLike, params: PPEParameters) -> float:
    """Predictive Performance Equation odds (paper Equation 6)."""

    values = _ages_oldest_to_newest(ages)
    scale = _positive(params.odds_scale, "odds_scale")
    frequency_exponent = _finite(params.frequency_exponent, "frequency_exponent")
    components = ppe_components(values, params)
    return float(
        scale
        * values.size**frequency_exponent
        * components.effective_age ** (-components.decay)
    )


@dataclass(frozen=True)
class MCMParameters:
    odds_scale: float
    time_scale: float
    time_ratio: float
    total_weight: float
    weight_ratio: float
    n_traces: int = 100


@dataclass(frozen=True)
class MCMState:
    trace_strengths: np.ndarray
    time_constants: np.ndarray
    trace_weights: np.ndarray
    weighted_strength: float


def mcm_state(ages: ArrayLike, params: MCMParameters) -> MCMState:
    """Run the released environmental MCM state update through prediction.

    The released MATLAB code computes all presentation increments from the
    pre-update state and truncates negative increments at zero.  Both details
    were omitted from the printed equations and are made explicit here.
    """

    values = _ages_oldest_to_newest(ages)
    time_scale = _positive(params.time_scale, "time_scale")
    time_ratio = _positive(params.time_ratio, "time_ratio")
    total_weight = _positive(params.total_weight, "total_weight")
    weight_ratio = _positive(params.weight_ratio, "weight_ratio")
    if time_ratio <= 1:
        raise ValueError("time_ratio must exceed 1 for successively slower traces")
    if total_weight >= 1:
        raise ValueError("total_weight must be below 1")
    if weight_ratio >= 1:
        raise ValueError("weight_ratio must be below 1 for decreasing weights")
    if int(params.n_traces) != params.n_traces or params.n_traces <= 0:
        raise ValueError("n_traces must be a positive integer")

    indices = np.arange(1, int(params.n_traces) + 1, dtype=float)
    time_constants = time_scale * time_ratio**indices
    unscaled_weights = weight_ratio**indices
    trace_weights = total_weight * unscaled_weights / np.sum(unscaled_weights)
    cumulative_weights = np.cumsum(trace_weights)

    strengths = np.ones_like(indices)
    presentation_lags = values[:-1] - values[1:]
    for lag in presentation_lags:
        strengths *= np.exp(-lag / time_constants)
        cumulative_strength = np.cumsum(trace_weights * strengths)
        increments = np.maximum(0.0, 1.0 - cumulative_strength / cumulative_weights)
        strengths += increments

    strengths *= np.exp(-values[-1] / time_constants)
    weighted_strength = float(np.sum(trace_weights * strengths))
    return MCMState(
        trace_strengths=strengths,
        time_constants=time_constants,
        trace_weights=trace_weights,
        weighted_strength=weighted_strength,
    )


def mcm_odds(ages: ArrayLike, params: MCMParameters) -> float:
    """Environmental MCM odds (paper Equation 7)."""

    scale = _positive(params.odds_scale, "odds_scale")
    strength = min(mcm_state(ages, params).weighted_strength, 0.999999)
    return float(scale * strength / (1.0 - strength))


@dataclass(frozen=True)
class AMPEParameters:
    desirability_scale: float
    decay_scale: float
    prior_age: float
    prior_range: float


@dataclass(frozen=True)
class AMPEComponents:
    frequency: int
    range: float
    currency: float
    effective_interval: float
    decay: float
    desirability: float
    log_odds: float


def ampe_components(
    ages: ArrayLike,
    params: AMPEParameters,
    *,
    range_value: float | None = None,
) -> AMPEComponents:
    """Return all AMPE quantities in paper Equations 10-13."""

    values = _ages_oldest_to_newest(ages)
    if np.any(values[:-1] <= values[1:]):
        raise ValueError("environmental AMPE requires distinct occurrence times")
    desirability_scale = _positive(params.desirability_scale, "desirability_scale")
    decay_scale = _positive(params.decay_scale, "decay_scale")
    prior_age = _positive(params.prior_age, "prior_age")
    prior_range = _positive(params.prior_range, "prior_range")
    observed_range = (
        inclusive_range(values)
        if range_value is None
        else _positive(range_value, "range_value")
    )

    currency = harmonic_mean(np.append(values, prior_age)) + 1.0
    effective_interval = (observed_range + prior_range) / 2.0
    decay = decay_scale / effective_interval
    desirability = desirability_scale * values.size / effective_interval
    log_odds = np.log(desirability) - decay * np.log(currency)
    return AMPEComponents(
        frequency=int(values.size),
        range=float(observed_range),
        currency=float(currency),
        effective_interval=float(effective_interval),
        decay=float(decay),
        desirability=float(desirability),
        log_odds=float(log_odds),
    )


def ampe_odds(
    ages: ArrayLike,
    params: AMPEParameters,
    *,
    range_value: float | None = None,
) -> float:
    """AMPE environmental odds (paper Equations 9-13)."""

    log_odds = ampe_components(ages, params, range_value=range_value).log_odds
    return float(np.exp(np.clip(log_odds, -745.0, 709.0)))


@dataclass(frozen=True)
class AMPERecallParameters:
    decay_scale: float
    prior_age: float
    prior_range: float
    threshold: float
    noise_scale: float


def ampe_recall_probability(
    ages: ArrayLike,
    params: AMPERecallParameters,
    *,
    range_value: float | None = None,
) -> float:
    """Human-recall AMPE (paper Equation 15).

    ``threshold`` is the paper's eta: the retrieval threshold after absorbing
    the environmental desirability scale.
    """

    values = _behavioral_ages_oldest_to_newest(ages)
    decay_scale = _positive(params.decay_scale, "decay_scale")
    prior_age = _positive(params.prior_age, "prior_age")
    prior_range = _positive(params.prior_range, "prior_range")
    noise_scale = _positive(params.noise_scale, "noise_scale")
    threshold = _finite(params.threshold, "threshold")
    observed_range = (
        float(np.max(values) - np.min(values) + 1.0)
        if range_value is None
        else _positive(range_value, "range_value")
    )

    # MATLAB's harmmean returns zero if any released behavioral age is zero.
    # The paper's +1 then makes currency exactly one for such a schedule.
    augmented_ages = np.append(values, prior_age)
    currency = (
        1.0
        if np.any(augmented_ages == 0)
        else harmonic_mean(augmented_ages) + 1.0
    )
    effective_interval = (observed_range + prior_range) / 2.0
    decay = decay_scale / effective_interval
    activation = (
        np.log(values.size)
        - np.log(effective_interval)
        - decay * np.log(currency)
    )
    return float(logistic((activation - threshold) / noise_scale))


def odds_to_reaction_time(
    odds: float | np.ndarray,
    intercept: float,
    scale: float,
    exponent: float,
) -> float | np.ndarray:
    """Paper Equation 16: ``RT = intercept + scale * odds**(-exponent)``."""

    values = np.asarray(odds, dtype=float)
    if np.any(values <= 0) or np.any(~np.isfinite(values)):
        raise ValueError("reaction-time odds must be finite and positive")
    intercept = _finite(intercept, "intercept")
    scale = _nonnegative(scale, "scale")
    exponent = _nonnegative(exponent, "exponent")
    result = intercept + scale * values ** (-exponent)
    return float(result) if result.ndim == 0 else result


@dataclass(frozen=True)
class AndersonMilsonParameters:
    desirability_shape: float
    desirability_scale: float
    mean_decay: float
    mean_revival_interval: float


@dataclass(frozen=True)
class AndersonMilsonSimulation:
    occurrences: np.ndarray
    probabilities: np.ndarray
    prediction_values: np.ndarray
    odds: np.ndarray
    elapsed_since_revival: np.ndarray
    revivals: np.ndarray
    initial_desirability: np.ndarray
    decay: np.ndarray
    occurrence_mapping: OccurrenceMapping


@dataclass(frozen=True)
class AndersonMilsonConditionalPredictions:
    """Monte Carlo estimates conditional on the paper's history summary.

    Each row is keyed by frequency, most-recent age, and second-most-recent
    age.  The last value is zero for a singleton history.
    """

    history_summaries: np.ndarray
    mean_prediction_values: np.ndarray
    counts: np.ndarray


def scale_anderson_milson_predictions(
    conditional_predictions: float | np.ndarray,
    output_scale: float,
    *,
    semantics: OutputScaleSemantics = "odds",
) -> float | np.ndarray:
    """Apply A&M's output scale *after* conditional prediction.

    ``semantics="odds"`` implements the paper's statement that the scale
    multiplies odds. ``semantics="released_probability"`` reproduces the
    MATLAB shortcut, which multiplies already-averaged raw scores and can
    therefore return a value above one.
    """

    values = np.asarray(conditional_predictions, dtype=float)
    scale = _positive(output_scale, "output_scale")
    if semantics == "odds":
        if np.any(~np.isfinite(values)) or np.any((values < 0) | (values > 1)):
            raise ValueError("odds scaling requires values within [0, 1]")
        odds = np.divide(
            values,
            1.0 - values,
            out=np.full_like(values, np.inf),
            where=values < 1.0,
        )
        result = odds_to_probability(scale * odds)
    elif semantics == "released_probability":
        if np.any(~np.isfinite(values)) or np.any(values < 0):
            raise ValueError(
                "released prediction values must be finite and nonnegative"
            )
        result = scale * values
    else:
        raise ValueError("semantics must be 'odds' or 'released_probability'")
    return float(result) if np.ndim(result) == 0 else result


def released_geometric_mean_scale(
    observed_probabilities: ArrayLike,
    predicted_values: ArrayLike,
) -> float:
    """Recover the release's post-hoc probability multiplier.

    This equals ``exp(mean(log(observed)) - mean(log(predicted)))`` and is not
    the arithmetic-mean or odds calibration described in the paper.
    """

    observed = np.asarray(list(observed_probabilities), dtype=float)
    predicted = np.asarray(list(predicted_values), dtype=float)
    if (
        observed.ndim != 1
        or predicted.ndim != 1
        or observed.shape != predicted.shape
    ):
        raise ValueError("observed and predicted must be same-length vectors")
    if observed.size == 0:
        raise ValueError("probability vectors must not be empty")
    if (
        np.any(~np.isfinite(observed))
        or np.any(~np.isfinite(predicted))
        or np.any(observed <= 0)
        or np.any(predicted <= 0)
        or np.any(observed > 1)
    ):
        raise ValueError(
            "geometric-mean calibration requires observed probabilities in "
            "(0, 1] and finite positive predicted values"
        )
    return float(np.exp(np.mean(np.log(observed)) - np.mean(np.log(predicted))))


def simulate_anderson_milson(
    n_items: int,
    n_steps: int,
    params: AndersonMilsonParameters,
    *,
    decay_kind: DecayKind,
    occurrence_mapping: OccurrenceMapping = "odds",
    rng: int | np.random.Generator | None = None,
) -> AndersonMilsonSimulation:
    """Simulate the A&M latent environment with explicit discrete timing.

    Desirability is sampled as Gamma(shape, scale), decay as Exponential with
    the supplied mean, and revivals as a Poisson process with mean interval
    ``mean_revival_interval`` observed at unit event boundaries. A sampled text
    immediately following a revival has age one, so the printed power law
    ``age**(-d)`` is finite and equals one. The exponential law is
    ``exp(-d*age)`` on the same age convention.

    By default latent desirability is treated as odds and converted to a
    Bernoulli probability. Set ``occurrence_mapping="released_probability"``
    to reproduce the MATLAB shortcut: Bernoulli draws use the value directly
    (implicitly capped at one), while conditional predictions average its raw,
    potentially above-one value. The fitted output scale is deliberately
    absent: it belongs after histories have been conditioned and averaged, via
    :func:`scale_anderson_milson_predictions`.
    """

    if int(n_items) != n_items or n_items <= 0:
        raise ValueError("n_items must be a positive integer")
    if int(n_steps) != n_steps or n_steps <= 0:
        raise ValueError("n_steps must be a positive integer")
    if decay_kind not in ("exponential", "power"):
        raise ValueError("decay_kind must be 'exponential' or 'power'")
    if occurrence_mapping not in ("odds", "released_probability"):
        raise ValueError(
            "occurrence_mapping must be 'odds' or 'released_probability'"
        )

    shape = _positive(params.desirability_shape, "desirability_shape")
    desirability_scale = _positive(
        params.desirability_scale, "desirability_scale"
    )
    mean_decay = _positive(params.mean_decay, "mean_decay")
    mean_revival = _positive(
        params.mean_revival_interval, "mean_revival_interval"
    )
    generator = (
        rng if isinstance(rng, np.random.Generator) else np.random.default_rng(rng)
    )

    initial_desirability = generator.gamma(shape, desirability_scale, int(n_items))
    decays = generator.exponential(mean_decay, int(n_items))
    revival_probability = -np.expm1(-1.0 / mean_revival)
    elapsed = generator.geometric(revival_probability, int(n_items)).astype(float)

    odds = np.empty((int(n_steps), int(n_items)), dtype=float)
    probabilities = np.empty_like(odds)
    prediction_values = np.empty_like(odds)
    occurrences = np.empty_like(odds, dtype=bool)
    elapsed_history = np.empty_like(odds)
    revival_history = np.empty_like(odds, dtype=bool)

    for step in range(int(n_steps)):
        if step == 0:
            revived = elapsed == 1.0
        else:
            revived = generator.random(int(n_items)) < revival_probability
            elapsed[revived] = 1.0
            elapsed[~revived] += 1.0
        if decay_kind == "exponential":
            retention = np.exp(-decays * elapsed)
        else:
            retention = elapsed ** (-decays)
        odds[step] = initial_desirability * retention
        if occurrence_mapping == "odds":
            probabilities[step] = odds_to_probability(odds[step])
            prediction_values[step] = probabilities[step]
        else:
            probabilities[step] = np.minimum(odds[step], 1.0)
            # MATLAB averages the raw value, even when it exceeds one; only
            # the Bernoulli comparison clips implicitly.
            prediction_values[step] = odds[step]
        occurrences[step] = generator.random(int(n_items)) < probabilities[step]
        elapsed_history[step] = elapsed
        revival_history[step] = revived

    return AndersonMilsonSimulation(
        occurrences=occurrences,
        probabilities=probabilities,
        prediction_values=prediction_values,
        odds=odds,
        elapsed_since_revival=elapsed_history,
        revivals=revival_history,
        initial_desirability=initial_desirability,
        decay=decays,
        occurrence_mapping=occurrence_mapping,
    )


def anderson_milson_conditional_predictions(
    simulation: AndersonMilsonSimulation,
    *,
    history_length: int = 1000,
    max_frequency: int | None = 225,
) -> AndersonMilsonConditionalPredictions:
    """Estimate next-event probability conditional on a history summary.

    This is the operational prediction stage omitted by a bare A&M generator.
    It averages each target's prediction value within cells keyed by
    ``(N, most_recent_age, second_most_recent_age)``. This is a probability in
    the coherent odds mode and the raw, potentially above-one MATLAB score in
    released mode. Histories with no occurrence, or above ``max_frequency``
    when supplied, are omitted. The function deliberately does not reproduce
    the paper's subsequent binning.
    """

    if int(history_length) != history_length or history_length <= 0:
        raise ValueError("history_length must be a positive integer")
    if max_frequency is not None and (
        int(max_frequency) != max_frequency or max_frequency <= 0
    ):
        raise ValueError("max_frequency must be a positive integer or None")

    occurrences = np.asarray(simulation.occurrences, dtype=bool)
    prediction_values = np.asarray(simulation.prediction_values, dtype=float)
    if occurrences.ndim != 2 or prediction_values.shape != occurrences.shape:
        raise ValueError("simulation occurrence and prediction arrays must align")
    n_steps, n_items = occurrences.shape
    if history_length >= n_steps:
        raise ValueError("history_length must be smaller than the simulation")

    probability_sums: defaultdict[tuple[int, int, int], float] = defaultdict(float)
    counts: defaultdict[tuple[int, int, int], int] = defaultdict(int)
    for item in range(n_items):
        positions = np.flatnonzero(occurrences[:, item])
        for target in range(int(history_length), n_steps):
            left = np.searchsorted(positions, target - history_length, side="left")
            right = np.searchsorted(positions, target, side="left")
            history_positions = positions[left:right]
            frequency = int(history_positions.size)
            if frequency == 0 or (
                max_frequency is not None and frequency > max_frequency
            ):
                continue
            most_recent_age = int(target - history_positions[-1])
            second_most_recent_age = (
                int(target - history_positions[-2]) if frequency > 1 else 0
            )
            key = (frequency, most_recent_age, second_most_recent_age)
            probability_sums[key] += float(prediction_values[target, item])
            counts[key] += 1

    keys = sorted(counts)
    return AndersonMilsonConditionalPredictions(
        history_summaries=np.asarray(keys, dtype=int).reshape(-1, 3),
        mean_prediction_values=np.asarray(
            [probability_sums[key] / counts[key] for key in keys], dtype=float
        ),
        counts=np.asarray([counts[key] for key in keys], dtype=int),
    )

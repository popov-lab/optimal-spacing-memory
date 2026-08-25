"""Very small qualitative checks for the Anderson et al. (2023) models.

This is intentionally not a reproduction or model-fitting script. It checks
only that each implementation returns finite probabilities, declines when the
same history is moved farther into the past, that AMPE shows its intended
spacing crossover, and that short A&M simulations produce valid probabilities.
Released parameter values are used only as convenient nominal scales; where
the corrected equations differ from the MATLAB evaluator they are not refits.

Usage: python scripts/sanity_check_anderson_2023.py
"""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from anderson_2023_models import (  # noqa: E402
    ACTRParameters,
    AMPEParameters,
    AndersonMilsonParameters,
    GPEParameters,
    MCMParameters,
    PPEParameters,
    PavlikAndersonParameters,
    actr_odds,
    ampe_odds,
    gpe_odds,
    mcm_odds,
    odds_to_probability,
    pavlik_anderson_odds,
    ppe_odds,
    simulate_anderson_milson,
)


MODELS = {
    "GPE": lambda ages: gpe_odds(
        ages,
        GPEParameters(
            alpha=0.02081011,
            frequency_exponent=0.58161261,
            decay=0.61681564,
        ),
    ),
    "ACT-R": lambda ages: actr_odds(
        ages, ACTRParameters(alpha=0.04010094, decay=0.79718902)
    ),
    "P&A": lambda ages: pavlik_anderson_odds(
        ages,
        PavlikAndersonParameters(
            alpha=0.05296011,
            minimum_decay=0.75811324,
            activation_sensitivity=0.45391362,
        ),
    ),
    "PPE": lambda ages: ppe_odds(
        ages,
        PPEParameters(
            alpha=0.0180,
            frequency_exponent=0.6178,
            recency_weight=8.6986,
            minimum_decay=0.5358,
            spacing_sensitivity=0.1862,
        ),
    ),
    "MCM": lambda ages: mcm_odds(
        ages,
        MCMParameters(
            alpha=0.0288,
            time_scale=0.0316,
            time_ratio=1.1112,
            total_weight=0.7041,
            weight_ratio=0.9784,
        ),
    ),
    "AMPE": lambda ages: ampe_odds(
        ages,
        AMPEParameters(
            alpha=214.1079,
            decay_scale=1401.1350,
            prior_age=15.1745,
            prior_range=1564.9791,
        ),
    ),
}


def check_closed_form_models() -> None:
    recent = [40, 5]
    remote = [140, 105]
    print("Closed-form models: same two-occurrence history shifted 100 events")
    print("model       p(recent)    p(remote)")
    for name, predictor in MODELS.items():
        p_recent = odds_to_probability(predictor(recent))
        p_remote = odds_to_probability(predictor(remote))
        assert 0 < p_remote < p_recent < 1
        print(f"{name:<8} {p_recent:>11.6f} {p_remote:>12.6f}")


def check_ampe_crossover() -> None:
    params = AMPEParameters(
        alpha=214.1079,
        decay_scale=1401.1350,
        prior_age=15.1745,
        prior_range=1564.9791,
    )

    def probability(gap: int, delay: int) -> float:
        return odds_to_probability(ampe_odds([gap + delay, delay], params))

    short_massed = probability(gap=1, delay=1)
    short_spaced = probability(gap=500, delay=1)
    long_massed = probability(gap=1, delay=500)
    long_spaced = probability(gap=500, delay=500)
    assert short_massed > short_spaced
    assert long_spaced > long_massed
    print("\nAMPE spacing crossover")
    print(f"delay 1:   gap 1 = {short_massed:.6f}, gap 500 = {short_spaced:.6f}")
    print(f"delay 500: gap 1 = {long_massed:.6f}, gap 500 = {long_spaced:.6f}")


def check_environment_simulations() -> None:
    fits = {
        "exponential": AndersonMilsonParameters(
            desirability_shape=0.164,
            gamma_scale=0.139,
            mean_decay=0.035,
            mean_revival_interval=333,
        ),
        "power": AndersonMilsonParameters(
            desirability_shape=0.199,
            gamma_scale=0.482,
            mean_decay=4.076,
            mean_revival_interval=800,
        ),
    }
    print("\nA&M simulations: 500 items x 200 events")
    print("decay           mean p   occurrence rate   p(at revival)")
    for kind, params in fits.items():
        simulation = simulate_anderson_milson(
            500, 200, params, decay_kind=kind, rng=20230820
        )
        assert simulation.probabilities.shape == (200, 500)
        assert np.all((simulation.probabilities >= 0) & (simulation.probabilities <= 1))
        at_revival = simulation.revivals
        revival_probability = float(np.mean(simulation.probabilities[at_revival]))
        mean_probability = float(np.mean(simulation.probabilities))
        occurrence_rate = float(np.mean(simulation.occurrences))
        assert revival_probability > mean_probability
        print(
            f"{kind:<12} {mean_probability:>10.6f}"
            f" {occurrence_rate:>17.6f} {revival_probability:>15.6f}"
        )


if __name__ == "__main__":
    check_closed_form_models()
    check_ampe_crossover()
    check_environment_simulations()
    print("\nAll qualitative checks passed.")

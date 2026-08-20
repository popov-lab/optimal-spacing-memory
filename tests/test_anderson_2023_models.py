"""Unit tests for the Anderson et al. (2023) model implementations."""

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from anderson_2023_models import (  # noqa: E402
    ACTRParameters,
    AMPEParameters,
    AMPERecallParameters,
    AndersonMilsonParameters,
    AndersonMilsonSimulation,
    GPEParameters,
    MCMParameters,
    PPEParameters,
    PavlikAndersonParameters,
    actr_odds,
    ampe_components,
    ampe_odds,
    ampe_recall_probability,
    anderson_milson_conditional_predictions,
    gpe_odds,
    inclusive_range,
    mcm_odds,
    mcm_state,
    odds_to_probability,
    odds_to_reaction_time,
    pavlik_anderson_component_decays,
    ppe_components,
    released_geometric_mean_scale,
    scale_anderson_milson_predictions,
    simulate_anderson_milson,
)


class TestHelpers(unittest.TestCase):
    def test_odds_conversion(self) -> None:
        self.assertEqual(odds_to_probability(0), 0)
        self.assertAlmostEqual(odds_to_probability(3), 0.75)
        self.assertEqual(odds_to_probability(np.inf), 1)

    def test_inclusive_range(self) -> None:
        self.assertEqual(inclusive_range([7]), 1)
        self.assertEqual(inclusive_range([2, 1]), 2)
        self.assertEqual(inclusive_range([13, 5, 9]), 9)


class TestClosedFormModels(unittest.TestCase):
    def test_gpe_and_actr_exact_values(self) -> None:
        gpe = GPEParameters(
            odds_scale=0.5, frequency_exponent=1.0, decay=2.0
        )
        actr = ACTRParameters(odds_scale=0.1, decay=1.0)
        self.assertAlmostEqual(gpe_odds([10, 2], gpe), 0.25)
        self.assertAlmostEqual(actr_odds([10, 2], actr), 0.06)

    def test_pavlik_anderson_decay_recursion(self) -> None:
        params = PavlikAndersonParameters(
            odds_scale=0.1,
            minimum_decay=0.5,
            activation_sensitivity=0.2,
        )
        decays = pavlik_anderson_component_decays([10, 6, 2], params)
        self.assertAlmostEqual(decays[0], 0.5)
        self.assertAlmostEqual(decays[1], 0.5 + 0.2 * 4**-0.5)
        expected_third = 0.5 + 0.2 * (
            8**-decays[0] + 4**-decays[1]
        )
        self.assertAlmostEqual(decays[2], expected_third)

    def test_ppe_singleton_and_weighted_age(self) -> None:
        params = PPEParameters(
            odds_scale=0.02,
            frequency_exponent=0.6,
            recency_weight=1.0,
            minimum_decay=0.5,
            spacing_sensitivity=0.2,
        )
        singleton = ppe_components([8], params)
        self.assertEqual(singleton.effective_age, 8)
        self.assertEqual(singleton.decay, 0.5)
        repeated = ppe_components([4, 2], params)
        self.assertAlmostEqual(repeated.effective_age, 8 / 3)
        self.assertAlmostEqual(repeated.decay, 0.5 + 0.2 / np.log(2 + np.e))

    def test_mcm_state_and_probability(self) -> None:
        params = MCMParameters(
            odds_scale=0.1,
            time_scale=1.0,
            time_ratio=2.0,
            total_weight=0.5,
            weight_ratio=0.5,
            n_traces=2,
        )
        state = mcm_state([30, 5], params)
        self.assertAlmostEqual(np.sum(state.trace_weights), 0.5)
        self.assertGreater(state.weighted_strength, 0)
        probability = odds_to_probability(mcm_odds([30, 5], params))
        self.assertGreater(probability, 0)
        self.assertLess(probability, 1)

        # Independent two-trace benchmark for a lag of 2 and test age of 1.
        first_decay = np.array([np.exp(-1), np.exp(-0.5)])
        weights = np.array([1 / 3, 1 / 6])
        means = np.cumsum(weights * first_decay) / np.cumsum(weights)
        updated = first_decay + np.maximum(0, 1 - means)
        expected = np.sum(weights * updated * np.exp(-1 / np.array([2, 4])))
        self.assertAlmostEqual(mcm_state([3, 1], params).weighted_strength, expected)

    def test_ampe_components_and_recall(self) -> None:
        environmental = AMPEParameters(
            desirability_scale=214,
            decay_scale=1401,
            prior_age=15.18,
            prior_range=1565,
        )
        components = ampe_components([101, 1], environmental)
        self.assertEqual(components.range, 101)
        self.assertEqual(components.frequency, 2)
        self.assertGreater(ampe_odds([101, 1], environmental), 0)

        recall = AMPERecallParameters(
            decay_scale=82.3,
            prior_age=2,
            prior_range=103,
            threshold=-6.51,
            noise_scale=0.63,
        )
        probability = ampe_recall_probability([20, 1], recall)
        self.assertGreater(probability, 0)
        self.assertLess(probability, 1)

        # Exact Appendix.m formula check for a released Bahrick age-0 pattern.
        zero_age_params = AMPERecallParameters(
            decay_scale=500,
            prior_age=500,
            prior_range=500,
            threshold=-8.01,
            noise_scale=0.38,
        )
        zero_age_probability = ampe_recall_probability(
            [0, 0, 0], zero_age_params, range_value=1
        )
        self.assertAlmostEqual(zero_age_probability, 0.999920096750237)

    def test_ampe_spacing_crossover_with_released_fit(self) -> None:
        params = AMPEParameters(
            desirability_scale=214.1079,
            decay_scale=1401.1350,
            prior_age=15.1745,
            prior_range=1564.9791,
        )

        def probability(gap: int, delay: int) -> float:
            return odds_to_probability(ampe_odds([gap + delay, delay], params))

        self.assertGreater(probability(1, 1), probability(500, 1))
        self.assertGreater(probability(500, 500), probability(1, 500))

    def test_reaction_time_mapping(self) -> None:
        self.assertAlmostEqual(odds_to_reaction_time(4, 200, 80, 0.5), 240)


class TestAndersonMilsonSimulation(unittest.TestCase):
    def test_small_simulation_is_valid_and_revival_boosts_probability(self) -> None:
        params = AndersonMilsonParameters(
            desirability_shape=0.199,
            desirability_scale=0.482,
            mean_decay=4.076,
            mean_revival_interval=50,
        )
        simulation = simulate_anderson_milson(
            50, 100, params, decay_kind="power", rng=123
        )
        self.assertEqual(simulation.probabilities.shape, (100, 50))
        self.assertTrue(np.all(simulation.probabilities >= 0))
        self.assertTrue(np.all(simulation.probabilities <= 1))
        at_revival = simulation.revivals
        self.assertTrue(np.any(at_revival))
        self.assertGreater(
            float(np.mean(simulation.probabilities[at_revival])),
            float(np.mean(simulation.probabilities)),
        )

    def test_released_probability_mapping_and_posthoc_scale(self) -> None:
        params = AndersonMilsonParameters(
            desirability_shape=0.5,
            desirability_scale=0.2,
            mean_decay=1.0,
            mean_revival_interval=20,
        )
        intended = simulate_anderson_milson(
            5, 10, params, decay_kind="exponential", rng=9
        )
        released = simulate_anderson_milson(
            5,
            10,
            params,
            decay_kind="exponential",
            occurrence_mapping="released_probability",
            rng=9,
        )
        np.testing.assert_allclose(intended.odds, released.odds)
        np.testing.assert_allclose(
            intended.probabilities, odds_to_probability(intended.odds)
        )
        np.testing.assert_allclose(
            intended.prediction_values, intended.probabilities
        )
        np.testing.assert_allclose(
            released.probabilities, np.minimum(released.odds, 1)
        )
        np.testing.assert_allclose(released.prediction_values, released.odds)
        self.assertAlmostEqual(
            scale_anderson_milson_predictions(0.5, 2, semantics="odds"),
            2 / 3,
        )
        self.assertAlmostEqual(
            scale_anderson_milson_predictions(
                0.25, 0.5, semantics="released_probability"
            ),
            0.125,
        )
        self.assertAlmostEqual(
            scale_anderson_milson_predictions(
                0.75, 2, semantics="released_probability"
            ),
            1.5,
        )
        self.assertAlmostEqual(
            released_geometric_mean_scale([0.2, 0.8], [0.1, 0.4]), 2
        )

    def test_conditional_prediction_stage(self) -> None:
        simulation = AndersonMilsonSimulation(
            occurrences=np.array([[1], [0], [1], [0]], dtype=bool),
            probabilities=np.array([[0.1], [0.2], [0.3], [0.4]]),
            prediction_values=np.array([[0.1], [0.2], [0.3], [0.4]]),
            odds=np.array([[1 / 9], [0.25], [3 / 7], [2 / 3]]),
            elapsed_since_revival=np.ones((4, 1)),
            revivals=np.zeros((4, 1), dtype=bool),
            initial_desirability=np.array([1.0]),
            decay=np.array([1.0]),
            occurrence_mapping="odds",
        )
        result = anderson_milson_conditional_predictions(
            simulation, history_length=2
        )
        np.testing.assert_array_equal(
            result.history_summaries, np.array([[1, 1, 0], [1, 2, 0]])
        )
        np.testing.assert_allclose(result.mean_prediction_values, [0.4, 0.3])
        np.testing.assert_array_equal(result.counts, [1, 1])


if __name__ == "__main__":
    unittest.main()

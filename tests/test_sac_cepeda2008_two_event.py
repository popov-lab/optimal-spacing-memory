"""Tests for the two-event, session-specific Cepeda (2008) SAC fit."""

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from fit_sac_cepeda2008_two_event import (  # noqa: E402
    ZERO_ISI_DAYS,
    continuous_optimum_isi,
    forgetting_function,
    forgetting_strength_from_log_tau,
    load_data,
    residual_vector,
    spacing_strength,
)


class TestCepeda2008TwoEventData(unittest.TestCase):
    def test_loader_includes_both_functions(self) -> None:
        data = load_data()
        self.assertEqual(len(data.forgetting_rows), 11)
        self.assertEqual(len(data.spacing_rows), 26)
        self.assertTrue(
            all(row["function"] == "forgetting" for row in data.forgetting_rows)
        )
        self.assertTrue(
            all(row["function"] == "spacing" for row in data.spacing_rows)
        )

    def test_nominal_zero_is_corrected_only_for_modeling(self) -> None:
        data = load_data()
        self.assertEqual(data.observed_lag[0], 0.0)
        self.assertEqual(data.model_lag[0], ZERO_ISI_DAYS)
        spacing_zero = data.observed_isi == 0.0
        self.assertEqual(int(spacing_zero.sum()), 4)
        np.testing.assert_allclose(data.model_isi[spacing_zero], ZERO_ISI_DAYS)
        np.testing.assert_allclose(
            data.model_isi[~spacing_zero], data.observed_isi[~spacing_zero]
        )


class TestCepeda2008TwoEventModel(unittest.TestCase):
    def test_spacing_formula_matches_explicit_trace_sum(self) -> None:
        isi = 11.0
        ri = 35.0
        delta_2 = 0.37
        d = 0.18
        tau = 0.06
        strength = float(spacing_strength(isi, ri, delta_2, d, tau))
        f_isi = float(forgetting_function(isi, d, tau))
        f_ri = float(forgetting_function(ri, d, tau))
        f_total = float(forgetting_function(isi + ri, d, tau))
        increment_1 = 1.0
        increment_2 = delta_2 * (1.0 - increment_1 * f_isi)
        explicit = increment_1 * f_total + increment_2 * f_ri
        self.assertAlmostEqual(strength, explicit, places=13)

    def test_no_decay_leaves_no_capacity_for_second_study(self) -> None:
        strength = float(
            spacing_strength(isi=21.0, ri=70.0, delta_2=0.3, d=0.0, tau=1.0)
        )
        self.assertAlmostEqual(strength, 1.0)

    def test_forgetting_strength_does_not_depend_on_delta_2(self) -> None:
        lags = np.array([ZERO_ISI_DAYS, 1.0, 35.0])
        expected = forgetting_function(lags, d=0.2, tau=0.3)
        actual = forgetting_strength_from_log_tau(
            lags, d=0.2, log_tau=float(np.log(0.3))
        )
        np.testing.assert_allclose(actual, expected)

    def test_protocols_use_the_requested_observations(self) -> None:
        data = load_data()
        z = np.array(
            [0.0, np.log(0.2), np.log(0.1), 0.3, np.log(0.05)]
        )
        self.assertEqual(residual_vector(z, data, "spacing_only").size, 26)
        self.assertEqual(residual_vector(z, data, "joint").size, 37)

    def test_closed_form_optimum_is_a_local_maximum(self) -> None:
        d = 0.16
        tau = 0.08
        ri = 70.0
        delta_2 = 0.45
        optimum = continuous_optimum_isi(d, tau, ri, delta_2)
        self.assertGreater(optimum, 0.0)
        center = float(spacing_strength(optimum, ri, delta_2, d, tau))
        left = float(spacing_strength(0.99 * optimum, ri, delta_2, d, tau))
        right = float(spacing_strength(1.01 * optimum, ri, delta_2, d, tau))
        self.assertGreater(center, left)
        self.assertGreater(center, right)


if __name__ == "__main__":
    unittest.main()

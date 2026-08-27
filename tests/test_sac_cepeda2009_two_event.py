"""Tests for the two-event, session-specific Cepeda (2009) SAC fit."""

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from fit_sac_cepeda2009_two_event import (  # noqa: E402
    PANEL_KEYS,
    SHORTEST_ISI_DAYS,
    continuous_optimum_isi_from_log_tau,
    forgetting_function,
    load_panels,
    residual_vector,
    spacing_strength,
)


class TestCepeda2009TwoEventData(unittest.TestCase):
    def test_loader_has_six_points_per_function_and_panel(self) -> None:
        panels = load_panels()
        self.assertEqual(tuple(panels), PANEL_KEYS)
        for key in PANEL_KEYS:
            panel = panels[key]
            self.assertEqual(len(panel.forgetting_rows), 6)
            self.assertEqual(len(panel.spacing_rows), 6)

    def test_nominal_zero_corrections_match_the_experiments(self) -> None:
        panels = load_panels()
        for key in PANEL_KEYS:
            panel = panels[key]
            self.assertEqual(panel.observed_lag[0], 0.0)
            self.assertEqual(panel.observed_isi[0], 0.0)
            self.assertEqual(panel.model_lag[0], SHORTEST_ISI_DAYS[key])
            self.assertEqual(panel.model_isi[0], SHORTEST_ISI_DAYS[key])
        self.assertAlmostEqual(
            24.0 * 60.0 * SHORTEST_ISI_DAYS["a"], 5.0
        )
        self.assertAlmostEqual(
            24.0 * 60.0 * SHORTEST_ISI_DAYS["b"], 20.0
        )


class TestCepeda2009TwoEventModel(unittest.TestCase):
    def test_spacing_formula_matches_explicit_trace_sum(self) -> None:
        isi = 28.0
        ri = 168.0
        delta_2 = 0.42
        d = 0.19
        tau = 0.7
        strength = float(spacing_strength(isi, ri, delta_2, d, tau))
        f_isi = float(forgetting_function(isi, d, tau))
        f_ri = float(forgetting_function(ri, d, tau))
        f_total = float(forgetting_function(isi + ri, d, tau))
        increment_2 = delta_2 * (1.0 - f_isi)
        self.assertAlmostEqual(
            strength, f_total + increment_2 * f_ri, places=13
        )

    def test_no_decay_leaves_no_capacity_for_second_study(self) -> None:
        strength = float(
            spacing_strength(isi=28.0, ri=168.0, delta_2=0.3, d=0.0, tau=1.0)
        )
        self.assertAlmostEqual(strength, 1.0)

    def test_new_rule_changes_only_the_isi_independent_offset(self) -> None:
        delta_2 = 0.37
        d = 0.21
        tau = 0.4
        ri = 168.0
        f_ri = float(forgetting_function(ri, d, tau))
        expected_difference = (delta_2 - 1.0) * f_ri
        for isi in (0.02, 1.0, 28.0, 168.0):
            f_isi = float(forgetting_function(isi, d, tau))
            f_total = float(forgetting_function(isi + ri, d, tau))
            old_normalized = f_total + f_ri - delta_2 * f_isi * f_ri
            new = float(spacing_strength(isi, ri, delta_2, d, tau))
            self.assertAlmostEqual(new - old_normalized, expected_difference)

    def test_old_and_new_rules_have_the_same_numeric_isi_slope(self) -> None:
        delta_2 = 0.44
        d = 0.18
        tau = 0.6
        ri = 168.0
        isi = 28.0
        step = 1e-4

        def old_normalized(value: float) -> float:
            f_value = float(forgetting_function(value, d, tau))
            f_ri = float(forgetting_function(ri, d, tau))
            f_total = float(forgetting_function(value + ri, d, tau))
            return f_total + f_ri - delta_2 * f_value * f_ri

        old_slope = (
            old_normalized(isi + step) - old_normalized(isi - step)
        ) / (2.0 * step)
        new_slope = (
            float(spacing_strength(isi + step, ri, delta_2, d, tau))
            - float(spacing_strength(isi - step, ri, delta_2, d, tau))
        ) / (2.0 * step)
        self.assertAlmostEqual(new_slope, old_slope, places=10)

    def test_protocols_use_the_requested_observations(self) -> None:
        panels = load_panels()
        z = np.array(
            [
                np.log(0.2),
                np.log(0.2),
                np.log(0.2),
                np.log(0.5),
                np.log(0.5),
                np.log(0.5),
                0.3,
                np.log(0.05),
                0.0,
            ]
        )
        self.assertEqual(residual_vector(z, panels, "spacing_only").size, 18)
        self.assertEqual(residual_vector(z, panels, "joint").size, 36)

    def test_closed_form_optimum_is_a_local_maximum(self) -> None:
        d = 0.17
        tau = 0.8
        log_tau = float(np.log(tau))
        ri = 168.0
        delta_2 = 0.45
        optimum = continuous_optimum_isi_from_log_tau(
            d, log_tau, ri, delta_2
        )
        self.assertGreater(optimum, 0.0)
        center = float(spacing_strength(optimum, ri, delta_2, d, tau))
        left = float(spacing_strength(0.99 * optimum, ri, delta_2, d, tau))
        right = float(spacing_strength(1.01 * optimum, ri, delta_2, d, tau))
        self.assertGreater(center, left)
        self.assertGreater(center, right)


if __name__ == "__main__":
    unittest.main()

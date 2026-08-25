"""Tests for the time-resolved Cepeda (2008) SAC implementation."""

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from fit_sac_cepeda2008_four_event import (  # noqa: E402
    WITHIN_SESSION_GAP_DAYS,
    ZERO_ISI_DAYS,
    four_event_times,
    load_spacing_data,
    sac_forgetting,
    sac_strength_for_schedule,
    sac_strength_from_times,
)


class TestCepeda2008Data(unittest.TestCase):
    def test_loader_excludes_new_forgetting_observations(self) -> None:
        data = load_spacing_data()
        self.assertEqual(len(data.rows), 26)
        self.assertTrue(all(row["function"] == "spacing" for row in data.rows))
        self.assertTrue(all(row["panel"] == "d" for row in data.rows))
        self.assertTrue(np.all(data.ri > 0.0))

    def test_nominal_zero_isi_is_corrected_only_for_modeling(self) -> None:
        data = load_spacing_data()
        zero = data.observed_isi == 0.0
        self.assertEqual(int(zero.sum()), 4)
        np.testing.assert_allclose(data.model_isi[zero], ZERO_ISI_DAYS)
        np.testing.assert_allclose(
            data.model_isi[~zero], data.observed_isi[~zero]
        )


class TestSACRecursion(unittest.TestCase):
    def test_two_events_equal_the_explicit_general_form(self) -> None:
        isi = 5.0
        ri = 35.0
        delta = 0.42
        d = 0.17
        tau = 0.08
        recursive = sac_strength_from_times(
            np.array([0.0, isi]), isi + ri, delta, d, tau
        )
        f_isi = float(sac_forgetting(isi, d, tau))
        f_ri = float(sac_forgetting(ri, d, tau))
        f_total = float(sac_forgetting(isi + ri, d, tau))
        explicit = (
            delta * f_ri
            + delta * f_total
            - delta**2 * f_isi * f_ri
        )
        self.assertAlmostEqual(recursive, explicit, places=13)

    def test_no_decay_reduces_to_repeated_delta_rule_learning(self) -> None:
        delta = 0.3
        strength = sac_strength_from_times(
            np.array([0.0, 1.0, 2.0, 3.0]),
            4.0,
            delta,
            d=0.0,
            tau=1.0,
        )
        self.assertAlmostEqual(strength, 1.0 - (1.0 - delta) ** 4)

    def test_four_event_schedule_matches_supplied_absolute_times(self) -> None:
        isi = 7.0
        ri = 35.0
        events, test_time = four_event_times(isi, ri)
        h = WITHIN_SESSION_GAP_DAYS
        np.testing.assert_allclose(events, [0.0, h, h + isi, 2.0 * h + isi])
        self.assertAlmostEqual(test_time, 3.0 * h + isi + ri)
        np.testing.assert_allclose(np.diff(events), [h, isi, h])
        self.assertAlmostEqual(test_time - events[-1], h + ri)

    def test_batched_schedule_strength_is_finite_and_bounded(self) -> None:
        strengths = sac_strength_for_schedule(
            np.array([ZERO_ISI_DAYS, 1.0, 105.0]),
            np.array([7.0, 35.0, 350.0]),
            delta=0.4,
            d=0.15,
            tau=0.03,
            schedule="four_event",
        )
        self.assertTrue(np.all(np.isfinite(strengths)))
        self.assertTrue(np.all(strengths >= 0.0))
        self.assertTrue(np.all(strengths <= 1.0))

    def test_small_tau_evaluation_remains_finite(self) -> None:
        values = sac_forgetting(
            np.array([0.0, 1e-6, 1.0, 350.0]), d=0.1, tau=1e-300
        )
        self.assertEqual(values[0], 1.0)
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertTrue(np.all(values >= 0.0))


if __name__ == "__main__":
    unittest.main()

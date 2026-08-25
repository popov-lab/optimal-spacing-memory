"""Integration checks for the Anderson et al. CSV preparation pipeline."""

from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_anderson_2023_data import csv_value, main  # noqa: E402


SOURCE = ROOT / "external" / "anderson_2023" / "matlab"


class TestCsvConventions(unittest.TestCase):
    def test_booleans_are_not_serialized_as_integers(self) -> None:
        self.assertEqual(csv_value(True), "true")
        self.assertEqual(csv_value(False), "false")


@unittest.skipUnless(SOURCE.exists(), "authors' MATLAB release not present")
class TestModelingTables(unittest.TestCase):
    def test_small_and_medium_tables_export_and_validate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            result = main(
                [
                    "--source-dir", str(SOURCE),
                    "--output-dir", str(output),
                    "--skip-corpus",
                    "--skip-exact-history",
                    "--skip-appendix-c",
                ]
            )
            self.assertEqual(result, 0)
            with (output / "manifest.csv").open(encoding="utf-8", newline="") as handle:
                manifest = {row["table_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(len(manifest), 21)
            expected_rows = {
                "environmental_fit_targets": 672,
                "frequency_recency_bins": 3072,
                "behavioral_conditions": 353,
                "behavioral_schedule_variants": 409,
                "behavioral_presentations": 2010,
                "range_spacing_histories": 128332,
                "range_frequency_histories": 582736,
                "published_environmental_surfaces": 6048,
            }
            for table_id, rows in expected_rows.items():
                self.assertEqual(int(manifest[table_id]["rows"]), rows)
            with (output / "environmental_model_parameters.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                parameters = list(csv.DictReader(handle))
            alpha_models = {
                row["model_id"]
                for row in parameters
                if row["parameter_name"] == "alpha"
            }
            self.assertEqual(
                alpha_models,
                {"gpe", "actr", "pavlik_anderson", "ppe", "mcm", "ampe"},
            )
            self.assertEqual(
                main(["--output-dir", str(output), "--validate-only"]),
                0,
            )
            conflict = output / "corpus_sources.csv"
            conflict.write_text("placeholder\n", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "mixed snapshot"):
                main(
                    [
                        "--source-dir", str(SOURCE),
                        "--output-dir", str(output),
                        "--skip-corpus", "--skip-exact-history", "--skip-appendix-c",
                    ]
                )
            conflict.unlink()
            extra = output / "unexpected.csv"
            extra.write_text("x\n1\n", encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "inventory differs"):
                main(["--output-dir", str(output), "--validate-only"])
            extra.unlink()
            interrupted = output / ".corpus_occurrences.csv.partial"
            interrupted.write_text("incomplete\n", encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "interrupted-write artifacts"):
                main(["--output-dir", str(output), "--validate-only"])
            interrupted.unlink()
            manifest_path = output / "manifest.csv"
            with manifest_path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                fieldnames = reader.fieldnames
                rows = list(reader)
            self.assertIsNotNone(fieldnames)
            rows[0]["rows"] = str(int(rows[0]["rows"]) + 1)
            with manifest_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(AssertionError, "Row count mismatch"):
                main(["--output-dir", str(output), "--validate-only"])


if __name__ == "__main__":
    unittest.main()

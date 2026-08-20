"""Tests for the dependency-light MATLAB 7.3 reader."""

from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from matlab_v73_reader import (  # noqa: E402
    MatlabV73File,
    iter_reddit_cells,
    read_reddit_metadata,
)


REDDIT_DATA = (
    ROOT / "external" / "anderson_2023" / "matlab" / "redditData.mat"
)


@unittest.skipUnless(REDDIT_DATA.exists(), "authors' redditData.mat not present")
class TestAndersonRedditData(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.mat = MatlabV73File(REDDIT_DATA)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.mat.close()

    def test_top_level_schema(self) -> None:
        self.assertEqual(
            self.mat.variable_names,
            (
                "arraysApr23",
                "arraysMay5",
                "countsApr23",
                "countsMay5",
                "subredditsApr23",
                "subredditsMay5",
                "vocabApr23",
                "vocabMay5",
            ),
        )
        self.assertEqual(self.mat.info("arraysApr23").matlab_shape, (501, 1))
        self.assertEqual(self.mat.info("countsApr23").matlab_shape, (20000, 1))

    def test_cell_references_and_matlab_axis_order(self) -> None:
        references = self.mat.references("arraysApr23")
        self.assertEqual(references.shape, (501,))
        self.assertTrue(np.all(references > 0))

        first = self.mat.read_cell("arraysApr23", 0)
        self.assertEqual(first.shape, (2394, 100))
        self.assertEqual(first.dtype, np.dtype("<f8"))
        np.testing.assert_array_equal(
            first[0, :6], [1, 3484, 76, 3696, 108, 321]
        )
        self.assertTrue(np.all(first == np.floor(first)))
        self.assertGreaterEqual(float(first.min()), 0)
        self.assertLessEqual(float(first.max()), 20000)

    def test_count_vectors(self) -> None:
        april = self.mat.read_dataset("countsApr23", matlab_order=True).ravel()
        may = self.mat.read_dataset("countsMay5", matlab_order=True).ravel()
        self.assertEqual(april.shape, (20000,))
        self.assertEqual(may.shape, (20000,))
        self.assertEqual(float(april.sum()), 8_351_677)
        self.assertEqual(float(may.sum()), 8_177_928)
        self.assertTrue(np.all(april >= 0))
        self.assertTrue(np.all(may >= 0))

    def test_string_cells(self) -> None:
        self.assertEqual(
            self.mat.read_string_cells("subredditsApr23")[:3],
            ["funny", "AskReddit", "gaming"],
        )
        self.assertEqual(
            self.mat.read_string_cells("vocabApr23")[:5],
            ["i", "you", "its", "like", "just"],
        )

    def test_streaming_api(self) -> None:
        stream = iter_reddit_cells(REDDIT_DATA, "arraysApr23")
        index, first = next(stream)
        stream.close()
        self.assertEqual(index, 0)
        self.assertIsInstance(first, np.ndarray)
        self.assertEqual(first.shape, (2394, 100))

        strings = iter_reddit_cells(REDDIT_DATA, "vocabApr23")
        self.assertEqual(next(strings), (0, "i"))
        strings.close()

    def test_partition_metadata_without_comment_loading(self) -> None:
        metadata = read_reddit_metadata(REDDIT_DATA)
        april = metadata["2021-04-23"]
        may = metadata["2021-05-05"]
        self.assertEqual(april.source_count, 501)
        self.assertEqual(may.source_count, 501)
        self.assertEqual(april.event_count, 573_836)
        self.assertEqual(may.event_count, 559_346)
        self.assertEqual(april.eligible_source_count, 223)
        self.assertEqual(may.eligible_source_count, 216)
        self.assertEqual(april.event_width, 100)
        self.assertEqual(may.event_width, 100)
        self.assertEqual(april.vocabulary_size, 20_000)
        self.assertEqual(may.vocabulary_size, 20_000)
        self.assertEqual(april.retained_item_count, 8_351_677)
        self.assertEqual(may.retained_item_count, 8_177_928)

    def test_streamed_comment_ids_reconcile_to_released_counts(self) -> None:
        for arrays_name, counts_name in (
            ("arraysApr23", "countsApr23"),
            ("arraysMay5", "countsMay5"),
        ):
            observed = np.zeros(20_001, dtype=np.int64)
            for comments in self.mat.iter_cells(arrays_name):
                self.assertEqual(comments.shape[1], 100)
                self.assertTrue(np.all(comments == np.floor(comments)))
                positive = comments > 0
                values = comments[positive].astype(np.int64)
                observed += np.bincount(values, minlength=20_001)

                # Zeros are padding only; retained IDs never resume after the
                # first zero in an event.
                seen_zero = np.maximum.accumulate(comments == 0, axis=1)
                self.assertFalse(np.any(seen_zero & positive))

                # The released comment rows contain unique term IDs.
                sorted_rows = np.sort(comments, axis=1)
                duplicates = (sorted_rows[:, 1:] == sorted_rows[:, :-1]) & (
                    sorted_rows[:, 1:] != 0
                )
                self.assertFalse(np.any(duplicates))

            released = self.mat.read_dataset(counts_name).ravel().astype(np.int64)
            np.testing.assert_array_equal(observed[1:], released)


if __name__ == "__main__":
    unittest.main()

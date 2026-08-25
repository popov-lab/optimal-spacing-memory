#!/usr/bin/env python3
"""Create deterministic, analysis-ready CSVs from Anderson et al. (2023).

The author release mixes raw corpus streams, binned sufficient statistics,
compressed history cells, laboratory condition means, and cached predictions.
This exporter keeps those layers separate and preserves MATLAB's 1-based
indices wherever an index is part of the released representation.

The script never modifies the vendored release.  CSV files are written
atomically, validated while they are produced, and summarized in a manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import shutil
import sys
import tempfile
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, TextIO

import numpy as np
from scipy.io import loadmat


SCHEMA_VERSION = "1.0.0"
MANIFEST_FIELDS = [
    "schema_version", "table_id", "filename", "rows", "bytes", "sha256",
    "source_files", "notes",
]

KNOWN_GENERATED_CSVS = frozenset(
    {
        "appendix_c_decay_replicates.csv",
        "appendix_c_frequency_distributions.csv",
        "appendix_c_no_revival_trials.csv",
        "appendix_c_revival_pairs.csv",
        "behavioral_conditions.csv",
        "behavioral_experiments.csv",
        "behavioral_fit_starts.csv",
        "behavioral_plot_points.csv",
        "behavioral_presentations.csv",
        "behavioral_schedule_variants.csv",
        "bin_definitions.csv",
        "corpus_occurrences.csv",
        "corpus_sources.csv",
        "corpus_texts.csv",
        "corpus_vocabulary.csv",
        "environmental_fit_targets.csv",
        "environmental_model_parameters.csv",
        "exact_history_cells.csv",
        "frequency_recency_bins.csv",
        "hick_frequency_recency_surface.csv",
        "hick_micro_fit.csv",
        "hick_microenvironment_odds.csv",
        "hick_observations.csv",
        "published_environmental_surfaces.csv",
        "range_frequency_histories.csv",
        "range_spacing_histories.csv",
        "spacing_figure_cells.csv",
        "spacing_recency_bins.csv",
        "three_occurrence_range_bins.csv",
        "three_occurrence_range_groups.csv",
    }
)
CORPUS_CSVS = frozenset(
    {"corpus_occurrences.csv", "corpus_sources.csv", "corpus_texts.csv", "corpus_vocabulary.csv"}
)
APPENDIX_C_CSVS = frozenset(
    {
        "appendix_c_decay_replicates.csv",
        "appendix_c_frequency_distributions.csv",
        "appendix_c_no_revival_trials.csv",
        "appendix_c_revival_pairs.csv",
    }
)


@dataclass(frozen=True)
class WrittenTable:
    filename: str
    table_id: str
    rows: int
    bytes: int
    sha256: str
    source_files: str
    notes: str = ""


class TableWriter:
    """RFC-4180-compatible CSV writer with atomic replacement and row count."""

    def __init__(
        self,
        output_dir: Path,
        filename: str,
        fieldnames: Sequence[str],
        *,
        table_id: str,
        source_files: str,
        notes: str = "",
    ) -> None:
        self.output_dir = output_dir
        self.filename = filename
        self.fieldnames = list(fieldnames)
        self.table_id = table_id
        self.source_files = source_files
        self.notes = notes
        self.rows = 0
        self._tmp_path: Path | None = None
        self._file: TextIO | None = None
        self._writer: csv.DictWriter[str] | None = None

    def __enter__(self) -> "TableWriter":
        self.output_dir.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(
            dir=self.output_dir, prefix=f".{self.filename}.", suffix=".tmp"
        )
        os.close(fd)
        self._tmp_path = Path(name)
        self._file = self._tmp_path.open("w", encoding="utf-8", newline="")
        self._writer = csv.DictWriter(
            self._file,
            fieldnames=self.fieldnames,
            extrasaction="raise",
            lineterminator="\n",
        )
        self._writer.writeheader()
        return self

    def writerow(self, row: dict[str, Any]) -> None:
        assert self._writer is not None
        self._writer.writerow({key: csv_value(row.get(key)) for key in self.fieldnames})
        self.rows += 1

    def writerows(self, rows: Iterable[dict[str, Any]]) -> None:
        for row in rows:
            self.writerow(row)

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        assert self._file is not None and self._tmp_path is not None
        self._file.flush()
        os.fsync(self._file.fileno())
        self._file.close()
        if exc_type is None:
            self._tmp_path.replace(self.output_dir / self.filename)
        else:
            self._tmp_path.unlink(missing_ok=True)

    def record(self) -> WrittenTable:
        path = self.output_dir / self.filename
        return WrittenTable(
            filename=self.filename,
            table_id=self.table_id,
            rows=self.rows,
            bytes=path.stat().st_size,
            sha256=sha256_file(path),
            source_files=self.source_files,
            notes=self.notes,
        )


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    # bool is a subclass of int, so this check must precede the integer case.
    if isinstance(value, (np.bool_, bool)):
        return "true" if bool(value) else "false"
    if isinstance(value, (np.floating, float)):
        x = float(value)
        if math.isnan(x):
            return ""
        if math.isinf(x):
            raise ValueError("CSV export encountered an infinite value")
        return format(x, ".17g")
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, str) and ("\n" in value or "\r" in value):
        raise ValueError("CSV export encountered an embedded line break")
    return value


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def csv_file_stats(path: Path, chunk_size: int = 1 << 20) -> tuple[int, int, str]:
    """Return byte count, data-row count, and SHA-256 for a generated CSV.

    Generated fields contain no embedded CR/LF characters, so each record is
    exactly one LF-terminated physical line and the header consumes one line.
    """

    digest = hashlib.sha256()
    byte_count = 0
    line_count = 0
    last_byte = b""
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
            byte_count += len(chunk)
            line_count += chunk.count(b"\n")
            last_byte = chunk[-1:]
    if byte_count == 0 or last_byte != b"\n" or line_count < 1:
        raise AssertionError(f"Generated CSV is empty or not LF-terminated: {path}")
    return byte_count, line_count - 1, digest.hexdigest()


def verify_source_checksums(source: Path, required: Sequence[str]) -> None:
    """Verify every input against the immutable-release checksum manifest."""

    checksum_path = source.parent / "SHA256SUMS"
    if not checksum_path.is_file():
        raise FileNotFoundError(f"Missing release checksum manifest: {checksum_path}")
    expected: dict[str, str] = {}
    for line in checksum_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        expected[relative.removeprefix("./")] = digest
    for name in required:
        if name not in expected:
            raise AssertionError(f"No release checksum recorded for {name}")
    actual_files = {path.name for path in source.iterdir() if path.is_file()}
    if actual_files != set(expected):
        missing = sorted(set(expected) - actual_files)
        extra = sorted(actual_files - set(expected))
        raise AssertionError(
            f"Vendored release inventory differs from SHA256SUMS; missing={missing}, extra={extra}"
        )
    # Verify the whole immutable snapshot, not only the numeric arrays.  The
    # exporter derives some semantics from the accompanying MATLAB scripts.
    for name, expected_digest in sorted(expected.items()):
        observed = sha256_file(source / name)
        if observed != expected_digest:
            raise AssertionError(
                f"Source checksum mismatch for {name}: {observed} != {expected_digest}"
            )


def load_classic(path: Path) -> dict[str, Any]:
    return {
        key: value
        for key, value in loadmat(path, squeeze_me=True, struct_as_record=False).items()
        if not key.startswith("__")
    }


def matlab_cells(value: Any) -> list[Any]:
    array = np.asarray(value, dtype=object)
    return list(array.ravel(order="F"))


def matlab_vector(value: Any) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim == 0:
        return array.reshape(1)
    return array.ravel(order="F")


def matlab_datenum_to_iso(value: float) -> str:
    # MATLAB serial dates are 366 days ahead of Python's proleptic ordinal.
    day = int(value)
    fraction = float(value) - day
    dt = datetime.fromordinal(day) + timedelta(days=fraction) - timedelta(days=366)
    dt = dt.replace(tzinfo=timezone.utc)
    # Released values nominally have one-second precision; remove FP micro-noise.
    dt = (dt + timedelta(microseconds=500_000)).replace(microsecond=0)
    return dt.isoformat().replace("+00:00", "Z")


def bin_rows(bounds: np.ndarray, means: np.ndarray, scheme: str) -> Iterator[dict[str, Any]]:
    for index, midpoint in enumerate(means, start=1):
        yield {
            "bin_scheme": scheme,
            "bin_index": index,
            "lower_inclusive": int(bounds[index - 1]) + 1,
            "upper_inclusive": int(bounds[index]),
            "representative_value": midpoint,
        }


def safe_probability(hits: int | float, trials: int | float) -> float | None:
    return float(hits) / float(trials) if trials else None


def export_bin_definitions(source: Path, output: Path) -> WrittenTable:
    data = load_classic(source / "base32.mat")
    bounds = matlab_vector(data["bounds"])
    means = matlab_vector(data["means"])
    bounds5 = matlab_vector(data["bounds5"])
    if len(bounds) != 33 or len(means) != 32 or len(bounds5) != 6:
        raise AssertionError("Unexpected base32 bin shapes")
    writer = TableWriter(
        output,
        "bin_definitions.csv",
        ["bin_scheme", "bin_index", "lower_inclusive", "upper_inclusive", "representative_value"],
        table_id="bin_definitions",
        source_files="base32.mat",
    )
    with writer:
        writer.writerows(bin_rows(bounds, means, "square32"))
        mid5 = (bounds5[:-1] + 1 + bounds5[1:]) / 2
        writer.writerows(bin_rows(bounds5, mid5, "spacing5"))
    return writer.record()


def corpus_aggregate_arrays(source: Path) -> dict[str, dict[str, np.ndarray]]:
    combined = load_classic(source / "Combined.mat")
    twitter = load_classic(source / "twitterData.mat")
    result: dict[str, dict[str, np.ndarray]] = {}
    for corpus, data in (("combined", combined), ("twitter", twitter)):
        result[corpus] = {
            name: np.asarray(data[name]) for name in ("hitsA", "countsA", "hits2", "counts2")
        }
    result["reddit_derived"] = {
        name: result["combined"][name].astype(np.int64)
        - result["twitter"][name].astype(np.int64)
        for name in ("hitsA", "countsA", "hits2", "counts2")
    }
    if any(np.any(value < 0) for value in result["reddit_derived"].values()):
        raise AssertionError("Combined - Twitter yielded a negative Reddit aggregate cell")
    return result


def export_environment_bins(source: Path, output: Path) -> list[WrittenTable]:
    base = load_classic(source / "base32.mat")
    bounds = matlab_vector(base["bounds"])
    means = matlab_vector(base["means"])
    aggregates = corpus_aggregate_arrays(source)
    records: list[WrittenTable] = []

    fr_writer = TableWriter(
        output,
        "frequency_recency_bins.csv",
        [
            "corpus", "recency_bin", "recency_lower", "recency_upper", "recency_midpoint",
            "frequency_bin", "frequency_lower", "frequency_upper", "frequency_midpoint",
            "hits", "trials", "probability", "in_fitted_frequency_domain", "meets_fit_count_threshold",
        ],
        table_id="frequency_recency_bins",
        source_files="Combined.mat;twitterData.mat;base32.mat",
        notes="reddit_derived is exactly Combined minus Twitter because redditData.mat omits these summaries",
    )
    with fr_writer:
        for corpus, arrays in aggregates.items():
            for fbin in range(32):
                for rbin in range(32):
                    hits = int(arrays["hitsA"][rbin, fbin])
                    trials = int(arrays["countsA"][rbin, fbin])
                    fr_writer.writerow({
                        "corpus": corpus,
                        "recency_bin": rbin + 1,
                        "recency_lower": int(bounds[rbin]) + 1,
                        "recency_upper": int(bounds[rbin + 1]),
                        "recency_midpoint": means[rbin],
                        "frequency_bin": fbin + 1,
                        "frequency_lower": int(bounds[fbin]) + 1,
                        "frequency_upper": int(bounds[fbin + 1]),
                        "frequency_midpoint": means[fbin],
                        "hits": hits,
                        "trials": trials,
                        "probability": safe_probability(hits, trials),
                        "in_fitted_frequency_domain": fbin < 15,
                        "meets_fit_count_threshold": fbin < 15 and trials >= 5000,
                    })
    records.append(fr_writer.record())

    spacing_writer = TableWriter(
        output,
        "spacing_recency_bins.csv",
        [
            "corpus", "recency_bin", "recency_lower", "recency_upper", "recency_midpoint",
            "spacing_bin", "spacing_lower", "spacing_upper", "spacing_midpoint",
            "hits", "trials", "probability",
        ],
        table_id="spacing_recency_bins",
        source_files="Combined.mat;twitterData.mat;base32.mat",
        notes="reddit_derived is exactly Combined minus Twitter",
    )
    with spacing_writer:
        for corpus, arrays in aggregates.items():
            for sbin in range(32):
                for rbin in range(32):
                    hits = int(arrays["hits2"][rbin, sbin])
                    trials = int(arrays["counts2"][rbin, sbin])
                    spacing_writer.writerow({
                        "corpus": corpus,
                        "recency_bin": rbin + 1,
                        "recency_lower": int(bounds[rbin]) + 1,
                        "recency_upper": int(bounds[rbin + 1]),
                        "recency_midpoint": means[rbin],
                        "spacing_bin": sbin + 1,
                        "spacing_lower": int(bounds[sbin]) + 1,
                        "spacing_upper": int(bounds[sbin + 1]),
                        "spacing_midpoint": means[sbin],
                        "hits": hits,
                        "trials": trials,
                        "probability": safe_probability(hits, trials),
                    })
    records.append(spacing_writer.record())

    grouped = [(0, 0), (1, 2), (3, 6), (7, 14), (15, 31)]
    figure_writer = TableWriter(
        output,
        "spacing_figure_cells.csv",
        [
            "corpus", "recency_bin", "recency_lower", "recency_upper", "recency_midpoint",
            "condition_index", "condition_type", "condition_lower", "condition_upper",
            "hits", "trials", "probability", "fit_included",
        ],
        table_id="spacing_figure_cells",
        source_files="Combined.mat;twitterData.mat;base32.mat;displayPair.m",
    )
    with figure_writer:
        for corpus, arrays in aggregates.items():
            for rbin in range(32):
                h = int(arrays["hitsA"][rbin, 0])
                n = int(arrays["countsA"][rbin, 0])
                figure_writer.writerow({
                    "corpus": corpus, "recency_bin": rbin + 1,
                    "recency_lower": int(bounds[rbin]) + 1,
                    "recency_upper": int(bounds[rbin + 1]), "recency_midpoint": means[rbin],
                    "condition_index": 1, "condition_type": "frequency_one_reference",
                    "condition_lower": 1, "condition_upper": 1, "hits": h, "trials": n,
                    "probability": safe_probability(h, n), "fit_included": n >= 5000,
                })
                for group_index, (lo, hi) in enumerate(grouped, start=2):
                    h = int(np.sum(arrays["hits2"][rbin, lo : hi + 1], dtype=np.int64))
                    n = int(np.sum(arrays["counts2"][rbin, lo : hi + 1], dtype=np.int64))
                    figure_writer.writerow({
                        "corpus": corpus, "recency_bin": rbin + 1,
                        "recency_lower": int(bounds[rbin]) + 1,
                        "recency_upper": int(bounds[rbin + 1]), "recency_midpoint": means[rbin],
                        "condition_index": group_index, "condition_type": "twice_occurring_spacing",
                        "condition_lower": int(bounds[lo]) + 1,
                        "condition_upper": int(bounds[hi + 1]), "hits": h, "trials": n,
                        "probability": safe_probability(h, n), "fit_included": n > 0,
                    })
    records.append(figure_writer.record())

    combined = load_classic(source / "Combined.mat")
    lags = np.asarray(combined["lags"], dtype=float)
    reconstructed = np.full((32, 21), np.nan, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        reconstructed[:, :15] = (
            aggregates["combined"]["hitsA"][:, :15]
            / aggregates["combined"]["countsA"][:, :15]
        )
        reconstructed[:, :15][aggregates["combined"]["countsA"][:, :15] < 5000] = np.nan
        reconstructed[:, 15] = reconstructed[:, 0]
        for group, (lo, hi) in enumerate(grouped, start=16):
            reconstructed[:, group] = (
                np.sum(aggregates["combined"]["hits2"][:, lo : hi + 1], axis=1)
                / np.sum(aggregates["combined"]["counts2"][:, lo : hi + 1], axis=1)
            )
    if not np.array_equal(reconstructed, lags, equal_nan=True):
        difference = np.nanmax(np.abs(reconstructed - lags))
        raise AssertionError(f"Released lags do not reconstruct exactly; max difference={difference}")
    target_writer = TableWriter(
        output,
        "environmental_fit_targets.csv",
        [
            "target_id", "matlab_linear_index", "panel", "recency_bin", "recency_lower",
            "recency_upper", "recency_midpoint", "condition_index", "condition_type",
            "condition_lower", "condition_upper", "probability", "fit_included",
        ],
        table_id="environmental_fit_targets",
        source_files="Combined.mat;base32.mat;displayPair.m",
    )
    with target_writer:
        for col in range(21):
            for row in range(32):
                if col < 15:
                    panel = "frequency_recency"
                    condition_type = "frequency"
                    lo, hi = int(bounds[col]) + 1, int(bounds[col + 1])
                elif col == 15:
                    panel = "spacing_recency"
                    condition_type = "frequency_one_reference"
                    lo = hi = 1
                else:
                    panel = "spacing_recency"
                    condition_type = "twice_occurring_spacing"
                    group_lo, group_hi = grouped[col - 16]
                    lo, hi = int(bounds[group_lo]) + 1, int(bounds[group_hi + 1])
                value = float(lags[row, col])
                target_writer.writerow({
                    "target_id": f"target_{row + 1:02d}_{col + 1:02d}",
                    "matlab_linear_index": row + 1 + 32 * col,
                    "panel": panel, "recency_bin": row + 1,
                    "recency_lower": int(bounds[row]) + 1,
                    "recency_upper": int(bounds[row + 1]), "recency_midpoint": means[row],
                    "condition_index": col + 1 if col < 15 else col - 14,
                    "condition_type": condition_type,
                    "condition_lower": lo, "condition_upper": hi,
                    "probability": value, "fit_included": np.isfinite(value),
                })
    if int(np.isfinite(lags).sum()) != 513:
        raise AssertionError("Expected 513 finite environmental fit targets")
    records.append(target_writer.record())
    return records


def export_exact_history_cells(source: Path, output: Path) -> WrittenTable:
    data = load_classic(source / "Combined.mat")
    counts = np.asarray(data["counts225"])
    if counts.shape != (1000, 1000, 225):
        raise AssertionError(f"Unexpected counts225 shape: {counts.shape}")
    writer = TableWriter(
        output,
        "exact_history_cells.csv",
        [
            "frequency", "most_recent_age", "second_most_recent_age",
            "raw_second_age_index", "spacing", "history_count",
        ],
        table_id="exact_history_cells",
        source_files="Combined.mat",
        notes="Sparse weighted cells; for N>2 this is not a complete occurrence schedule and has no cell-level hit count",
    )
    total = 0
    nonzero = 0
    with writer:
        for frequency in range(1, 226):
            plane = counts[:, :, frequency - 1]
            recent_idx, second_idx = np.nonzero(plane)
            values = plane[recent_idx, second_idx]
            if frequency == 1 and np.any(second_idx != 0):
                raise AssertionError("Singleton counts225 cells must use second-age index 1 as a sentinel")
            if frequency > 1 and np.any(second_idx <= recent_idx):
                raise AssertionError("For N>1, every second-most-recent age must exceed the most-recent age")
            nonzero += len(values)
            total += int(np.sum(values, dtype=np.int64))
            for recent0, second0, weight in zip(recent_idx, second_idx, values, strict=True):
                recent = int(recent0) + 1
                raw_second = int(second0) + 1
                singleton = frequency == 1
                writer.writerow({
                    "frequency": frequency,
                    "most_recent_age": recent,
                    "second_most_recent_age": None if singleton else raw_second,
                    "raw_second_age_index": raw_second,
                    "spacing": None if singleton else raw_second - recent,
                    "history_count": int(weight),
                })
    if nonzero != 10_404_456 or total != 2_384_082_074:
        raise AssertionError(f"counts225 invariant failed: nnz={nonzero}, total={total}")
    return writer.record()


def export_three_occurrence(source: Path, output: Path) -> list[WrittenTable]:
    base = load_classic(source / "base32.mat")
    bounds = matlab_vector(base["bounds"])
    means = matlab_vector(base["means"])
    datasets = {
        "twitter": load_classic(source / "threesTwitter.mat"),
        "reddit": load_classic(source / "threesReddit.mat"),
    }
    records: list[WrittenTable] = []
    raw_writer = TableWriter(
        output,
        "three_occurrence_range_bins.csv",
        [
            "corpus", "recency_bin", "recency_lower", "recency_upper", "recency_midpoint",
            "range_bin", "range_lower", "range_upper", "range_midpoint", "range_convention",
            "hits", "trials", "probability",
        ],
        table_id="three_occurrence_range_bins",
        source_files="threesTwitter.mat;threesReddit.mat;base32.mat",
    )
    with raw_writer:
        for corpus, data in datasets.items():
            counts, hits = np.asarray(data["counts"]), np.asarray(data["hits"])
            for qbin in range(32):
                for rbin in range(32):
                    h, n = int(hits[rbin, qbin]), int(counts[rbin, qbin])
                    raw_writer.writerow({
                        "corpus": corpus, "recency_bin": rbin + 1,
                        "recency_lower": int(bounds[rbin]) + 1,
                        "recency_upper": int(bounds[rbin + 1]), "recency_midpoint": means[rbin],
                        "range_bin": qbin + 1, "range_lower": int(bounds[qbin]) + 1,
                        "range_upper": int(bounds[qbin + 1]), "range_midpoint": means[qbin],
                        "range_convention": "exclusive_last_minus_first",
                        "hits": h, "trials": n, "probability": safe_probability(h, n),
                    })
    records.append(raw_writer.record())

    range_groups = [(2, 4), (5, 16), (17, 49), (50, 225), (226, 1000)]
    group_writer = TableWriter(
        output,
        "three_occurrence_range_groups.csv",
        [
            "corpus", "recency_bin", "recency_lower", "recency_upper", "recency_midpoint",
            "range_group", "range_lower", "range_upper", "range_convention",
            "hits", "trials", "probability", "fit_included",
        ],
        table_id="three_occurrence_range_groups",
        source_files="threesTwitter.mat;threesReddit.mat;base32.mat;displayThree.m",
        notes="Grouped hit/trial denominators are reconstructed exactly from the released 32-bin arrays",
    )
    grouped_bin_slices = [(1, 1), (2, 3), (4, 6), (7, 14), (15, 31)]
    with group_writer:
        for corpus, data in datasets.items():
            groups5 = np.asarray(data["groups5"], dtype=float)
            counts, hits = np.asarray(data["counts"]), np.asarray(data["hits"])
            for group in range(5):
                for rbin in range(32):
                    lo, hi = range_groups[group]
                    bin_lo, bin_hi = grouped_bin_slices[group]
                    grouped_hits = int(np.sum(hits[rbin, bin_lo : bin_hi + 1], dtype=np.int64))
                    grouped_trials = int(np.sum(counts[rbin, bin_lo : bin_hi + 1], dtype=np.int64))
                    probability = safe_probability(grouped_hits, grouped_trials)
                    released_probability = groups5[rbin, group]
                    agrees = (
                        not np.isfinite(released_probability)
                        if probability is None
                        else np.isclose(probability, released_probability, rtol=0, atol=0)
                    )
                    if not agrees:
                        raise AssertionError(
                            f"Grouped three-occurrence probability mismatch for {corpus}, row {rbin + 1}, group {group + 1}"
                        )
                    group_writer.writerow({
                        "corpus": corpus, "recency_bin": rbin + 1,
                        "recency_lower": int(bounds[rbin]) + 1,
                        "recency_upper": int(bounds[rbin + 1]), "recency_midpoint": means[rbin],
                        "range_group": group + 1, "range_lower": lo, "range_upper": hi,
                        "range_convention": "exclusive_last_minus_first",
                        "hits": grouped_hits, "trials": grouped_trials,
                        "probability": probability, "fit_included": grouped_trials > 0,
                    })
    records.append(group_writer.record())
    return records


def export_range_histories(source: Path, output: Path) -> list[WrittenTable]:
    records: list[WrittenTable] = []
    lag_data = load_classic(source / "results3to5.mat")
    lag_writer = TableWriter(
        output,
        "range_spacing_histories.csv",
        [
            "history_id", "frequency", "gap_1", "gap_2", "gap_3", "gap_4",
            "range_exclusive", "future_count", "future_hit", "corpus",
        ],
        table_id="range_spacing_histories",
        source_files="results3to5.mat",
        notes="Corpus identity was discarded before release; future_count is the number of occurrences in the observation window",
    )
    expected = {3: 65_908, 4: 38_887, 5: 23_537}
    sequence = 0
    hit_totals: dict[int, int] = {}
    future_totals: dict[int, int] = {}
    with lag_writer:
        for frequency in (3, 4, 5):
            matrix = np.asarray(lag_data[f"results{frequency}"])
            if matrix.shape != (expected[frequency], frequency):
                raise AssertionError(f"Unexpected results{frequency} shape: {matrix.shape}")
            hit_totals[frequency] = int(np.count_nonzero(matrix[:, -1]))
            future_totals[frequency] = int(np.sum(matrix[:, -1], dtype=np.int64))
            for row in matrix:
                sequence += 1
                gaps = [int(value) for value in row[: frequency - 1]]
                future_count = int(row[-1])
                lag_writer.writerow({
                    "history_id": f"lag_history_{sequence:06d}", "frequency": frequency,
                    "gap_1": gaps[0], "gap_2": gaps[1],
                    "gap_3": gaps[2] if frequency >= 4 else None,
                    "gap_4": gaps[3] if frequency >= 5 else None,
                    "range_exclusive": sum(gaps), "future_count": future_count,
                    "future_hit": future_count > 0, "corpus": "unknown_combined",
                })
    if hit_totals != {3: 28_103, 4: 20_065, 5: 13_651}:
        raise AssertionError(f"Unexpected range-history hit totals: {hit_totals}")
    if future_totals != {3: 52_175, 4: 41_081, 5: 30_834}:
        raise AssertionError(f"Unexpected range-history future-count totals: {future_totals}")
    records.append(lag_writer.record())

    range_data = load_classic(source / "allNs1000.mat")
    range_writer = TableWriter(
        output,
        "range_frequency_histories.csv",
        [
            "history_id", "corpus_id", "corpus_family", "snapshot", "source_id",
            "source_index", "range_inclusive", "frequency", "future_count", "future_hit",
        ],
        table_id="range_frequency_histories",
        source_files="allNs1000.mat",
        notes="source_id joins directly to corpus_sources.csv within corpus_id",
    )
    sequence = 0
    with range_writer:
        for variable, corpus in (("resultsTW", "twitter"), ("resultsRE", "reddit")):
            cells = matlab_cells(range_data[variable])
            for cell_index, cell in enumerate(cells, start=1):
                if corpus == "reddit":
                    snapshot = "apr23" if cell_index <= 501 else "may5"
                    source_index = cell_index if cell_index <= 501 else cell_index - 501
                    corpus_id = f"reddit_{snapshot}"
                else:
                    snapshot = "release"
                    source_index = cell_index
                    corpus_id = "twitter"
                source_id = f"{corpus_id}_source_{source_index:03d}"
                array = np.asarray(cell)
                if array.size == 0:
                    continue
                array = np.atleast_2d(array)
                if array.shape[1] != 3:
                    raise AssertionError(f"Unexpected {variable} cell shape {array.shape}")
                for values in array:
                    sequence += 1
                    future_count = int(values[2])
                    range_writer.writerow({
                        "history_id": f"range_history_{sequence:06d}", "corpus_id": corpus_id,
                        "corpus_family": corpus, "snapshot": snapshot, "source_id": source_id,
                        "source_index": source_index,
                        "range_inclusive": int(values[0]), "frequency": int(values[1]),
                        "future_count": future_count, "future_hit": future_count > 0,
                    })
    records.append(range_writer.record())
    return records


def export_appendix_c_data(source: Path, output: Path) -> list[WrittenTable]:
    """Export the released Appendix-C distribution, decay, and revival data."""

    data = load_classic(source / "FigureC1.mat")
    records: list[WrittenTable] = []

    probabilities = np.asarray(data["probs"], dtype=float)
    if probabilities.shape != (225, 3) or not np.allclose(probabilities.sum(axis=0), 1.0):
        raise AssertionError("Unexpected Appendix-C frequency distributions")
    frequency_writer = TableWriter(
        output,
        "appendix_c_frequency_distributions.csv",
        ["frequency", "series", "probability"],
        table_id="appendix_c_frequency_distributions",
        source_files="FigureC1.mat",
        notes="Cached Figure C1a distributions; referenced upstream construction files were not released",
    )
    series_names = ("environmental_data", "power_am_release_simulation", "fitted_negative_binomial")
    with frequency_writer:
        for series0, series in enumerate(series_names):
            for frequency in range(1, 226):
                frequency_writer.writerow({
                    "frequency": frequency, "series": series,
                    "probability": probabilities[frequency - 1, series0],
                })
    records.append(frequency_writer.record())

    decay_writer = TableWriter(
        output,
        "appendix_c_decay_replicates.csv",
        ["series", "replicate_index", "statistic", "delay_events", "count"],
        table_id="appendix_c_decay_replicates",
        source_files="FigureC1.mat",
        notes="Column 1 is the normalization count; columns 2-201 are future counts at event delays 1-200",
    )
    decay_specs = (
        ("environmental_data", "decayData", 1502),
        ("exponential_am_release_simulation", "decayExp", 5000),
        ("power_am_release_simulation", "decayPower", 5000),
    )
    with decay_writer:
        for series, variable, expected_rows in decay_specs:
            matrix = np.asarray(data[variable])
            if matrix.shape != (expected_rows, 201):
                raise AssertionError(f"Unexpected {variable} shape: {matrix.shape}")
            for replicate0, row in enumerate(matrix):
                decay_writer.writerow({
                    "series": series, "replicate_index": replicate0 + 1,
                    "statistic": "baseline_normalization_count", "delay_events": None,
                    "count": int(row[0]),
                })
                for delay, count in enumerate(row[1:], start=1):
                    decay_writer.writerow({
                        "series": series, "replicate_index": replicate0 + 1,
                        "statistic": "future_count", "delay_events": delay,
                        "count": int(count),
                    })
    if decay_writer.rows != 2_311_902:
        # 1,502 + 5,000 + 5,000 replicates, each with 201 columns.
        raise AssertionError(f"Unexpected Appendix-C decay row count: {decay_writer.rows}")
    records.append(decay_writer.record())

    revival_writer = TableWriter(
        output,
        "appendix_c_revival_pairs.csv",
        [
            "series", "corpus_id", "source_cell_index", "within_cell_index",
            "extra_occurrences_period_1", "extra_occurrences_period_2",
        ],
        table_id="appendix_c_revival_pairs",
        source_files="FigureC1.mat",
        notes="Pairs used for the environmental-data and old power-A&M curves in Figure C1c",
    )
    revival_specs = (
        ("environmental_data", "twitter", "revivalsTwitter", 500, 4_166_793),
        ("environmental_data", "reddit_apr23", "revivalsApr23", 501, 1_550_224),
        ("environmental_data", "reddit_may5", "revivalsMay5", 501, 1_375_454),
        ("power_am_old_release_simulation", "simulation", "revivalsModel", 10_000, 2_272_362),
    )
    revival_rows = 0
    with revival_writer:
        for series, corpus_id, variable, expected_cells, expected_rows in revival_specs:
            cells = matlab_cells(data[variable])
            if len(cells) != expected_cells:
                raise AssertionError(f"Unexpected {variable} cell count: {len(cells)}")
            variable_rows = 0
            for cell0, cell in enumerate(cells):
                matrix = np.asarray(cell)
                if matrix.size == 0:
                    continue
                matrix = np.atleast_2d(matrix)
                if matrix.shape[1] != 2:
                    raise AssertionError(f"Unexpected {variable} cell shape: {matrix.shape}")
                for within0, row in enumerate(matrix):
                    revival_writer.writerow({
                        "series": series, "corpus_id": corpus_id,
                        "source_cell_index": cell0 + 1, "within_cell_index": within0 + 1,
                        "extra_occurrences_period_1": int(row[0]),
                        "extra_occurrences_period_2": int(row[1]),
                    })
                variable_rows += matrix.shape[0]
            if variable_rows != expected_rows:
                raise AssertionError(f"Unexpected {variable} row count: {variable_rows}")
            revival_rows += variable_rows
    if revival_rows != 9_364_833:
        raise AssertionError(f"Unexpected Appendix-C revival-pair rows: {revival_rows}")
    records.append(revival_writer.record())

    no_revival_writer = TableWriter(
        output,
        "appendix_c_no_revival_trials.csv",
        [
            "simulation_item", "within_item_index", "extra_occurrences_period_1",
            "extra_occurrences_period_2", "released_column_3", "released_column_4",
            "released_column_5", "released_column_6", "released_column_7", "released_column_8",
        ],
        table_id="appendix_c_no_revival_trials",
        source_files="FigureC1.mat",
        notes="First two columns produce the new/no-revival curve in Figure C1c; meanings of columns 3-8 are undocumented and their generator is absent",
    )
    cells = matlab_cells(data["resultsNoRevival"])
    if len(cells) != 25_000:
        raise AssertionError(f"Unexpected resultsNoRevival cell count: {len(cells)}")
    with no_revival_writer:
        for item0, cell in enumerate(cells):
            matrix = np.asarray(cell)
            if matrix.size == 0:
                continue
            matrix = np.atleast_2d(matrix)
            if matrix.shape[1] != 8:
                raise AssertionError(f"Unexpected resultsNoRevival shape: {matrix.shape}")
            for within0, row in enumerate(matrix):
                no_revival_writer.writerow({
                    "simulation_item": item0 + 1, "within_item_index": within0 + 1,
                    "extra_occurrences_period_1": row[0], "extra_occurrences_period_2": row[1],
                    "released_column_3": row[2], "released_column_4": row[3],
                    "released_column_5": row[4], "released_column_6": row[5],
                    "released_column_7": row[6], "released_column_8": row[7],
                })
    if no_revival_writer.rows != 2_219_920:
        raise AssertionError(f"Unexpected resultsNoRevival row count: {no_revival_writer.rows}")
    records.append(no_revival_writer.record())
    return records


def experiment_group(index: int) -> str:
    return "one_day" if index <= 5 else "between_days" if index <= 11 else "mixed"


def as_age_vector(value: Any) -> np.ndarray:
    array = np.asarray(value)
    return array.reshape(1).astype(float) if array.ndim == 0 else array.ravel(order="F").astype(float)


def condition_coordinates(shape: tuple[int, ...], linear0: int) -> tuple[int, int]:
    if len(shape) == 1:
        return linear0 + 1, 1
    row, col = np.unravel_index(linear0, shape, order="F")
    return int(row) + 1, int(col) + 1


def build_behavior_plot_points(
    experiment: int,
    data: np.ndarray,
    x_axis: Any,
    labels: Any,
) -> Iterator[dict[str, Any]]:
    label_values = [str(x) for x in matlab_cells(labels)] if np.asarray(labels).size else []
    shape = data.shape if data.ndim > 1 else (data.size,)
    mapping: list[tuple[int, int, int]] = []  # source linear (1-based), display row, display series
    if experiment == 1:
        return
    if experiment == 5:  # Young special reshape
        mapping.extend((i, i, 2) for i in range(1, 19))
        mapping.extend((i, i - 18, 1) for i in range(19, 30))
    elif experiment == 7:  # Cepeda 2008 ragged four-series reshape
        mapping.extend((i, i, 1) for i in range(1, 7))
        mapping.extend((i, i - 6, 2) for i in range(7, 14))
        mapping.extend((i, i - 13, 3) for i in range(14, 20))
        mapping.extend((i, i - 19, 4) for i in range(20, 27))
    elif experiment == 11:  # Pavlik & Anderson extractSubset.m
        keep = np.array([
            [1,1,1,6,6,6,11,11],[16,16,16,21,21,21,26,26],[31,31,np.nan,37,37,np.nan,43,np.nan],
            [32,32,17,38,38,22,44,27],[49,np.nan,18,57,np.nan,23,np.nan,28],
            [50,33,19,58,39,24,45,29],[51,34,20,59,40,25,46,30],
            [52,35,np.nan,60,41,np.nan,47,np.nan],[np.nan,36,np.nan,np.nan,42,np.nan,48,np.nan],
            [53,np.nan,np.nan,61,np.nan,np.nan,np.nan,np.nan],
            [54,np.nan,np.nan,62,np.nan,np.nan,np.nan,np.nan],
            [55,np.nan,np.nan,63,np.nan,np.nan,np.nan,np.nan],
            [56,np.nan,np.nan,64,np.nan,np.nan,np.nan,np.nan],
        ])
        for row in range(keep.shape[0]):
            for col in range(keep.shape[1]):
                if np.isfinite(keep[row, col]):
                    mapping.append((int(keep[row, col]), row + 1, col + 1))
    else:
        for linear0 in range(data.size):
            row, col = condition_coordinates(shape, linear0)
            display_row = row if experiment != 6 or row <= 3 else row + 1
            mapping.append((linear0 + 1, display_row, col))

    x_array = np.asarray(x_axis)
    point = 0
    for condition_linear, display_row, series in mapping:
        point += 1
        if x_array.size == 0:
            display_x = None
        elif experiment == 6:
            # Appendix.m inserts a blank separator row and replaces xAxis by 1:11.
            display_x = display_row
        elif x_array.ndim == 2:
            display_x = x_array[display_row - 1, series - 1]
        else:
            display_x = matlab_vector(x_array)[display_row - 1]
        yield {
            "plot_point_id": f"exp{experiment:02d}_plot_{point:03d}",
            "experiment_id": f"exp{experiment:02d}",
            "condition_id": f"exp{experiment:02d}_condition_{condition_linear:03d}",
            "display_row": display_row, "display_series_index": series,
            "display_series_label": label_values[series - 1] if series <= len(label_values) else None,
            "display_x": display_x,
            "provenance": "Appendix.m displayExperiment/drawit mapping",
        }


def export_behavioral_data(source: Path, output: Path) -> list[WrittenTable]:
    data = load_classic(source / "experiments14.mat")
    names = matlab_cells(data["names14"])
    observations = matlab_cells(data["data14"])
    patterns = matlab_cells(data["patterns14"])
    gaps = matlab_cells(data["gaps14"])
    labels = matlab_cells(data["labels14"])
    axes = matlab_cells(data["xAxis14"])
    xlabels = matlab_cells(data["xLabels14"])
    starts = np.asarray(data["starts14"], dtype=float)
    records: list[WrittenTable] = []

    exp_writer = TableWriter(
        output,
        "behavioral_experiments.csv",
        [
            "experiment_id", "source_order", "released_name", "paper_group", "response_measure",
            "x_axis_label", "age_clock", "range_clock", "source_file", "notes",
        ],
        table_id="behavioral_experiments",
        source_files="experiments14.mat;Appendix.m",
    )
    with exp_writer:
        for i, name in enumerate(names, start=1):
            notes = []
            if i == 4:
                notes.append("the first matrix row (eight conditions) has no schedules and release fixes predictions to 1/3")
            if i == 11:
                notes.append("two schedule variants per observed condition are averaged with equal weight")
            exp_writer.writerow({
                "experiment_id": f"exp{i:02d}", "source_order": i, "released_name": str(name),
                "paper_group": experiment_group(i), "response_measure": "aggregate_probability",
                "x_axis_label": str(xlabels[i - 1]) if np.asarray(xlabels[i - 1]).size else None,
                "age_clock": "released_encoded_units", "range_clock": "released_gap_input",
                "source_file": "experiments14.mat", "notes": "; ".join(notes),
            })
    records.append(exp_writer.record())

    condition_writer = TableWriter(
        output,
        "behavioral_conditions.csv",
        [
            "condition_id", "experiment_id", "matlab_linear_index", "matlab_row", "matlab_column",
            "observed_probability", "fixed_prediction", "schedule_aggregation", "n_schedule_variants",
        ],
        table_id="behavioral_conditions",
        source_files="experiments14.mat;Appendix.m",
        notes="Aggregate condition means only; participant outcomes and denominators were not released",
    )
    variant_writer = TableWriter(
        output,
        "behavioral_schedule_variants.csv",
        [
            "schedule_id", "condition_id", "variant_index", "mixture_weight", "source_pattern_index",
            "range_input", "n_presentations", "encoding_provenance",
        ],
        table_id="behavioral_schedule_variants",
        source_files="experiments14.mat;Appendix.m",
    )
    presentation_writer = TableWriter(
        output,
        "behavioral_presentations.csv",
        ["schedule_id", "presentation_index", "age_events"],
        table_id="behavioral_presentations",
        source_files="experiments14.mat",
        notes="Ages preserve source order, zeros, repeats, and fractional values exactly",
    )
    plot_writer = TableWriter(
        output,
        "behavioral_plot_points.csv",
        [
            "plot_point_id", "experiment_id", "condition_id", "display_row",
            "display_series_index", "display_series_label", "display_x", "provenance",
        ],
        table_id="behavioral_plot_points",
        source_files="experiments14.mat;Appendix.m",
    )
    with condition_writer, variant_writer, presentation_writer, plot_writer:
        for exp0 in range(14):
            exp = exp0 + 1
            obs = np.asarray(observations[exp0], dtype=float)
            flat_obs = obs.ravel(order="F")
            shape = obs.shape if obs.ndim > 1 else (obs.size,)
            for linear0, value in enumerate(flat_obs):
                row, col = condition_coordinates(shape, linear0)
                fixed = 1 / 3 if exp == 4 and row == 1 else None
                variants = 0 if fixed is not None else 2 if exp == 11 else 1
                condition_writer.writerow({
                    "condition_id": f"exp{exp:02d}_condition_{linear0 + 1:03d}",
                    "experiment_id": f"exp{exp:02d}", "matlab_linear_index": linear0 + 1,
                    "matlab_row": row, "matlab_column": col, "observed_probability": value,
                    "fixed_prediction": fixed,
                    "schedule_aggregation": "fixed" if fixed is not None else "equal_weight_mean" if exp == 11 else "single",
                    "n_schedule_variants": variants,
                })

            pattern_cells = matlab_cells(patterns[exp0])
            gap_values = matlab_vector(gaps[exp0])
            if len(pattern_cells) != len(gap_values):
                raise AssertionError(f"Experiment {exp}: pattern/gap count mismatch")
            for pattern0, (pattern, gap) in enumerate(zip(pattern_cells, gap_values, strict=True)):
                if exp == 4:
                    pattern_shape = np.asarray(patterns[exp0], dtype=object).shape
                    prow, pcol = condition_coordinates(pattern_shape, pattern0)
                    data_linear = (prow + 1) + obs.shape[0] * (pcol - 1)
                    condition_linear = data_linear
                    variant = 1
                    weight = 1.0
                elif exp == 11:
                    condition_linear = pattern0 % 64 + 1
                    variant = pattern0 // 64 + 1
                    weight = 0.5
                else:
                    condition_linear = pattern0 + 1
                    variant = 1
                    weight = 1.0
                ages = as_age_vector(pattern)
                schedule_id = f"exp{exp:02d}_condition_{condition_linear:03d}_schedule_{variant}"
                variant_writer.writerow({
                    "schedule_id": schedule_id,
                    "condition_id": f"exp{exp:02d}_condition_{condition_linear:03d}",
                    "variant_index": variant, "mixture_weight": weight,
                    "source_pattern_index": pattern0 + 1, "range_input": gap,
                    "n_presentations": len(ages), "encoding_provenance": "patterns14/gaps14 verbatim",
                })
                for age_index, age in enumerate(ages, start=1):
                    presentation_writer.writerow({
                        "schedule_id": schedule_id, "presentation_index": age_index, "age_events": age,
                    })
            plot_writer.writerows(build_behavior_plot_points(exp, obs, axes[exp0], labels[exp0]))
    if condition_writer.rows != 353 or variant_writer.rows != 409 or presentation_writer.rows != 2010:
        raise AssertionError(
            "Behavioral row invariant failed: "
            f"conditions={condition_writer.rows}, variants={variant_writer.rows}, presentations={presentation_writer.rows}"
        )
    records.extend([condition_writer.record(), variant_writer.record(), presentation_writer.record(), plot_writer.record()])

    starts_writer = TableWriter(
        output,
        "behavioral_fit_starts.csv",
        ["experiment_id", "b", "t_prior", "range_prior", "threshold", "noise_scale", "parameter_status"],
        table_id="behavioral_fit_starts",
        source_files="experiments14.mat;Appendix.m",
        notes="These are optimizer starting values, not final fitted estimates",
    )
    with starts_writer:
        for exp, row in enumerate(starts, start=1):
            starts_writer.writerow({
                "experiment_id": f"exp{exp:02d}", "b": row[0], "t_prior": row[1],
                "range_prior": row[2], "threshold": row[3], "noise_scale": row[4],
                "parameter_status": "optimizer_start",
            })
    records.append(starts_writer.record())
    return records


def export_model_parameters(source: Path, output: Path) -> WrittenTable:
    data = load_classic(source / "modelParams.mat")
    mappings = {
        "paramsGPE": ("gpe", ["frequency_exponent", "decay", "alpha"]),
        "paramsACTR": ("actr", ["decay", "alpha"]),
        "paramsPA": ("pavlik_anderson", ["activation_sensitivity", "minimum_decay", "alpha"]),
        "paramsPPE": ("ppe", ["recency_weight_exponent", "frequency_exponent", "minimum_decay", "spacing_sensitivity", "alpha"]),
        "paramsMCM": ("mcm", ["time_constant_scale", "time_constant_ratio", "total_trace_weight", "trace_weight_ratio", "alpha"]),
        "paramsAMPE": ("ampe", ["alpha", "decay_scale", "time_prior", "range_prior"]),
        "paramsEXP": ("anderson_milson_exponential", ["gamma_shape", "gamma_scale", "mean_decay", "mean_revival_interval"]),
        "paramsPOWER": ("anderson_milson_power", ["gamma_shape", "gamma_scale", "mean_decay", "mean_revival_interval"]),
    }
    writer = TableWriter(
        output,
        "environmental_model_parameters.csv",
        ["model_id", "parameter_index", "parameter_name", "value", "status", "source_variable"],
        table_id="environmental_model_parameters",
        source_files="modelParams.mat",
        notes="Released cached-fit vectors named in MATLAB parameter order; A&M output multipliers are not stored in modelParams.mat",
    )
    with writer:
        for variable, (model, names) in mappings.items():
            values = matlab_vector(data[variable])
            if len(values) != len(names):
                raise AssertionError(f"Parameter mapping mismatch for {variable}")
            for index, (name, value) in enumerate(zip(names, values, strict=True), start=1):
                writer.writerow({
                    "model_id": model, "parameter_index": index, "parameter_name": name,
                    "value": value, "status": "released_cached_fit", "source_variable": variable,
                })
    return writer.record()


def export_cached_surfaces(source: Path, output: Path) -> WrittenTable:
    old = load_classic(source / "display46.mat")
    final = load_classic(source / "display461.mat")
    frequency = {
        "data": np.asarray(old["DF6"][0]), "gpe": np.asarray(old["DF6"][1]),
        "actr": np.asarray(old["DF6"][2]), "pavlik_anderson": np.asarray(old["DF6"][3]),
        "ppe": np.asarray(old["DF6"][4]), "mcm": np.asarray(old["DF6"][5]),
        "anderson_milson_exponential": np.asarray(final["DF4"][1]),
        "anderson_milson_power": np.asarray(final["DF4"][2]), "ampe": np.asarray(final["DF4"][3]),
    }
    spacing = {
        "data": np.asarray(old["DS6"][0]), "gpe": np.asarray(old["DS6"][1]),
        "actr": np.asarray(old["DS6"][2]), "pavlik_anderson": np.asarray(old["DS6"][3]),
        "ppe": np.asarray(old["DS6"][4]), "mcm": np.asarray(old["DS6"][5]),
        "anderson_milson_exponential": np.asarray(final["DS4"][1]),
        "anderson_milson_power": np.asarray(final["DS4"][2]), "ampe": np.asarray(final["DS4"][3]),
    }
    expected_metrics = {
        "gpe": (0.581341646, 0.886138199, 513),
        "actr": (0.824587107, 0.784265359, 513),
        "pavlik_anderson": (0.683574373, 0.847171009, 513),
        "ppe": (0.548635787, 0.898373344, 513),
        "mcm": (0.578253532, 0.888356495, 513),
        "anderson_milson_exponential": (0.573769198, 0.905926876, 513),
        "anderson_milson_power": (0.409184328, 0.943891527, 481),
        "ampe": (0.397780742, 0.946968777, 513),
    }
    observed_vector = np.concatenate((frequency["data"].ravel(), spacing["data"].ravel()))
    for model, (expected_rmse, expected_r2, expected_n) in expected_metrics.items():
        predicted_vector = np.concatenate((frequency[model].ravel(), spacing[model].ravel()))
        valid = np.isfinite(observed_vector) & np.isfinite(predicted_vector)
        observed_log = np.log(observed_vector[valid])
        predicted_log = np.log(predicted_vector[valid])
        rmse = float(np.sqrt(np.mean((observed_log - predicted_log) ** 2)))
        r2 = float(np.corrcoef(observed_log, predicted_log)[0, 1] ** 2)
        if int(valid.sum()) != expected_n or not np.isclose(rmse, expected_rmse, atol=5e-9) or not np.isclose(r2, expected_r2, atol=5e-9):
            raise AssertionError(
                f"Cached metric mismatch for {model}: rmse={rmse}, r2={r2}, n={valid.sum()}"
            )
    writer = TableWriter(
        output,
        "published_environmental_surfaces.csv",
        ["panel", "model_id", "recency_bin", "condition_index", "probability", "cache_source"],
        table_id="published_environmental_surfaces",
        source_files="display46.mat;display461.mat",
        notes="Uses display461 for final A&M/AMPE arrays and display46 for deterministic models; stale A&M arrays are excluded",
    )
    with writer:
        for panel, surfaces in (("frequency_recency", frequency), ("spacing_recency", spacing)):
            for model, matrix in surfaces.items():
                expected_cols = 15 if panel == "frequency_recency" else 6
                if matrix.shape != (32, expected_cols):
                    raise AssertionError(f"Unexpected cached surface shape for {panel}/{model}: {matrix.shape}")
                cache = "display461.mat" if model.startswith("anderson_milson") or model == "ampe" else "display46.mat"
                for col in range(expected_cols):
                    for row in range(32):
                        writer.writerow({
                            "panel": panel, "model_id": model, "recency_bin": row + 1,
                            "condition_index": col + 1, "probability": matrix[row, col],
                            "cache_source": cache,
                        })
    return writer.record()


def export_hick(source: Path, output: Path) -> list[WrittenTable]:
    data = load_classic(source / "Schneider.mat")
    observed = np.asarray(data["schneiderData"])
    odds = np.asarray(data["odds"], dtype=float)
    full = np.asarray(data["fullprobs"], dtype=float)
    params = matlab_vector(data["params6"])
    base = load_classic(source / "base32.mat")
    bounds, means = matlab_vector(base["bounds"]), matlab_vector(base["means"])
    records: list[WrittenTable] = []

    obs_writer = TableWriter(
        output,
        "hick_observations.csv",
        ["n_alternatives", "repetition_status", "rt_ms"],
        table_id="hick_observations", source_files="Schneider.mat",
    )
    with obs_writer:
        for col, status in enumerate(("repeated", "nonrepeated")):
            for row, alternatives in enumerate((2, 4, 6)):
                obs_writer.writerow({"n_alternatives": alternatives, "repetition_status": status, "rt_ms": observed[row, col]})
    records.append(obs_writer.record())

    odds_writer = TableWriter(
        output,
        "hick_microenvironment_odds.csv",
        ["n_alternatives", "repetition_status", "odds"],
        table_id="hick_microenvironment_odds", source_files="Schneider.mat",
    )
    with odds_writer:
        for col, status in enumerate(("repeated", "nonrepeated")):
            for row, alternatives in enumerate(range(2, 7)):
                odds_writer.writerow({"n_alternatives": alternatives, "repetition_status": status, "odds": odds[row, col]})
    records.append(odds_writer.record())

    surface_writer = TableWriter(
        output,
        "hick_frequency_recency_surface.csv",
        [
            "corpus", "recency_bin", "recency_lower", "recency_upper", "recency_midpoint",
            "frequency_bin", "frequency_lower", "frequency_upper", "frequency_midpoint",
            "probability", "history_count",
        ],
        table_id="hick_frequency_recency_surface",
        source_files="Schneider.mat;Combined.mat;base32.mat",
    )
    combined = load_classic(source / "Combined.mat")
    counts = np.asarray(combined["countsA"])
    with surface_writer:
        for fbin in range(32):
            for rbin in range(32):
                surface_writer.writerow({
                    "corpus": "combined", "recency_bin": rbin + 1,
                    "recency_lower": int(bounds[rbin]) + 1,
                    "recency_upper": int(bounds[rbin + 1]), "recency_midpoint": means[rbin],
                    "frequency_bin": fbin + 1, "frequency_lower": int(bounds[fbin]) + 1,
                    "frequency_upper": int(bounds[fbin + 1]), "frequency_midpoint": means[fbin],
                    "probability": full[rbin, fbin], "history_count": int(counts[rbin, fbin]),
                })
    records.append(surface_writer.record())

    fit_writer = TableWriter(
        output,
        "hick_micro_fit.csv",
        ["parameter_name", "value", "status"],
        table_id="hick_micro_fit", source_files="Schneider.mat",
        notes="Released local-optimizer result for the microenvironment fit",
    )
    with fit_writer:
        for name, value in zip(("intercept_ms", "scale_ms", "power"), params, strict=True):
            fit_writer.writerow({"parameter_name": name, "value": value, "status": "released_cached_fit"})
    records.append(fit_writer.record())
    return records


def export_twitter_corpus(source: Path, output: Path) -> list[WrittenTable]:
    data = load_classic(source / "twitterData.mat")
    matrices = matlab_cells(data["allTweets"])
    times = matlab_cells(data["twitterTimes"])
    names = matlab_cells(data["tweeters"])
    vocabulary = matlab_cells(data["twitter20000"])
    counts = matlab_vector(data["countTweets"])
    if not (len(matrices) == len(times) == len(names) == 500):
        raise AssertionError("Unexpected Twitter source count")
    records: list[WrittenTable] = []

    source_writer = TableWriter(
        output,
        "corpus_sources.csv",
        ["corpus_id", "source_id", "source_index", "source_name", "text_count", "eligible_1001"],
        table_id="corpus_sources", source_files="twitterData.mat;redditData.mat",
        notes="Reddit rows are appended by the v7.3 exporter",
    )
    vocab_writer = TableWriter(
        output,
        "corpus_vocabulary.csv",
        ["corpus_id", "word_id", "token", "retained_occurrence_count"],
        table_id="corpus_vocabulary", source_files="twitterData.mat;redditData.mat",
        notes="Word IDs are scoped to corpus_id",
    )
    text_writer = TableWriter(
        output,
        "corpus_texts.csv",
        [
            "corpus_id", "source_id", "source_text_index", "corpus_text_index",
            "matlab_datenum", "timestamp_utc", "unique_token_count",
        ],
        table_id="corpus_texts", source_files="twitterData.mat;redditData.mat",
    )
    occurrence_writer = TableWriter(
        output,
        "corpus_occurrences.csv",
        [
            "corpus_id", "source_id", "source_text_index", "corpus_text_index",
            "position_in_text", "word_id",
        ],
        table_id="corpus_occurrences", source_files="twitterData.mat;redditData.mat",
        notes="Within-text duplicate word IDs were already collapsed by the authors",
    )
    corpus_text_index = 0
    occurrence_count = 0
    with source_writer, vocab_writer, text_writer, occurrence_writer:
        for word_id, (token, count) in enumerate(zip(vocabulary, counts, strict=True), start=1):
            vocab_writer.writerow({
                "corpus_id": "twitter", "word_id": word_id, "token": str(token),
                "retained_occurrence_count": int(count),
            })
        for source0, (matrix, time_values, name) in enumerate(zip(matrices, times, names, strict=True)):
            source_index = source0 + 1
            source_id = f"twitter_source_{source_index:03d}"
            matrix = np.atleast_2d(np.asarray(matrix))
            time_values = matlab_vector(time_values)
            if matrix.shape[0] != len(time_values):
                raise AssertionError(f"Twitter source {source_index} matrix/time length mismatch")
            source_writer.writerow({
                "corpus_id": "twitter", "source_id": source_id, "source_index": source_index,
                "source_name": str(name), "text_count": matrix.shape[0], "eligible_1001": matrix.shape[0] >= 1001,
            })
            for text0, (row, datenum) in enumerate(zip(matrix, time_values, strict=True)):
                corpus_text_index += 1
                words = np.asarray(row)
                words = words[words != 0].astype(np.int64, copy=False)
                if len(words) != len(np.unique(words)):
                    raise AssertionError(f"Twitter source {source_index}, text {text0 + 1} contains duplicate IDs")
                text_writer.writerow({
                    "corpus_id": "twitter", "source_id": source_id, "source_text_index": text0 + 1,
                    "corpus_text_index": corpus_text_index, "matlab_datenum": datenum,
                    "timestamp_utc": matlab_datenum_to_iso(float(datenum)), "unique_token_count": len(words),
                })
                for position, word_id in enumerate(words, start=1):
                    occurrence_count += 1
                    occurrence_writer.writerow({
                        "corpus_id": "twitter", "source_id": source_id,
                        "source_text_index": text0 + 1, "corpus_text_index": corpus_text_index,
                        "position_in_text": position, "word_id": int(word_id),
                    })
        if corpus_text_index != 1_029_655 or occurrence_count != 7_812_969:
            raise AssertionError(
                f"Twitter invariant failed: texts={corpus_text_index}, occurrences={occurrence_count}"
            )

        # Append Reddit in the same open relational tables.  The helper is a
        # read-only, dependency-free parser for the limited MATLAB-v7.3/HDF5
        # constructs used in this release.
        try:
            from matlab_v73_reader import append_reddit_csv_rows
        except ImportError as error:
            raise RuntimeError(
                "scripts/matlab_v73_reader.py is required for redditData.mat"
            ) from error
        reddit_counts = append_reddit_csv_rows(
            source / "redditData.mat",
            source_writer=source_writer,
            vocabulary_writer=vocab_writer,
            text_writer=text_writer,
            occurrence_writer=occurrence_writer,
            starting_corpus_text_index=corpus_text_index,
        )
        if reddit_counts["texts"] != 1_133_182:
            raise AssertionError(f"Reddit text count mismatch: {reddit_counts}")
    records.extend([
        source_writer.record(), vocab_writer.record(), text_writer.record(), occurrence_writer.record()
    ])
    return records


def write_manifest(output: Path, tables: Sequence[WrittenTable]) -> WrittenTable:
    writer = TableWriter(
        output,
        "manifest.csv",
        MANIFEST_FIELDS,
        table_id="manifest", source_files="all listed source files",
        notes="The manifest does not list itself so that generation is deterministic",
    )
    with writer:
        for table in sorted(tables, key=lambda item: item.filename):
            writer.writerow({
                "schema_version": SCHEMA_VERSION, "table_id": table.table_id,
                "filename": table.filename, "rows": table.rows, "bytes": table.bytes,
                "sha256": table.sha256, "source_files": table.source_files, "notes": table.notes,
            })
    return writer.record()


def validate_manifest(output: Path) -> None:
    manifest = output / "manifest.csv"
    listed: set[str] = set()
    table_ids: set[str] = set()
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != MANIFEST_FIELDS:
            raise AssertionError(
                f"Manifest header mismatch: {reader.fieldnames} != {MANIFEST_FIELDS}"
            )
        for row in reader:
            if row["schema_version"] != SCHEMA_VERSION:
                raise AssertionError(
                    f"Unexpected schema version for {row['filename']}: "
                    f"{row['schema_version']} != {SCHEMA_VERSION}"
                )
            filename = row["filename"]
            table_id = row["table_id"]
            if filename in listed:
                raise AssertionError(f"Duplicate manifest filename: {filename}")
            if table_id in table_ids:
                raise AssertionError(f"Duplicate manifest table_id: {table_id}")
            if Path(filename).name != filename or not filename.endswith(".csv"):
                raise AssertionError(f"Unsafe or invalid manifest filename: {filename}")
            listed.add(filename)
            table_ids.add(table_id)
            path = output / filename
            if not path.is_file():
                raise AssertionError(f"Manifest file missing: {path}")
            observed_bytes, observed_rows, observed_sha256 = csv_file_stats(path)
            if int(row["bytes"]) != observed_bytes:
                raise AssertionError(f"Byte count mismatch: {path}")
            if int(row["rows"]) != observed_rows:
                raise AssertionError(f"Row count mismatch: {path}")
            if row["sha256"] != observed_sha256:
                raise AssertionError(f"Checksum mismatch: {path}")
    actual = {path.name for path in output.glob("*.csv")}
    expected = listed | {"manifest.csv"}
    if actual != expected:
        raise AssertionError(
            "Output CSV inventory differs from manifest; "
            f"unlisted={sorted(actual - expected)}, missing={sorted(expected - actual)}"
        )
    allowed_regular_files = expected | {"README.md"}
    unexpected_regular_files = sorted(
        path.name
        for path in output.iterdir()
        if path.is_file() and path.name not in allowed_regular_files
    )
    if unexpected_regular_files:
        raise AssertionError(
            "Output directory contains unexpected regular files, including possible "
            f"interrupted-write artifacts: {unexpected_regular_files}"
        )


def ensure_skip_destination_is_safe(
    output: Path, *, skip_corpus: bool, skip_exact_history: bool, skip_appendix_c: bool
) -> None:
    """Prevent a partial run from silently mixing with a full prior export."""

    deliberately_excluded: set[str] = set()
    if skip_corpus:
        deliberately_excluded.update(CORPUS_CSVS)
    if skip_exact_history:
        deliberately_excluded.add("exact_history_cells.csv")
    if skip_appendix_c:
        deliberately_excluded.update(APPENDIX_C_CSVS)
    conflicts = sorted(name for name in deliberately_excluded if (output / name).exists())
    if conflicts:
        raise FileExistsError(
            "Skip-mode destination contains tables that this run would exclude. "
            "Use a fresh/separate output directory instead of creating a mixed snapshot: "
            + ", ".join(conflicts)
        )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir", type=Path, default=Path("external/anderson_2023/matlab"),
        help="Directory containing the immutable author MATLAB release",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/derived/anderson_2023"),
        help="Destination for generated CSV files",
    )
    parser.add_argument("--skip-corpus", action="store_true", help="Exclude the four large raw-corpus tables; requires a destination without them")
    parser.add_argument("--skip-exact-history", action="store_true", help="Exclude the 10.4-million-row exact-history table; requires a destination without it")
    parser.add_argument("--skip-appendix-c", action="store_true", help="Exclude the four large Appendix-C diagnostic tables; requires a destination without them")
    parser.add_argument("--validate-only", action="store_true", help="Validate an existing output manifest")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    if args.validate_only:
        validate_manifest(output)
        print(f"Validated {output / 'manifest.csv'}")
        return 0
    required = [
        "base32.mat", "Combined.mat", "twitterData.mat", "redditData.mat", "experiments14.mat",
        "results3to5.mat", "allNs1000.mat", "threesTwitter.mat", "threesReddit.mat",
        "display46.mat", "display461.mat", "modelParams.mat", "Schneider.mat", "FigureC1.mat",
    ]
    missing = [name for name in required if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing author-release files: {', '.join(missing)}")
    verify_source_checksums(source, required)
    output.mkdir(parents=True, exist_ok=True)
    ensure_skip_destination_is_safe(
        output,
        skip_corpus=args.skip_corpus,
        skip_exact_history=args.skip_exact_history,
        skip_appendix_c=args.skip_appendix_c,
    )
    tables: list[WrittenTable] = []
    tables.append(export_bin_definitions(source, output))
    tables.extend(export_environment_bins(source, output))
    if not args.skip_exact_history:
        tables.append(export_exact_history_cells(source, output))
    tables.extend(export_three_occurrence(source, output))
    tables.extend(export_range_histories(source, output))
    if not args.skip_appendix_c:
        tables.extend(export_appendix_c_data(source, output))
    tables.extend(export_behavioral_data(source, output))
    tables.append(export_model_parameters(source, output))
    tables.append(export_cached_surfaces(source, output))
    tables.extend(export_hick(source, output))
    if not args.skip_corpus:
        tables.extend(export_twitter_corpus(source, output))
    generated_filenames = [table.filename for table in tables]
    if len(generated_filenames) != len(set(generated_filenames)):
        raise AssertionError("The export generated a duplicate table filename")
    generated_table_ids = [table.table_id for table in tables]
    if len(generated_table_ids) != len(set(generated_table_ids)):
        raise AssertionError("The export generated a duplicate table_id")
    expected_filenames = set(KNOWN_GENERATED_CSVS)
    if args.skip_corpus:
        expected_filenames.difference_update(CORPUS_CSVS)
    if args.skip_exact_history:
        expected_filenames.remove("exact_history_cells.csv")
    if args.skip_appendix_c:
        expected_filenames.difference_update(APPENDIX_C_CSVS)
    observed_filenames = set(generated_filenames)
    if observed_filenames != expected_filenames:
        raise AssertionError(
            "Generated table inventory differs from the registered mode; "
            f"unexpected={sorted(observed_filenames - expected_filenames)}, "
            f"missing={sorted(expected_filenames - observed_filenames)}"
        )
    write_manifest(output, tables)
    validate_manifest(output)
    print(f"Wrote and validated {len(tables)} data tables in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

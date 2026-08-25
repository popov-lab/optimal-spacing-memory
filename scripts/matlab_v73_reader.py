"""Small, dependency-light reader for the release's MATLAB 7.3 file.

MATLAB 7.3 MAT files are HDF5 containers.  SciPy deliberately does not read
them, and the execution environment used for this project does not always
provide ``h5py``.  This module implements only the old-style HDF5 structures
used by ``external/anderson_2023/matlab/redditData.mat``:

* superblock version 0 with an old symbol-table root group;
* object-header version 1;
* compact, contiguous, or version-3 chunked datasets;
* the legacy raw-data chunk B-tree; and
* unfiltered or DEFLATE-compressed chunks.

It is intentionally not a general HDF5 or MAT-file implementation.  The
reader is useful for a reproducible conversion pipeline because it has no
binary dependency beyond NumPy, preserves the source dimensions and object
references, and fails explicitly on unsupported structures.

Examples
--------
List top-level variables::

    python scripts/matlab_v73_reader.py \
        external/anderson_2023/matlab/redditData.mat

Read the first April subreddit as its MATLAB ``n x 100`` array::

    with MatlabV73File(path) as mat:
        comments = mat.read_cell("arraysApr23", 0, matlab_order=True)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import struct
from typing import BinaryIO, Iterator
import zlib

import numpy as np


HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
UNDEFINED_ADDRESS = (1 << 64) - 1


class UnsupportedHDF5Error(RuntimeError):
    """Raised when a file uses an HDF5 feature outside this small reader."""


@dataclass(frozen=True)
class DatasetInfo:
    """Structural metadata for one HDF5 dataset.

    ``hdf_shape`` is the dimension order stored by MATLAB's HDF5 writer.
    MATLAB array dimensions are reversed in this file; ``matlab_shape`` is
    therefore the reversed tuple.
    """

    reference: int
    hdf_shape: tuple[int, ...]
    matlab_shape: tuple[int, ...]
    numpy_dtype: np.dtype
    layout: str
    chunk_shape: tuple[int, ...] | None = None


@dataclass(frozen=True)
class RedditPartitionMetadata:
    """Metadata for one independently encoded Reddit day."""

    date: str
    arrays_variable: str
    counts_variable: str
    subreddits_variable: str
    vocabulary_variable: str
    source_names: tuple[str, ...]
    event_counts: tuple[int, ...]
    event_width: int
    vocabulary_size: int
    retained_item_count: int

    @property
    def source_count(self) -> int:
        return len(self.source_names)

    @property
    def event_count(self) -> int:
        return sum(self.event_counts)

    @property
    def eligible_source_count(self) -> int:
        """Sources long enough to define at least one 1,001-event window."""

        return sum(count >= 1_001 for count in self.event_counts)


@dataclass(frozen=True)
class _Message:
    type: int
    body: bytes
    flags: int


class MatlabV73File:
    """Read the subset of HDF5 used by Anderson et al.'s Reddit MAT file."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._file: BinaryIO = self.path.open("rb")
        self._signature_offset = self._find_signature()
        self._base_address, self._root_btree, self._root_heap = (
            self._read_superblock()
        )
        self._variables = self._read_root_variables()

    def __enter__(self) -> "MatlabV73File":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def close(self) -> None:
        if not self._file.closed:
            self._file.close()

    @property
    def variable_names(self) -> tuple[str, ...]:
        """Top-level MATLAB variable names, excluding MATLAB's refs group."""

        return tuple(sorted(self._variables))

    def reference(self, name: str) -> int:
        """Return the HDF5 object reference for a top-level variable."""

        try:
            return self._variables[name]
        except KeyError as error:
            raise KeyError(f"No top-level variable named {name!r}") from error

    def info(self, target: str | int) -> DatasetInfo:
        """Return dataset metadata without loading its values."""

        reference = self.reference(target) if isinstance(target, str) else target
        metadata = self._dataset_metadata(reference)
        shape = metadata["shape"]
        dtype = metadata["dtype"]
        layout_body = metadata["layout"]
        layout_class = layout_body[1]
        if layout_class == 0:
            layout = "compact"
            chunk_shape = None
        elif layout_class == 1:
            layout = "contiguous"
            chunk_shape = None
        elif layout_class == 2:
            layout = "chunked"
            dimensions = layout_body[2]
            # The final stored dimension is the element byte size, not a
            # logical array dimension.
            chunk_shape = tuple(
                self._u32(layout_body, 11 + 4 * i)
                for i in range(dimensions - 1)
            )
        else:
            raise UnsupportedHDF5Error(
                f"Unsupported data-layout class {layout_class}"
            )
        return DatasetInfo(
            reference=reference,
            hdf_shape=shape,
            matlab_shape=tuple(reversed(shape)),
            numpy_dtype=dtype,
            layout=layout,
            chunk_shape=chunk_shape,
        )

    def read_dataset(
        self, target: str | int, *, matlab_order: bool = False
    ) -> np.ndarray:
        """Read one numeric, character, or object-reference dataset.

        By default, the returned array uses the dimensions physically stored
        in HDF5.  Set ``matlab_order=True`` to reverse the axes into the shape
        seen by MATLAB (for example, ``100 x n`` becomes ``n x 100``).
        """

        reference = self.reference(target) if isinstance(target, str) else target
        metadata = self._dataset_metadata(reference)
        shape = metadata["shape"]
        dtype = metadata["dtype"]
        layout = metadata["layout"]
        layout_class = layout[1]

        if layout[0] != 3:
            raise UnsupportedHDF5Error(
                f"Only data-layout message version 3 is supported, got {layout[0]}"
            )

        if layout_class == 0:
            stored_size = self._u16(layout, 2)
            raw = layout[4 : 4 + stored_size]
            array = np.frombuffer(raw, dtype=dtype).reshape(shape).copy()
        elif layout_class == 1:
            address = self._u64(layout, 2)
            stored_size = self._u64(layout, 10)
            raw = self._read_at(self._physical(address), stored_size)
            array = np.frombuffer(raw, dtype=dtype).reshape(shape).copy()
        elif layout_class == 2:
            array = self._read_chunked(shape, dtype, layout)
        else:
            raise UnsupportedHDF5Error(
                f"Unsupported data-layout class {layout_class}"
            )

        if matlab_order and array.ndim > 1:
            array = array.transpose(tuple(reversed(range(array.ndim))))
        return array

    def references(self, name: str) -> np.ndarray:
        """Read a cell-array variable as a flat vector of object references."""

        array = self.read_dataset(name)
        if array.dtype != np.dtype("<u8"):
            raise TypeError(f"{name!r} is not an HDF5 object-reference dataset")
        return array.ravel()

    def iter_cells(
        self, name: str, *, matlab_order: bool = True
    ) -> Iterator[np.ndarray]:
        """Yield cell contents without retaining the whole cell array in RAM."""

        for reference in self.references(name):
            if reference == 0:
                yield np.empty((0, 0))
            else:
                yield self.read_dataset(
                    int(reference), matlab_order=matlab_order
                )

    def read_cell(
        self, name: str, index: int, *, matlab_order: bool = True
    ) -> np.ndarray:
        """Read one zero-based cell from a top-level MATLAB cell array."""

        references = self.references(name)
        reference = int(references[index])
        if reference == 0:
            return np.empty((0, 0))
        return self.read_dataset(reference, matlab_order=matlab_order)

    def read_char(self, target: int) -> str:
        """Decode one referenced MATLAB char array."""

        array = self.read_dataset(target)
        if array.dtype.itemsize != 2:
            raise TypeError(
                f"Expected a two-byte MATLAB char array, got {array.dtype}"
            )
        return array.astype("<u2", copy=False).tobytes().decode("utf-16le")

    def read_string_cells(self, name: str) -> list[str]:
        """Decode a top-level cell array of MATLAB character arrays."""

        values: list[str] = []
        for reference in self.references(name):
            values.append("" if reference == 0 else self.read_char(int(reference)))
        return values

    def _find_signature(self) -> int:
        offset = 0
        while offset <= 1 << 20:
            if self._read_at(offset, 8) == HDF5_SIGNATURE:
                return offset
            offset = 512 if offset == 0 else offset * 2
        raise UnsupportedHDF5Error("HDF5 signature not found in the first MiB")

    def _read_superblock(self) -> tuple[int, int, int]:
        block = self._read_at(self._signature_offset, 96)
        if block[8] != 0:
            raise UnsupportedHDF5Error(
                f"Only HDF5 superblock version 0 is supported, got {block[8]}"
            )
        if block[13] != 8 or block[14] != 8:
            raise UnsupportedHDF5Error(
                "Only eight-byte HDF5 addresses and lengths are supported"
            )
        base_address = self._u64(block, 24)
        root_entry = 56
        cache_type = self._u32(block, root_entry + 16)
        if cache_type != 1:
            raise UnsupportedHDF5Error(
                "Expected an old-style symbol-table root group"
            )
        root_btree = self._u64(block, root_entry + 24)
        root_heap = self._u64(block, root_entry + 32)
        return base_address, root_btree, root_heap

    def _read_root_variables(self) -> dict[str, int]:
        heap_header = self._read_at(self._physical(self._root_heap), 32)
        if heap_header[:4] != b"HEAP":
            raise UnsupportedHDF5Error("Invalid local-heap signature")
        heap_data = self._u64(heap_header, 24)
        variables: dict[str, int] = {}
        for symbol_node in self._group_btree_nodes(self._root_btree):
            header = self._read_at(self._physical(symbol_node), 8)
            if header[:4] != b"SNOD":
                raise UnsupportedHDF5Error("Invalid symbol-table node signature")
            entries = self._u16(header, 6)
            payload = self._read_at(self._physical(symbol_node) + 8, 40 * entries)
            for i in range(entries):
                entry = payload[40 * i : 40 * (i + 1)]
                name_offset = self._u64(entry, 0)
                object_reference = self._u64(entry, 8)
                name = self._read_c_string(self._physical(heap_data) + name_offset)
                if not name.startswith("#"):
                    variables[name] = object_reference
        return variables

    def _group_btree_nodes(self, reference: int) -> list[int]:
        header = self._read_at(self._physical(reference), 24)
        if header[:4] != b"TREE" or header[4] != 0:
            raise UnsupportedHDF5Error("Invalid group B-tree node")
        level = header[5]
        entries = self._u16(header, 6)
        payload = self._read_at(self._physical(reference) + 24, 16 * entries + 8)
        children = [self._u64(payload, 8 + 16 * i) for i in range(entries)]
        if level == 0:
            return children
        nodes: list[int] = []
        for child in children:
            nodes.extend(self._group_btree_nodes(child))
        return nodes

    def _object_messages(self, reference: int) -> list[_Message]:
        address = self._physical(reference)
        header = self._read_at(address, 16)
        if not header or header[0] != 1:
            version = None if not header else header[0]
            raise UnsupportedHDF5Error(
                f"Only object-header version 1 is supported, got {version} "
                f"at reference 0x{reference:x}"
            )
        message_count = self._u16(header, 2)
        first_chunk_size = self._u32(header, 8)
        chunks = [self._read_at(address + 16, first_chunk_size)]
        messages: list[_Message] = []
        chunk_index = 0
        while chunk_index < len(chunks) and len(messages) < message_count:
            chunk = chunks[chunk_index]
            offset = 0
            while offset + 8 <= len(chunk) and len(messages) < message_count:
                message_type = self._u16(chunk, offset)
                size = self._u16(chunk, offset + 2)
                flags = chunk[offset + 4]
                body_start = offset + 8
                body_end = body_start + size
                if body_end > len(chunk):
                    raise UnsupportedHDF5Error("Object-header message is truncated")
                body = chunk[body_start:body_end]
                messages.append(_Message(message_type, body, flags))
                if message_type == 16:
                    continuation = self._u64(body, 0)
                    continuation_size = self._u64(body, 8)
                    chunks.append(
                        self._read_at(
                            self._physical(continuation), continuation_size
                        )
                    )
                offset += 8 + self._align8(size)
            chunk_index += 1
        if len(messages) != message_count:
            raise UnsupportedHDF5Error(
                f"Expected {message_count} object messages, read {len(messages)}"
            )
        return messages

    def _dataset_metadata(self, reference: int) -> dict[str, object]:
        shape: tuple[int, ...] | None = None
        dtype: np.dtype | None = None
        layout: bytes | None = None
        for message in self._object_messages(reference):
            if message.type == 1:
                rank = message.body[1]
                shape = tuple(
                    self._u64(message.body, 8 + 8 * i) for i in range(rank)
                )
            elif message.type == 3:
                dtype = self._numpy_dtype(message.body)
            elif message.type == 8:
                layout = message.body
        if shape is None or dtype is None or layout is None:
            raise UnsupportedHDF5Error(
                f"Object 0x{reference:x} is not a supported dataset"
            )
        return {"shape": shape, "dtype": dtype, "layout": layout}

    def _numpy_dtype(self, body: bytes) -> np.dtype:
        datatype_class = body[0] & 0x0F
        element_size = self._u32(body, 4)
        if datatype_class == 0:  # integer; MATLAB chars are unsigned uint16
            flags = int.from_bytes(body[1:4], "little")
            signed = bool(flags & 0x08)
            code = "i" if signed else "u"
            return np.dtype(f"<{code}{element_size}")
        if datatype_class == 1:  # IEEE floating point
            if element_size not in (4, 8):
                raise UnsupportedHDF5Error(
                    f"Unsupported floating-point size {element_size}"
                )
            return np.dtype(f"<f{element_size}")
        if datatype_class == 7 and element_size == 8:  # object reference
            return np.dtype("<u8")
        raise UnsupportedHDF5Error(
            f"Unsupported datatype class {datatype_class}, size {element_size}"
        )

    def _read_chunked(
        self, shape: tuple[int, ...], dtype: np.dtype, layout: bytes
    ) -> np.ndarray:
        dimensions = layout[2]
        btree = self._u64(layout, 3)
        stored_chunk_shape = tuple(
            self._u32(layout, 11 + 4 * i) for i in range(dimensions)
        )
        if stored_chunk_shape[-1] != dtype.itemsize:
            raise UnsupportedHDF5Error(
                "Chunk layout's element-size dimension does not match datatype"
            )
        chunk_shape = stored_chunk_shape[:-1]
        if len(chunk_shape) != len(shape):
            raise UnsupportedHDF5Error("Chunk rank does not match dataset rank")

        result = np.zeros(shape, dtype=dtype)
        for stored_size, filter_mask, offsets, child in self._chunk_leaves(
            btree, dimensions
        ):
            raw = self._read_at(self._physical(child), stored_size)
            if filter_mask & ~1:
                raise UnsupportedHDF5Error(
                    f"Unsupported chunk filter mask 0x{filter_mask:x}"
                )
            if not (filter_mask & 1):
                raw = zlib.decompress(raw)
            expected_size = int(np.prod(chunk_shape)) * dtype.itemsize
            if len(raw) != expected_size:
                raise UnsupportedHDF5Error(
                    f"Chunk expands to {len(raw)} bytes; expected {expected_size}"
                )
            chunk = np.frombuffer(raw, dtype=dtype).reshape(chunk_shape)
            destination_slices: list[slice] = []
            source_slices: list[slice] = []
            for offset, chunk_length, array_length in zip(
                offsets[:-1], chunk_shape, shape
            ):
                copied = max(0, min(chunk_length, array_length - offset))
                destination_slices.append(slice(offset, offset + copied))
                source_slices.append(slice(0, copied))
            result[tuple(destination_slices)] = chunk[tuple(source_slices)]
        return result

    def _chunk_leaves(
        self, reference: int, dimensions: int
    ) -> list[tuple[int, int, tuple[int, ...], int]]:
        header = self._read_at(self._physical(reference), 24)
        if header[:4] != b"TREE" or header[4] != 1:
            raise UnsupportedHDF5Error("Invalid raw-data chunk B-tree node")
        level = header[5]
        entries = self._u16(header, 6)
        key_size = 8 + 8 * dimensions
        entry_size = key_size + 8
        payload = self._read_at(
            self._physical(reference) + 24, entry_size * entries + key_size
        )
        leaves: list[tuple[int, int, tuple[int, ...], int]] = []
        for i in range(entries):
            start = entry_size * i
            key = payload[start : start + key_size]
            child = self._u64(payload, start + key_size)
            if level == 0:
                offsets = tuple(
                    self._u64(key, 8 + 8 * j) for j in range(dimensions)
                )
                leaves.append(
                    (self._u32(key, 0), self._u32(key, 4), offsets, child)
                )
            else:
                leaves.extend(self._chunk_leaves(child, dimensions))
        return leaves

    def _physical(self, reference: int) -> int:
        if reference == UNDEFINED_ADDRESS:
            raise UnsupportedHDF5Error("Encountered an undefined HDF5 address")
        return self._base_address + reference

    def _read_at(self, offset: int, size: int) -> bytes:
        self._file.seek(offset)
        data = self._file.read(size)
        if len(data) != size:
            raise EOFError(f"Could not read {size} bytes at file offset {offset}")
        return data

    def _read_c_string(self, offset: int) -> str:
        self._file.seek(offset)
        chunks: list[bytes] = []
        while True:
            chunk = self._file.read(64)
            if not chunk:
                raise EOFError("Unterminated HDF5 heap string")
            null = chunk.find(b"\0")
            if null >= 0:
                chunks.append(chunk[:null])
                break
            chunks.append(chunk)
        return b"".join(chunks).decode("utf-8")

    @staticmethod
    def _align8(size: int) -> int:
        return (size + 7) // 8 * 8

    @staticmethod
    def _u16(data: bytes, offset: int = 0) -> int:
        return struct.unpack_from("<H", data, offset)[0]

    @staticmethod
    def _u32(data: bytes, offset: int = 0) -> int:
        return struct.unpack_from("<I", data, offset)[0]

    @staticmethod
    def _u64(data: bytes, offset: int = 0) -> int:
        return struct.unpack_from("<Q", data, offset)[0]


def iter_reddit_cells(
    path: str | Path,
    top_level_name: str,
    *,
    matlab_order: bool = True,
) -> Iterator[tuple[int, np.ndarray | str]]:
    """Stream a Reddit cell array as ``(zero_based_index, value)`` pairs.

    Numeric comment cells are yielded one at a time, so the full collection
    of roughly 17 million retained item IDs is never resident in memory.
    Two-byte integer cells (the vocabulary and subreddit variables) are
    decoded as MATLAB UTF-16 character arrays and yielded as Python strings.
    """

    with MatlabV73File(path) as mat:
        for index, reference in enumerate(mat.references(top_level_name)):
            if reference == 0:
                yield index, ""
                continue
            info = mat.info(int(reference))
            if info.numpy_dtype.kind == "u" and info.numpy_dtype.itemsize == 2:
                yield index, mat.read_char(int(reference))
            else:
                yield index, mat.read_dataset(
                    int(reference), matlab_order=matlab_order
                )


def read_reddit_metadata(
    path: str | Path,
) -> dict[str, RedditPartitionMetadata]:
    """Read source-level metadata without loading any comment matrices.

    The returned mapping is keyed by ISO date.  Reddit term IDs are local to
    these partitions because each day has its own vocabulary array.
    """

    specifications = {
        "2021-04-23": (
            "arraysApr23",
            "countsApr23",
            "subredditsApr23",
            "vocabApr23",
        ),
        "2021-05-05": (
            "arraysMay5",
            "countsMay5",
            "subredditsMay5",
            "vocabMay5",
        ),
    }
    metadata: dict[str, RedditPartitionMetadata] = {}
    with MatlabV73File(path) as mat:
        for date, names in specifications.items():
            arrays_name, counts_name, sources_name, vocabulary_name = names
            array_references = mat.references(arrays_name)
            cell_shapes = [
                mat.info(int(reference)).matlab_shape
                for reference in array_references
            ]
            widths = {shape[1] for shape in cell_shapes}
            if len(widths) != 1:
                raise ValueError(
                    f"{arrays_name} has inconsistent comment widths: {widths}"
                )
            source_names = tuple(mat.read_string_cells(sources_name))
            counts = mat.read_dataset(counts_name).ravel()
            vocabulary_size = mat.references(vocabulary_name).size
            metadata[date] = RedditPartitionMetadata(
                date=date,
                arrays_variable=arrays_name,
                counts_variable=counts_name,
                subreddits_variable=sources_name,
                vocabulary_variable=vocabulary_name,
                source_names=source_names,
                event_counts=tuple(shape[0] for shape in cell_shapes),
                event_width=widths.pop(),
                vocabulary_size=int(vocabulary_size),
                retained_item_count=int(counts.sum()),
            )
    return metadata


def append_reddit_csv_rows(
    path: str | Path,
    *,
    source_writer: object,
    vocabulary_writer: object,
    text_writer: object,
    occurrence_writer: object,
    starting_corpus_text_index: int,
) -> dict[str, int]:
    """Append both Reddit partitions to the pipeline's relational CSVs.

    The writer arguments need only expose a ``writerow(mapping)`` method, as
    :class:`csv.DictWriter` and the pipeline's atomic ``TableWriter`` do.  A
    single ``n x 100`` comment matrix is decoded at a time.  Event rows with
    no retained vocabulary items are deliberately retained because deleting
    them would change every later event-distance calculation.

    Reddit's two term dictionaries are independent.  They are therefore
    exported with distinct ``corpus_id`` values (``reddit_apr23`` and
    ``reddit_may5``); a numeric ``word_id`` is never meaningful without that
    partition key.
    """

    if starting_corpus_text_index < 0:
        raise ValueError("starting_corpus_text_index must be nonnegative")
    for writer_name, writer in (
        ("source_writer", source_writer),
        ("vocabulary_writer", vocabulary_writer),
        ("text_writer", text_writer),
        ("occurrence_writer", occurrence_writer),
    ):
        if not callable(getattr(writer, "writerow", None)):
            raise TypeError(f"{writer_name} must expose writerow(mapping)")

    specifications = (
        (
            "reddit_apr23",
            "arraysApr23",
            "countsApr23",
            "subredditsApr23",
            "vocabApr23",
            573_836,
            8_351_677,
        ),
        (
            "reddit_may5",
            "arraysMay5",
            "countsMay5",
            "subredditsMay5",
            "vocabMay5",
            559_346,
            8_177_928,
        ),
    )
    totals = {
        "sources": 0,
        "vocabulary_rows": 0,
        "texts": 0,
        "occurrences": 0,
        "empty_texts": 0,
        "at_storage_limit_texts": 0,
        "ending_corpus_text_index": starting_corpus_text_index,
    }
    corpus_text_index = starting_corpus_text_index

    with MatlabV73File(path) as mat:
        for (
            corpus_id,
            arrays_name,
            counts_name,
            sources_name,
            vocabulary_name,
            expected_texts,
            expected_occurrences,
        ) in specifications:
            source_names = mat.read_string_cells(sources_name)
            array_references = mat.references(arrays_name)
            vocabulary_references = mat.references(vocabulary_name)
            released_counts = (
                mat.read_dataset(counts_name).ravel().astype(np.int64)
            )
            if len(source_names) != len(array_references):
                raise AssertionError(
                    f"{corpus_id}: source-name/cell count mismatch"
                )
            if len(vocabulary_references) != len(released_counts):
                raise AssertionError(
                    f"{corpus_id}: vocabulary/count-vector length mismatch"
                )

            for word0, (reference, retained_count) in enumerate(
                zip(vocabulary_references, released_counts, strict=True)
            ):
                if reference == 0:
                    raise AssertionError(
                        f"{corpus_id}: null vocabulary cell at index {word0 + 1}"
                    )
                vocabulary_writer.writerow(
                    {
                        "corpus_id": corpus_id,
                        "word_id": word0 + 1,
                        "token": mat.read_char(int(reference)),
                        "retained_occurrence_count": int(retained_count),
                    }
                )
                totals["vocabulary_rows"] += 1

            observed_counts = np.zeros(20_001, dtype=np.int64)
            partition_texts = 0
            partition_occurrences = 0
            for source0, (source_name, reference) in enumerate(
                zip(source_names, array_references, strict=True)
            ):
                source_index = source0 + 1
                source_id = f"{corpus_id}_source_{source_index:03d}"
                comments = mat.read_dataset(
                    int(reference), matlab_order=True
                )
                if comments.ndim != 2 or comments.shape[1] != 100:
                    raise AssertionError(
                        f"{source_id}: expected n x 100 comments, got "
                        f"{comments.shape}"
                    )
                if not np.all(np.isfinite(comments)) or not np.all(
                    comments == np.floor(comments)
                ):
                    raise AssertionError(f"{source_id}: noninteger word ID")
                if np.any((comments < 0) | (comments > 20_000)):
                    raise AssertionError(f"{source_id}: word ID outside 0..20000")

                positive = comments > 0
                seen_zero = np.maximum.accumulate(comments == 0, axis=1)
                if np.any(seen_zero & positive):
                    raise AssertionError(
                        f"{source_id}: a positive ID occurs after zero padding"
                    )
                sorted_rows = np.sort(comments, axis=1)
                duplicate = (sorted_rows[:, 1:] == sorted_rows[:, :-1]) & (
                    sorted_rows[:, 1:] != 0
                )
                if np.any(duplicate):
                    raise AssertionError(
                        f"{source_id}: duplicate word ID within a comment"
                    )

                text_count = comments.shape[0]
                source_writer.writerow(
                    {
                        "corpus_id": corpus_id,
                        "source_id": source_id,
                        "source_index": source_index,
                        "source_name": source_name,
                        "text_count": text_count,
                        "eligible_1001": text_count >= 1_001,
                    }
                )
                totals["sources"] += 1

                values = comments[positive].astype(np.int64)
                observed_counts += np.bincount(values, minlength=20_001)
                for text0, row in enumerate(comments):
                    words = row[row > 0].astype(np.int64, copy=False)
                    corpus_text_index += 1
                    partition_texts += 1
                    totals["texts"] += 1
                    if len(words) == 0:
                        totals["empty_texts"] += 1
                    if len(words) == 100:
                        totals["at_storage_limit_texts"] += 1
                    text_writer.writerow(
                        {
                            "corpus_id": corpus_id,
                            "source_id": source_id,
                            "source_text_index": text0 + 1,
                            "corpus_text_index": corpus_text_index,
                            "matlab_datenum": None,
                            "timestamp_utc": None,
                            "unique_token_count": len(words),
                        }
                    )
                    for position, word_id in enumerate(words, start=1):
                        occurrence_writer.writerow(
                            {
                                "corpus_id": corpus_id,
                                "source_id": source_id,
                                "source_text_index": text0 + 1,
                                "corpus_text_index": corpus_text_index,
                                "position_in_text": position,
                                "word_id": int(word_id),
                            }
                        )
                    partition_occurrences += len(words)
                    totals["occurrences"] += len(words)

            if partition_texts != expected_texts:
                raise AssertionError(
                    f"{corpus_id}: texts={partition_texts}, "
                    f"expected {expected_texts}"
                )
            if partition_occurrences != expected_occurrences:
                raise AssertionError(
                    f"{corpus_id}: occurrences={partition_occurrences}, "
                    f"expected {expected_occurrences}"
                )
            if not np.array_equal(observed_counts[1:], released_counts):
                raise AssertionError(
                    f"{corpus_id}: streamed IDs do not reproduce released counts"
                )
            totals[f"{corpus_id}_texts"] = partition_texts
            totals[f"{corpus_id}_occurrences"] = partition_occurrences

    totals["ending_corpus_text_index"] = corpus_text_index
    return totals


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mat_file", type=Path)
    arguments = parser.parse_args()
    with MatlabV73File(arguments.mat_file) as mat:
        for name in mat.variable_names:
            info = mat.info(name)
            print(
                f"{name:<18} MATLAB {info.matlab_shape!s:<15} "
                f"HDF5 {info.hdf_shape!s:<15} {info.numpy_dtype} "
                f"{info.layout}"
            )


if __name__ == "__main__":
    _main()

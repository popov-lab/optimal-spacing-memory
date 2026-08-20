# Repository guidance for coding agents

Keep changes consistent with the scientific and data-provenance constraints
below. These are the non-obvious repository facts that are expensive or risky
to rediscover.

## Data and provenance

- Treat `external/anderson_2023/matlab/` as an immutable vendor snapshot. Do
  not fix formatting, whitespace, or source defects in place. Put adapters,
  corrected implementations, and commentary in `scripts/`, `src/`, or `docs/`.
- Verify vendor integrity from `external/anderson_2023/matlab/` with
  `shasum -a 256 -c ../SHA256SUMS`.
- Sixteen `.mat` and `.docx` files use Git LFS. Before investigating missing or
  malformed vendor data, run `git lfs pull`, `git lfs checkout`, and
  `git lfs fsck`. A pointer file is not the source artifact.
- Generated CSVs under `data/derived/anderson_2023/` are deterministic and
  ignored except for `manifest.csv`. Never edit them manually or commit them
  without an explicit project-wide large-data decision.
- A full CSV build is large: 30 tables, 51,599,837 rows, and 2,458,559,538
  bytes. For a development build, combine the exporter's `--skip-corpus`,
  `--skip-exact-history`, and `--skip-appendix-c` flags with a new, empty
  output directory. Do not mix partial and complete snapshots.
- `scripts/matlab_v73_reader.py` is an intentional, dependency-light reader
  for the old MATLAB HDF5 structures in `redditData.mat`. SciPy cannot read
  MATLAB v7.3 files, and this reader intentionally avoids an `h5py`
  dependency; do not replace it casually with a general HDF5 loader.

## Modeling boundaries

- State which semantics an analysis uses: `release_compat`, `paper_intended`,
  or `corrected_common_engine`. A release quirk is not automatically an
  implementation bug, and a corrected result is not an exact replay.
- The behavioral release contains aggregate condition means, not participant
  responses, trial denominators, or standard errors. Do not invent binomial
  denominators or participant-level uncertainty.
- The authors did not release the producer for `Combined.mat/counts225` or the
  original tokenizer. Rebuilding it from raw event streams is a replication
  target, not an exact conversion step.
- Exact stochastic A&M replay is impossible because seeds and parameter-search
  history were not released. Compare stochastic outputs with declared seeds
  and uncertainty rather than requiring bit identity.
- The ten cached Hick micro-window odds and the 32-by-32 environmental surface
  are distinct targets. Do not infer the missing micro-window construction
  from the larger surface.
- Do not run `external/anderson_2023/matlab/temp.m` as part of reproduction; it
  is an unrelated, extremely heavy matrix-multiplication benchmark.

## Documentation and verification

- In Markdown, wrap inline math in single dollar signs (for example, `$x$`)
  and put display math between double-dollar delimiter lines. Do not use
  backslash-parenthesis or backslash-bracket math delimiters in Markdown.
  Native LaTeX files such as `notes.tex` keep normal LaTeX conventions.
- The repository has no pinned Python environment. The full test and CSV
  pipeline requires NumPy and SciPy; use `python3` on systems without a
  `python` alias.
- Before publishing model changes, run:

  ```bash
  python3 -m unittest discover -s tests -v
  python3 scripts/sanity_check_anderson_2023.py
  git lfs fsck
  ```

- For data-pipeline changes, also build the affected tables in a new output
  directory and validate that directory. Before claiming full-snapshot
  compatibility, regenerate the complete snapshot and run
  `python3 scripts/prepare_anderson_2023_data.py --validate-only`.
- `--validate-only` expects all tables recorded in that directory's manifest.
  Keep partial development builds separate from the complete snapshot.
- Keep large generated CSVs and local handoff/export files out of commits.

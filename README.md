# Optimal spacing and memory

This repository studies how the spacing between learning episodes affects later memory. It contains the working manuscript and simulations for a simple spacing model, data recovered from the Cepeda et al. literature, and Python implementations and converted data for the environmental memory models in Anderson et al. (2023).

## Repository map

- `notes.tex` is the working manuscript on learning, forgetting, and optimal study gaps.
- `data/` contains the small, tracked datasets recovered from published spacing-effect figures.
- `src/` contains plotting and fitting code, the original R simulation, and independent Python implementations of the Anderson et al. models.
- `results/` contains fitted parameters, predictions, and modeling notes.
- `external/anderson_2023/` is an immutable copy of the authors' MATLAB release, with provenance and SHA-256 checksums.
- `scripts/prepare_anderson_2023_data.py` converts that release into validated, relational CSV tables.
- `docs/` contains the Anderson et al. model reference, focused implementation notes, data dictionary, CSV schemas, and minimal replication plan.
- `tests/` contains model, conversion-pipeline, and MATLAB v7.3 reader tests.

## Python environment

The scripts require Python 3.10 or newer. Create and activate a repository-local
virtual environment, then install the runtime dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run scripts from the repository root. For example:

```bash
python src/fit_sac_mcm_2008_response_mapping.py
python src/fit_mcm_cepeda2008.py --starts 8
python src/fit_mcm_cepeda2009.py --starts 8
python src/plot_mcm_replication.py
python src/plot_optimal_gaps.py
python src/plot_spacing_recall.py
```

The fitting scripts write CSV files to `results/`; the plotting scripts write
figures to `figures/`. Some fit commands use many optimization starts by
default and can take several minutes.

## Clone and retrieve Git LFS files

Install [Git LFS](https://git-lfs.com/) before cloning. A normal clone will then retrieve the 16 `.mat` and `.docx` objects (about 357 MB) automatically:

```bash
git lfs install
git clone https://github.com/popov-lab/optimal-spacing-memory.git
cd optimal-spacing-memory
git lfs pull
git lfs fsck
```

If the repository was cloned before Git LFS was installed, the affected files will contain small text pointers instead of MATLAB or Word data. From the repository root, recover them with:

```bash
git lfs install
git lfs pull
git lfs checkout
git lfs fsck
```

Do this before diagnosing a `.mat` parsing error or a suspiciously small file under `external/anderson_2023/matlab/`.

## Python verification

Run the fast verification suite from the repository root:

```bash
python3 -m unittest discover -s tests -v
python3 scripts/sanity_check_anderson_2023.py
git lfs fsck
```

The custom MATLAB v7.3 reader deliberately does not require MATLAB or `h5py`.

## Generated Anderson et al. CSV data

Only `data/derived/anderson_2023/manifest.csv` is tracked. The 30 generated CSV tables are deterministic but intentionally ignored because a complete export contains 51,599,837 rows and occupies 2,458,559,538 bytes.

Rebuild or validate the complete snapshot with:

```bash
python3 scripts/prepare_anderson_2023_data.py
python3 scripts/prepare_anderson_2023_data.py --validate-only
```

For development, the exporter also provides `--skip-corpus`, `--skip-exact-history`, and `--skip-appendix-c`. Use those flags with a new, empty output directory; the exporter rejects mixed partial and complete snapshots. See [the data reference](docs/anderson_2023_data_reference.md) and [CSV schema](docs/anderson_2023_csv_schema.md) before interpreting a table.

## Anderson et al. replication notes

The vendored MATLAB directory is a byte-for-byte source snapshot; fixes and adapters belong elsewhere. The Python module implements the corrected equations used for new work and exposes a separate Anderson PPE mapping and A&M release option only where they are needed to understand the released results.

The release does not include the code that produced `Combined.mat/counts225` or random seeds for the stochastic A&M fits. The directly relevant paper-to-code differences are listed in the [implementation notes](docs/anderson_2023_implementation_notes.md); the next validation steps are deliberately limited in the [replication plan](docs/anderson_2023_replication_plan.md).

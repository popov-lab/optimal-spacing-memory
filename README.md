# Optimal spacing memory

## Python environment

The scripts in `src/` require Python 3.10 or newer. Create and activate a
repository-local virtual environment, then install the runtime dependencies:

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

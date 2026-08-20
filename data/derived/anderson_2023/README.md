# Anderson et al. (2023) generated CSV data

These files are deterministic derivatives of the immutable author release in `external/anderson_2023/matlab/`.  Do not edit generated CSVs manually.

Rebuild and verify them from the repository root with:

```bash
python scripts/prepare_anderson_2023_data.py
python scripts/prepare_anderson_2023_data.py --validate-only
```

`manifest.csv` records each generated table's row count, byte size, SHA-256, and source lineage.  The schemas, table meanings, validation invariants, and known release limitations are documented in `docs/anderson_2023_data_reference.md` and `docs/anderson_2023_csv_schema.md`.

# Anderson et al. (2023) author release

This directory is a byte-for-byte import of the MATLAB code and data released with:

> Anderson, J. R., Betts, S., Byrne, M. D., Schooler, L. J., & Stanley, C. (2023). The environmental basis of memory. *Psychological Review, 130*, 1137-1166. https://doi.org/10.1037/rev0000409

## Provenance

- Source page: https://act-r.psy.cmu.edu/?p=32939&post_type=publications
- Original archive: https://act-r.psy.cmu.edu/wordpress/wp-content/uploads/2021/10/paper_matlab.zip
- Imported: 2026-08-20
- Upstream `Last-Modified`: 2023-03-02 14:47:05 GMT
- Upstream ETag: `"2001a5e-a060af3-5f5ebe6fccc40"`
- Archive size: 168,168,179 bytes
- Archive SHA-256: `b3618e8cfa724fc9e321a464f32317920e4a82ee9a96101292efdc5a294a318c`

The 53 files in `matlab/` were copied unchanged. `SHA256SUMS` records their individual hashes. Only packaging metadata was omitted: the ZIP's `__MACOSX` resource-fork directory, `.DS_Store`, and the transient Word lock file `~$readme.docx`.

The supplied `.mat` and `.docx` files are tracked with Git LFS because `redditData.mat` exceeds GitHub's ordinary per-file limit.

Treat `matlab/` as an immutable vendor snapshot; put fixes and adapters elsewhere. The release has no run-all entry point or recorded MATLAB version. Its scripts require the Statistics and Machine Learning and Parallel Computing toolboxes, and the preprocessing scripts as written also call Wavelet Toolbox's `wrev`. Some fitting and simulation entry points are computationally large and unseeded. In particular, `temp.m` is an unrelated, extremely heavy matrix-multiplication benchmark and should not be run as part of reproduction.

## Relationship to this repository

- `src/anderson_2023_models.py` is the independent Python implementation.
- `scripts/prepare_anderson_2023_data.py` creates validated relational CSVs.
- `scripts/matlab_v73_reader.py` streams the released MATLAB 7.3 Reddit data without requiring MATLAB or an HDF5 Python package.
- `docs/anderson_2023_data_reference.md` documents every generated table.
- `docs/anderson_2023_model_reference.md` transcribes the model equations.
- `docs/anderson_2023_release_audit.md` compares the paper, this release, and the independent implementation.

The author release has no standalone license file. Retain its attribution and consult the source page or authors before redistributing it outside the research purposes for which it was shared.

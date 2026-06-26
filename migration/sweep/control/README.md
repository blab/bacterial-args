# control — classical single-tree phylogeography

The classical baseline to compare GAIA against: mask recombination (Gubbins), build one tree
(IQ-TREE), time-scale it (`augur refine`), and reconstruct ancestral demes (`augur traits` /
TreeTime discrete-trait analysis). Ancestral-location accuracy is then scored the same way as the
GAIA side, per pairwise MRCA.

## Contents

| Path | What |
|------|------|
| `ctrl/Snakefile` | **Primary** control pipeline: Gubbins → rescale → `augur refine` → `augur traits`; writes per-run location CSVs via `get_control_df`. |
| `snake/Snakefile` | Parallel variant of the same pipeline. |
| `preprocess_control.ipynb` | Exploratory notebook prototyping the control (tree loading, height calc, traits parsing). |
| `control_gaia.Rmd` | Scores/plots the control output alongside GAIA. |
| `ctrl/make_plots_ctrl.ipynb` | Control plots. |

## Status / caveats

- **`all.py` was removed** in cleanup — it was a broken early-draft Snakefile (Python code inside
  `shell:` blocks, mismatched I/O names) superseded by `ctrl/Snakefile` and `snake/Snakefile`.
- **MRCA height formula unified.** `ctrl/Snakefile` and `snake/Snakefile` now both use
  `height = (t.distance(a,b) / 2) / mu` (MRCA time above tips, comparable to `t.time(u)` on the
  sim/inf side). `preprocess_control.ipynb` still contains the older
  `(total_branch_length - distance) / mu`, which is **not** comparable (`total_branch_length()` is
  the sum of all branches, not tree height). Confirm the formula before re-running — see
  [`../../../KNOWN_ISSUES.md`](../../../KNOWN_ISSUES.md).
- Gubbins is not on conda-forge for macOS arm64; run the control on the Linux cluster.

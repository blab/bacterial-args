# migration — recombination-aware phylogeography

Given a simulated, geographically structured bacterial population, how well does ancestral-location
reconstruction on the ARG (**GAIA**) recover the truth, compared with a classical single-tree
**DTA** control — and where does the graph help?

## Pipeline

1. **Simulate + infer** ([`sweep/make_geo_dfs.ipynb`](sweep/make_geo_dfs.ipynb)) — an `msprime`
   demography (source-sink island model) produces a tree sequence with known per-node demes;
   `tsinfer`+`tsdate` re-infer it. For every pair of tips in every local tree, the deme of their
   MRCA is recorded (`locations/sim_*.csv`, `locations/inf_*.csv`) and the tree sequences are
   dumped (`trees/{sim,inf}_*.trees`).
2. **GAIA reconstruction** ([`sweep/sweep_gaia.Rmd`](sweep/sweep_gaia.Rmd)) — runs GAIA's discrete
   maximum-parsimony reconstruction (`treeseq_discrete_mpr`) on the tree sequences with a migration
   cost matrix, and scores each inferred MRCA location against the simulated truth (proportion
   correct per genomic bin / per node-height bin).
3. **Classical control** ([`sweep/control/`](sweep/control/)) — Gubbins → IQ-TREE → `augur refine`
   → `augur traits` on a single tree; and ([`sweep/merge_trees/`](sweep/merge_trees/)) a
   "multi-tree" DTA that runs TreeTime discrete-trait reconstruction over the local trees joined at
   a zero-length root.
4. **Analyses** — accuracy vs node depth (`sweep/downsample.ipynb` +
   `sweep/plot_logistic_gaia_downsample.Rmd`), clade recovery (`sweep/clades/`), and migration-event
   counts / tip histories (`sweep/merge_trees/`).

## Contents

| Path | What |
|------|------|
| [`sweep/`](sweep/) | **Canonical** pipeline (source-sink model). See its sub-directories below. |
| [`sweep/make_geo_dfs.ipynb`](sweep/make_geo_dfs.ipynb) | Geographic simulation + ARG inference. |
| [`sweep/sweep_gaia.Rmd`](sweep/sweep_gaia.Rmd) | GAIA reconstruction + per-bin accuracy scoring. |
| [`sweep/downsample.ipynb`](sweep/downsample.ipynb) | Node-height binning/downsampling for the depth analysis. |
| [`sweep/control/`](sweep/control/) | Classical Gubbins + IQ-TREE + DTA control. See [README](sweep/control/README.md). |
| [`sweep/merge_trees/`](sweep/merge_trees/) | Multi-tree DTA + migration-count/tip-history diagnostics. See [README](sweep/merge_trees/README.md). |
| [`sweep/clades/`](sweep/clades/) | Clade-recovery (exact + Jaccard) analysis. See [README](sweep/clades/README.md). |
| [`sweep/pooled/`, `sweep/clades/pooled/`](sweep/pooled/) | Pooled re-analyses of the same scores. |
| [`in_terminal/`](in_terminal/) | **Earlier, divergent copy** of the sim+GAIA workflow. Not canonical (see caveat). |
| [`tree_styling.ipynb`](tree_styling.ipynb) | Tree visualization helpers. |

Sweep parameters (source-sink): 3 demes, 5–20 samples/deme, `L=3e6`,
`PM_GRID`/`MU_GRID`/`MIGRATION_RATES` each `linspace(...,5)`, `N_REPS=5`.

## Status / caveats

- **`sweep/` is canonical; `in_terminal/` is an earlier copy with different demography** (source
  `Ne=4e3` vs `6e3`, sink `Ne=1e3` vs `1.5e2`, `sink_src=/10` vs `/20`, `sink_sink=*2` vs `*10`,
  background `1e-6` vs `0`) and a correspondingly different GAIA cost matrix. Don't mix results
  across the two.
- The **GAIA-vs-DTA comparison is confounded** along several axes at once (graph vs single tree,
  parsimony vs ML, and cost-matrix specification). See [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md).
- The **equal-migration model is disabled** in `make_geo_dfs.ipynb` (commented out); only the
  source-sink model runs.
- Recombination is simulated as crossover, not gene conversion — see
  [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md).

# bacterial-args

Using ancestral recombination graphs (ARGs) to perform **recombination-aware bacterial
phylogeography**.

## Motivation

Bacteria recombine frequently, so different regions of a genome can have different evolutionary
histories. Classical phylogeographic pipelines (e.g. Gubbins to mask recombinant regions →
IQ-TREE → discrete trait analysis) collapse this to a single tree and discard the masked data.
This project asks whether a **graph-based** approach — inferring an ARG / tree sequence with the
[tskit](https://tskit.dev) suite (`msprime`, `tsinfer`, `tsdate`) and reconstructing ancestral
geography on it with [GAIA](https://github.com/blueraleigh/gaia) — recovers migration history more
accurately, and **under which mutation/recombination regimes** the graph buys you anything over a
single tree.

The work is simulation-based: we simulate bacteria-like ancestries with known truth, re-infer them
from sequence data, and score the inference.

## Repository layout

```
sweep/        ARG-inference validation: can tsinfer+tsdate reconstruct the simulated ARG?
migration/    Phylogeography: GAIA on the (true / inferred) ARG vs a classical tree+DTA control.
environment.yml      conda env for the Python side (msprime / tsinfer / tsdate / treetime / augur)
install_packages.R   R packages for the GAIA + plotting side
KNOWN_ISSUES.md      current caveats and things to confirm before trusting outputs
```

- **[sweep/](sweep/)** — simulate tree sequences across a grid of mutation rate (μ) and
  recombination/mutation ratio (ρ/μ), re-infer with `tsinfer`+`tsdate`, and validate by comparing
  per-100kb-bin pairwise-MRCA (patristic) distances between the simulated and inferred trees
  (summarized as R², against a permuted-bin null). Corresponds to Figure 2 of the writeup.
- **[migration/](migration/)** — add geographic structure (demes + migration), then compare
  ancestral-location reconstruction by GAIA on the ARG against a classical merged-tree + DTA
  control. Corresponds to Figures 4–5.

## Setup

```bash
# Python side (simulation, ARG inference, DTA control)
conda env create -f environment.yml
conda activate bacterial-args

# R side (GAIA reconstruction + plotting)
Rscript install_packages.R
```

The `.Rmd` files locate inputs with `here::here(...)`, so they can be knit from anywhere in the
repo once the `here` package is installed.

> **Heads-up:** the heavy sweeps were run on a SLURM cluster. Committed CSV/figure outputs predate
> recent code fixes (notably the simulation seed bug) and should be regenerated before use — see
> **[KNOWN_ISSUES.md](KNOWN_ISSUES.md)**.

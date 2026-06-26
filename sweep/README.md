# sweep — ARG-inference validation

Can `tsinfer` + `tsdate` reconstruct a simulated bacterial ARG from sequence data, and how does
accuracy depend on the recombination/mutation ratio (ρ/μ)?

## Pipeline

1. **Simulate** a tree sequence with `msprime` (`sim_ancestry` + `sim_mutations`) for a given
   (μ, ρ) — base parameters `Ne=5000`, `L=3e6`, 60 haploid samples across 5 timepoints.
2. **Export** to VCF → VCZ (`vcf2zarr`).
3. **Infer** the ARG back with `tsinfer` (`generate_ancestors` → proxy samples → `match_ancestors`
   → `match_samples`) and date it with `tsdate`.
4. **Validate**: for every pair of samples in every local tree, take the pairwise-MRCA time
   (patristic distance). At each 100kb genomic bin, correlate simulated vs inferred distances
   (R² of `log1p` values). Compare against a **null** where inferred bins are randomly permuted.

Output tables `inf_all.csv` (inferred R² per bin) and `null_all.csv` (null R² per bin) are read by
the plotting notebooks. Expectation: high R² at low ρ/μ, degrading as recombination outpaces
mutation and each inter-breakpoint segment carries too few mutations to reconstruct.

## Contents

| Path | What |
|------|------|
| [`base/example/example.ipynb`](base/example/) | Simulate, infer, validate for a single (μ, ρ/μ). |
| [`base/make_dfs.ipynb`](base/) | The sweep over μ × ρ/μ. Produces `out/inf_all.csv`, `out/null_all.csv`. |
| [`base/make_plots.ipynb`](base/) | Figure-2 panels from the base sweep. |
| [`replicates/make_dfs_5reps.ipynb`](replicates/) | **Canonical** sweep with 5 independent replicates per cell. |
| [`replicates/make_plots_simple.ipynb`](replicates/) | Plots from the replicate sweep. |

Parameter grids (both sweeps): `MUS = [1e-10, 1e-9, 1e-8, 1e-7]`,
`PM_GRID = [1e-3, 1e-2, 0.1, 0.3, 1.0, 3.0]` (ρ/μ), 100kb bins.

## Status / caveats

- **Use `replicates/` for anything quantitative.** The base `make_dfs.ipynb` historically had a
  seed bug (every cell simulated from a single global `seed = 50`); it is now fixed to thread the
  per-combo seed, matching the replicates notebook, but it still runs only `N_REPS = 1`.
- **Committed `base/out/*.csv` and the PNG figures are stale** — generated before the seed fix, and
  the saved run only covered 2 of the 24 cells. Regenerate on the cluster.
- A `stagger` variant (inference under misspecified parameters) is described in the writeup but is
  **not implemented** in this repo.
- Recombination is simulated as crossover, not bacterial gene conversion — see
  [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md).

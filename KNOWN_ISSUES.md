# Known issues & caveats

Status as of the cleanup PR. Split into what was fixed here vs. what still needs attention before
the results can be trusted. None of the heavy pipelines were re-run during cleanup (msprime/R are
cluster-scale and not available locally), so **committed outputs are not regenerated**.

## Fixed in this cleanup

- **Simulation seed bug** (`sweep/base/make_dfs.ipynb`). `sim()` read a module-level `seed = 50`
  for every simulation, so the per-combo seed computed in `sim_sweep_parallel()` was ignored and
  the whole grid shared one random origin. `seed` is now threaded
  `run_cell` → `simulate_one` → `sim()`, matching the already-correct
  `sweep/replicates/make_dfs_5reps.ipynb`.
- **Hard-coded `/home/nahmed/...` paths.** All *active* `setwd()` / `root.dir` / absolute
  `read.csv()` calls in the `.Rmd` files now use `here::here(...)`. (A few inert references remain
  in old notebook *outputs* and one commented `os.chdir`; they clear on re-run.)
- **Dead code removed.** `migration/sweep/control/all.py` (a broken early-draft Snakefile with
  Python inside `shell:` blocks) was deleted; the working pipelines are `control/ctrl/Snakefile`
  and `control/snake/Snakefile`.
- **Control MRCA-height formula unified.** Both control Snakefiles now use
  `height = (t.distance(a,b) / 2) / mu`.
- **Reproducibility scaffolding added**: `environment.yml`, `install_packages.R`, `.gitignore`.

## Open — confirm before trusting results

- **Stale committed outputs (regenerate on cluster).** `sweep/base/out/{inf_all,null_all}.csv` and
  the Figure-2 PNGs predate the seed fix, and the saved sweep runs only covered 2 of the intended
  parameter cells. The `migration/` location CSVs / figures predate downstream fixes too.
- **Crossover vs. gene conversion** *(deferred to its own PR)*. Every `msprime.sim_ancestry(...)`
  uses `recombination_rate=` (crossover). Bacteria recombine via short **gene-conversion** tracts
  (`gene_conversion_rate` + `gene_conversion_tract_length`), which produces a structurally different
  ARG. The commented `#def sim(..., gcr, gcrl)` and `#later: tract lengths` show this was intended.
  Switching the model will change essentially all results, so it is being handled separately.
- **GAIA-vs-DTA comparison is confounded.** "Control sometimes wins" mixes three variables at once:
  substrate (ARG vs single tree), reconstruction criterion (GAIA Sankoff parsimony vs TreeTime ML /
  `augur traits`), and cost-model specification. A clean test holds the reconstruction method fixed
  and varies only the substrate. See `migration/sweep/merge_trees/` diagnostics.
- **`migration/sweep/` vs `migration/in_terminal/` diverge.** Different demography (source/sink Ne,
  migration ratios, background rate) and GAIA cost matrices. `sweep/` is canonical; don't mix.
- **`preprocess_control.ipynb` height formula** still uses `(total_branch_length - distance)/mu`,
  which is not comparable to `t.time(u)` (total_branch_length is the sum of all branches, not tree
  height). The Snakefiles were corrected; confirm and align the notebook if it is still used.
- **Non-reproducible sample counts** (`migration/sweep/make_geo_dfs.ipynb`). Per-deme sample sizes
  use `random.randint(...)` from the unseeded global `random` module, so they vary run-to-run even
  at fixed `seed` (which only seeds msprime).
- **Perfect separation in logistic fits** (`merge_trees/random_sampling/`, `clades/`). Cells with
  100%/0% accuracy make node-height-threshold metrics unreliable and not comparable across cells.
- **Not implemented:** the `stagger` (inference under misspecified parameters) variant described in
  the writeup; the equal-migration demography in `make_geo_dfs.ipynb` is present but commented out.

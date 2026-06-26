# merge_trees — multi-tree DTA and migration diagnostics

A second flavour of the classical control, plus the diagnostics used to investigate *why* results
were sometimes incoherent (cases where the single-tree control matched or beat GAIA).

## Multi-tree DTA

`multi-tree-dta.py` takes the inferred local trees (one Newick per genomic bin, named
`index<NN>*.nwk`) plus a TSV of tip→deme states, joins them into a single "multi-tree" with a
zero-length polytomy root, and runs **TreeTime** discrete-trait reconstruction once over the whole
thing. Outputs per-node states/confidence/entropy (`.states.tsv`, node-data JSON, an Auspice
dataset, and an across-tree agreement histogram). `Snakefile` drives it; `merge_fxns.ipynb` holds
helper logic.

```bash
python multi-tree-dta.py --trees <dir-of-nwk> --metadata states.tsv --output-prefix out/run
```

> Note: this is "tree + DTA over many local trees", distinct from GAIA's parsimony reconstruction
> over the tree sequence. Comparing the two conflates *substrate* (graph vs trees) with
> *method* (parsimony vs ML) — see [`../../../KNOWN_ISSUES.md`](../../../KNOWN_ISSUES.md).

## Diagnostics

| Path | What it shows |
|------|----------------|
| `random_sampling/inspect.ipynb` | The "**runs where control does better**" investigation. Logistic fits of location accuracy vs node height for `control` (tree+DTA), `inf` (inferred ARG+GAIA), and `sim` (true ARG+GAIA). Finds an inverted, non-monotonic pattern (tree better on deep nodes, network better on shallow nodes). **Root cause left unresolved** by the original analysis. |
| `random_sampling/count_migrations.ipynb` | Counts migration events per route (true vs sim vs inf vs control) across all runs. The single-tree control collapses many migrations into few; networks track counts closer to truth — but more events ≠ better per-node accuracy. |
| `random_sampling/plot_random_sampling.ipynb` | Accuracy curves and "node depth at 80% accuracy" across the parameter grid. Note the logistic fits hit perfect-separation in several cells, so that threshold metric is **not comparable** across conditions. |
| `example/count_migrations_example.ipynb` | Single-run worked example of the migration counting. |
| `example/tip_history_example.ipynb` | Traces one tip's inferred geographic ancestry back in time across genomic positions (tree vs network). Diagnostic, not a metric. |
| `plot_even_sampling.ipynb` | Even-sampling variant plots. |

## Caveats

- The diagnostics document a real tension (the graph does not uniformly beat the tree) but do not
  resolve it; see [`../../../KNOWN_ISSUES.md`](../../../KNOWN_ISSUES.md) for the likely confounds.
- Some notebook *outputs* still show old absolute paths from prior cluster runs; these clear on
  re-execution.

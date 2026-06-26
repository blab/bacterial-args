# clades — clade recovery + GAIA location accuracy

Does the inferred ARG recover the same clades (monophyletic tip sets) as the simulated ARG, and is
the GAIA-reconstructed location of each clade's MRCA correct?

## Pipeline

1. **Find clades** (`find_clades.ipynb`, `find_clades_jaccard.ipynb`) — enumerate internal nodes
   in the simulated and inferred trees per genomic bin and match them:
   - *exact*: identical tip set at the same genomic position;
   - *Jaccard*: tip-set overlap above a threshold (for the inferred-vs-simulated case).
2. **Score with GAIA** (`clades_gaia.Rmd`, `clades_gaia_jaccard.Rmd`) — run GAIA's MPR on the
   matched clades and check whether the inferred MRCA location equals the true one.
3. **Plot** (`plot_clades_gaia.Rmd`, `plot_clades_gaia_jaccard.Rmd`, and `pooled/`) — accuracy vs
   node height, logistic fits, node-height thresholds. Figures in `figures/`.

## Caveats

- Several logistic fits hit **perfect separation** (cells with 100% or 0% accuracy), so confidence
  intervals and node-height thresholds from those fits are unreliable / not comparable across cells.
- Exact-match scoring conflates topology error with tip-assignment error; the Jaccard variant
  relaxes this.
- `pooled/` re-aggregates the same scores; check which aggregation a given figure used.
- See [`../../../KNOWN_ISSUES.md`](../../../KNOWN_ISSUES.md) for cross-cutting caveats.

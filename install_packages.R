# R packages for the GAIA phylogeography and plotting side of bacterial-args.
# Run once with:  Rscript install_packages.R
#
# The .Rmd files use here::here() for paths, so run them from anywhere inside
# the repository (no setwd() needed).

cran <- c(
  "here",      # project-relative paths
  "ape",       # tree handling
  "dplyr",     # data wrangling
  "tidyr",
  "ggplot2",   # plotting
  "cowplot",
  "remotes"    # for installing gaia from GitHub
)

to_install <- cran[!(cran %in% rownames(installed.packages()))]
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

# GAIA (geographic ancestor inference on tree sequences) is not on CRAN.
# Install from GitHub. See https://github.com/blueraleigh/gaia for details.
if (!requireNamespace("gaia", quietly = TRUE)) {
  remotes::install_github("blueraleigh/gaia")
}

cat("Done. Installed:", paste(c(cran, "gaia"), collapse = ", "), "\n")

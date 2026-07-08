## setup_environment.R -- install everything needed to reproduce the PACE
## manuscript analyses. Run once from the repository root:
##
##   Rscript setup_environment.R
##
## Verified against R 4.5.0. A C/C++ toolchain is required: the PACE package's
## Rcpp contamination solver compiles from source on install. The melanoma
## notebook reads an .h5ad via zellkonverter, which needs a working
## basilisk/Python; the first run sets that up automatically.

options(repos = c(CRAN = "https://cloud.r-project.org"))

## ---- 1. CRAN packages -------------------------------------------------------
## The notebooks' own dependencies. PACE's numerical/plotting dependencies
## (Rcpp, RcppEigen, mashr, ashr, ...) are pulled in when the package installs
## in step 3.
cran <- c(
  "remotes", "BiocManager", "knitr", "rmarkdown",
  "data.table", "dplyr", "tidyr", "forcats", "stringr",
  "ggplot2", "ggrepel", "patchwork", "cowplot", "ggrastr",
  "Matrix", "FNN", "dbscan"
)
to_install <- setdiff(cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

## ---- 2. Bioconductor packages ----------------------------------------------
bioc <- c(
  "SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment",
  "zellkonverter",                         # melanoma .h5ad reader
  "Statial"                                # Statial = the SpatioMark framework
)
to_install <- setdiff(bioc, rownames(installed.packages()))
if (length(to_install)) BiocManager::install(to_install, update = FALSE, ask = FALSE)

## ---- 3. The PACE method (this manuscript's package) -------------------------
remotes::install_github("ecool50/PACE")

## ---- 4. Method-comparison packages (Figure 3) ------------------------------
## SpatioMark (Statial) installs above from Bioconductor. niche-DE installs from
## its upstream source; confirm the repository against the method's documentation
## before running the comparison notebooks (nichede_bc.qmd, spatiomark_bc.qmd,
## spillover_composite_bc.qmd). DenoIST ships as a repo-local helper
## (pace_mv/helpers/denoist_correction.R) and needs no install.
if (!requireNamespace("nicheDE", quietly = TRUE))
  message("Install niche-DE from its source, e.g. remotes::install_github(\"kaishumason/NicheDE\")")

cat("\nEnvironment setup complete. Next: Rscript fetch_data.R\n")
cat("Quarto is required to render the .qmd notebooks (https://quarto.org).\n")

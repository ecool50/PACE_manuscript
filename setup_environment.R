## setup_environment.R -- install everything needed to reproduce the PACE
## manuscript analyses. Run once from the repository root:
##
##   Rscript setup_environment.R
##
## Verified against R 4.5.0. A C/C++ toolchain is required (Rcpp/RcppEigen/TMB
## and the streaming solver compile from source). The melanoma notebook reads an
## .h5ad via zellkonverter, which needs a working basilisk/Python; the first run
## sets that up automatically.

options(repos = c(CRAN = "https://cloud.r-project.org"))

## ---- 1. CRAN packages -------------------------------------------------------
cran <- c(
  "remotes", "BiocManager", "here", "knitr", "rmarkdown",
  "data.table", "dplyr", "tidyr", "tibble", "forcats", "stringr",
  "ggplot2", "scales", "cowplot", "ggrepel", "ggrastr", "patchwork", "viridisLite",
  "Matrix", "Rcpp", "RcppEigen", "TMB",
  "mashr", "ashr", "glmmTMB", "MASS", "FNN", "dbscan", "spatstat.geom",
  "Rfast", "akima"
)
to_install <- setdiff(cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

## ---- 2. Bioconductor packages ----------------------------------------------
bioc <- c(
  "SpatialExperiment", "SummarizedExperiment", "SingleCellExperiment",
  "S4Vectors", "BiocParallel", "scran", "DropletUtils", "zellkonverter",
  "limma", "Statial"                       # Statial = the SpatioMark framework
)
to_install <- setdiff(bioc, rownames(installed.packages()))
if (length(to_install)) BiocManager::install(to_install, update = FALSE, ask = FALSE)

## ---- 3. The PACE method (this manuscript's package) -------------------------
remotes::install_github("ecool50/PACE")

## ---- 4. Method-comparison packages (Figure 4.3) ----------------------------
## niche-DE and CellChat are installed from their upstream sources. Confirm the
## exact repositories against each method's documentation before running the
## comparison notebooks (nichede_bc.qmd, spatiomark_bc.qmd, spillover_composite_bc.qmd).
if (!requireNamespace("nicheDE", quietly = TRUE))
  message("Install niche-DE from its source, e.g. remotes::install_github(\"kaishumason/NicheDE\")")
if (!requireNamespace("CellChat", quietly = TRUE))
  message("Install CellChat from its source, e.g. remotes::install_github(\"jinworks/CellChat\")")

cat("\nEnvironment setup complete. Next: Rscript fetch_data.R\n")
cat("Quarto is required to render the .qmd notebooks (https://quarto.org).\n")

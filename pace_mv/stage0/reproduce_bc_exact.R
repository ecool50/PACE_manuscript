## reproduce_bc_exact.R -- REGRESSION TEST, not a usage template.
##
## This is the ONE place cohort-specific values appear, and only because the job
## is to reproduce a *locked* artifact (streaming/mvpql_bc_streaming_exact.rds)
## bit-for-bit: the canonical 9-cell-type ORDER, the 4 flagship pairs, and the
## historical `_exact` solver config (alpha_warmup off, no early stop) are all
## passed EXPLICITLY for that reason.
##
## To apply PACE to a new dataset you do NOT write a script like this -- you call
## pace_fit_streaming(Y, coldata, celltype_col=, image_col=) with defaults (cell
## types/coords/sample read from the data). See the header of pace_core.R.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(Matrix); library(BiocParallel)
  library(FNN); library(dbscan); library(methods)
})

REPO <- "/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"
source(file.path(REPO, "scripts/zzz.R"))
.load_all(file.path(REPO, "scripts"))
source(file.path(REPO, "scripts/helpers/pace_mvpql.R"))
source(file.path(REPO, "scripts/helpers/pace_mvpql_joint_multi.R"))
source(file.path(REPO, "streaming/helpers/pace_mvpql_streaming.R"))
source(file.path(REPO, "stage0/pace_core.R"))

ymv <- readRDS(file.path(REPO, "data/breast_cancer/Y_df_for_mcsd.rds"))

## locked values (reproduction only) -----------------------------------------
TYPES <- c("T_Cell", "Stromal", "Macrophage", "Tumour", "Endothelial",
           "B_Cell", "Myoepithelial", "Mast", "Dendritic_Cell")   ## canonical ORDER
PAIRS <- list(c("Stromal", "Tumour"), c("Macrophage", "Tumour"),
              c("Endothelial", "Tumour"), c("Myoepithelial", "Tumour"))

res <- pace_fit_streaming(
  Y = ymv$Y, df = ymv$df, types = TYPES,                 ## explicit order to match the lock
  celltype_col = "cellType", image_col = "slide", coord_cols = c("x", "y"),
  h_bio = 30, h_tech = 5, eps = 90,
  contamination = "percell_hc", dispersion = "nb2",
  drop_sparse_neff = 30, within_image = TRUE,
  edge_correct = TRUE, data_informed_tau = TRUE,
  n_iter = 32L, threads = 8L, chunk_size = 128L, tau_shrinkage = "adaptive",
  alpha_warmup = 100000, early_stop_tol = 0, min_iter = 12L,  ## `_exact` = no approximations
  fuse = FALSE, return_mu = TRUE, verbose = TRUE)

shr  <- pace_shrink(res$fit, TYPES)
dec  <- pace_decompose(res$fit, res$df, res$Y, TYPES, res$X_fixed)
mcsd <- pace_top_drivers(res$fit, shr, dec, TYPES, PAIRS)      ## explicit flagship pairs

OUT <- Sys.getenv("R_OUT", unset = file.path(REPO, "stage0/bc_streaming_stage0.rds"))
saveRDS(list(fit = res$fit, shrunken_long = shr, gene_set = colnames(res$Y),
             decomposition = dec, mcsd_canonical = mcsd,
             h_tech = 5, h_bio = 30, K_tech = res$K_tech, K_bio = res$K_bio),
        OUT)
cat(sprintf("Saved %s\n", OUT))

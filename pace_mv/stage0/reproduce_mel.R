## reproduce_mel.R -- REGRESSION TEST: Melanoma through the SAME pace_fit_streaming.
##
## Proves the generalised method core handles a cohort that BC does not exercise:
## a disease/condition contrast (Responder PD vs non-PD), a second imageID RE
## block (patient-level neighbour slopes), per-image kernels (images = patients),
## the 7-CT collapse, and NB1. The cohort DATA PREP below (h5ad load, response
## coding, Tumor_a..f collapse, gene filter, image ordering) is application code;
## the METHOD is the single pace_fit_streaming() call -- the same function BC uses.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(forcats); library(tidyr)
  library(SingleCellExperiment); library(Matrix); library(BiocParallel)
  library(zellkonverter); library(scran); library(FNN); library(dbscan); library(methods)
})

REPO <- "/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"
source(file.path(REPO, "scripts/zzz.R"))
.load_all(file.path(REPO, "scripts"))
source(file.path(REPO, "scripts/helpers/pace_mvpql.R"))
source(file.path(REPO, "scripts/helpers/pace_mvpql_joint_multi.R"))
source(file.path(REPO, "streaming/helpers/pace_mvpql_streaming.R"))
source(file.path(REPO, "stage0/pace_core.R"))

## ============================================================================
## COHORT DATA PREP (ported verbatim from build_mel_streaming.R; application code)
## ============================================================================
TYPES   <- c("Tumour", "Endothelial", "Fibroblast", "Macrophage",
             "B_Cell", "T_CD8_memory", "Treg")
TUM_RAW <- c("Tumor_a","Tumor_b","Tumor_c","Tumor_d","Tumor_e","Tumor_f")

sce <- readH5AD(file = file.path(REPO, "data/simvi_melanoma/Melanoma_5612.h5ad"), reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
px_size <- 0.12028
set.seed(1994)

## pd_nonpd response coding: keep all; PD vs non-PD(SD+CR+PR).
sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR"))]
sce$BEST_RESPONSE_BY_SCAN <- ifelse(
  as.character(sce$BEST_RESPONSE_BY_SCAN) == "PD", "PD", "nonPD")
RESP_REF <- "nonPD"; RESP_CASE <- "PD"
RESP_TERM <- paste0("Responder", RESP_CASE)   # ResponderPD

df_raw <- cbind(colData(sce), t(as.matrix(assay(sce, "counts"))),
                reducedDim(sce, "spatial")) |>
  as.data.frame() |>
  dplyr::mutate(imageID = paste(SPID, fov, sep = "_"),
                V1 = V1 * px_size, V2 = V2 * px_size) |>
  dplyr::rename(x = V1, y = V2) |>
  dplyr::mutate(celltype = make.names(celltype),
                BEST_RESPONSE_BY_SCAN = droplevels(as.factor(BEST_RESPONSE_BY_SCAN)),
                cell_ID = paste0(cell_ID, "_", fov))
df_raw$celltype <- dplyr::recode(df_raw$celltype,
  'B.cell' = "B_Cell", 'T.CD8.memory' = "T_CD8_memory",
  'endothelial' = "Endothelial", 'fibroblast' = "Fibroblast", 'macrophage' = "Macrophage")
df_raw$celltype <- ifelse(df_raw$celltype %in% TUM_RAW, "Tumour", df_raw$celltype)
df_raw <- df_raw |> dplyr::filter(celltype %in% TYPES)
df_raw$imageID <- as.character(df_raw$imageID)
df_raw$cell_id <- as.character(df_raw$cellID_str)
hvg_genes <- make.names(rownames(sce))

## order by imageID (so the per-image kernel + RE row order match the canonical)
df_raw <- df_raw[order(df_raw$imageID), ]

## det >= 5% in any focal cell type
hvg_in_df <- intersect(hvg_genes, names(df_raw))
ct_chr <- as.character(df_raw$celltype)
ct_idx_lst <- lapply(TYPES, function(c) which(ct_chr == c))
max_det <- vapply(hvg_in_df, function(g) {
  v <- df_raw[[g]] > 0
  max(vapply(ct_idx_lst, function(idx) if (length(idx) == 0L) 0 else mean(v[idx]), numeric(1)))
}, numeric(1))
hvg_in_df <- hvg_in_df[max_det >= 0.05]
df_raw$nCount <- rowSums(as.matrix(df_raw[, hvg_in_df, drop = FALSE]))
df_raw <- df_raw[df_raw$nCount > 0, ]

df <- df_raw
df$Responder <- forcats::fct_relevel(as.factor(df$BEST_RESPONSE_BY_SCAN), RESP_REF)
df$.resp_dummy <- as.integer(df$Responder == RESP_CASE)
Y <- as.matrix(df[, hvg_in_df, drop = FALSE]); storage.mode(Y) <- "integer"
cat(sprintf("[prep] %d cells x %d genes; %d images; Responder: %s\n",
            nrow(df), ncol(Y), length(unique(df$imageID)),
            paste(names(table(df$Responder)), table(df$Responder), sep="=", collapse=", ")))

## ============================================================================
## THE METHOD: one call -- same function as BC, with Mel's structure as arguments.
## ============================================================================
N_ITER <- as.integer(Sys.getenv("R_N_ITER", unset = "32"))
res <- pace_fit_streaming(
  Y = Y, df = df, types = TYPES,
  celltype_col = "celltype", image_col = "imageID", coord_cols = c("x", "y"),
  h_bio = 30, h_tech = 5, eps = 90,
  contamination = "percell_hc", dispersion = "nb1",
  condition_col = "Responder",          ## <-- disease contrast (BC has none)
  kernel_per_image = TRUE,              ## <-- images = patients (BC = one section)
  image_re = "condition_slopes",        ## <-- second imageID RE block (PAT_DZ)
  drop_sparse_neff = 0, within_image = FALSE,   ## Mel canonical: neither
  edge_correct = TRUE, data_informed_tau = TRUE,
  n_iter = N_ITER, threads = 8L, chunk_size = 128L, tau_shrinkage = "adaptive",
  alpha_warmup = 100000, early_stop_tol = 0, min_iter = 12L,  ## deterministic reference
  fuse = FALSE, return_mu = TRUE, verbose = TRUE)

cat("\n[design check]\n")
for (b in res$fit$re_meta$blocks)
  cat(sprintf("  block %-9s K_terms=%2d K_groups=%2d\n", b$group_col, b$K_terms, b$K_groups))

if (N_ITER >= 8L) {
  shr  <- pace_shrink(res$fit, TYPES, resp_term = RESP_TERM)
  dec  <- pace_decompose(res$fit, res$df, res$Y, TYPES, res$X_fixed, resp_term = RESP_TERM)
  PAIRS <- list(c("Macrophage","Tumour"), c("Fibroblast","Tumour"),
                c("Endothelial","Tumour"), c("T_CD8_memory","Tumour"),
                c("Treg","Tumour"), c("B_Cell","Tumour"))
  mcsd <- pace_top_drivers(res$fit, shr, dec, TYPES, PAIRS,
                           resp_term = RESP_TERM, resp_dummy = df$.resp_dummy)
  OUT <- Sys.getenv("R_OUT", unset = file.path(REPO, "stage0/mel_streaming_stage0.rds"))
  saveRDS(list(fit = res$fit, shrunken_long = shr, gene_set = colnames(res$Y),
               decomposition = dec, mcsd_canonical = mcsd,
               h_tech = 5, h_bio = 30, K_tech = res$K_tech, K_bio = res$K_bio),
          OUT)
  cat(sprintf("Saved %s\n", OUT))
}

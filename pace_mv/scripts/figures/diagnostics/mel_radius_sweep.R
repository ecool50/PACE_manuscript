## Supplement (melanoma): near-range-radius sensitivity of the PACE contamination
## correction -- the CosMx counterpart of the breast-cancer Supp Fig S2a/b.
##
## The E^tech (contamination) kernel truncates at 3 * h_tech, so the "near-range
## radius" the manuscript reports is R = 3 * h_tech. The Xenium-calibrated
## canonical is h_tech = 5 um  =>  R = 15 um. We refit the PACE melanoma arm at
## R in {5, 10, 15, 20} um (h_tech = R / 3) and, per radius, extract the flagship
## SPP1 macrophage<-tumour slopes and a contamination-vs-genuine gene count. The
## goal is to show the SPP1 result and the contamination control behave sensibly
## as R changes, defending the 15 um choice on the CosMx platform.
##
## This script MIRRORS the canonical melanoma notebook notebooks/analysis_mel.qmd
## exactly (same data prep, same pace_fit_streaming() call, 6 cell types with the
## T-cell subtypes pooled into T_Cell). It does NOT use build_mel_streaming.R,
## whose 7-type scheme (T_CD8_memory + Treg kept separate) would not reproduce the
## canonical fit. The parallel BC template (scratchpad/s1_radius_sweep.R) likewise
## mirrors its qmd via pace_fit_streaming(), not a builder.
##
## R = 15 IS THE CANONICAL FIT: canonical_fit.rds already is the h_tech = 5 result,
## so we do NOT refit R15. We refit only R5 / R10 / R20 and read the R15 column
## (U / se_U) straight from canonical_fit.rds$fit. The canonical file is READ-ONLY.
##
## CONSISTENCY CHECK (in place of a reproduction gate): each R5/R10/R20 refit must
## share the canonical fit's gene set (colnames U) and its neighbour-slope rownames
## (rownames without "(Intercept)"), so the radii are comparable. We also confirm
## the canonical SPP1 macrophage<-tumour slopes match the manuscript
## (baseline ~ +0.007; baseline + ResponderPD interaction ~ -0.06). Any mismatch
## stops the script.

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(forcats)
  library(Matrix)
  library(BiocParallel)
  library(FNN)
  library(dbscan)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

source("scripts/zzz.R")                              # auto-loads scripts/ (core + helpers + C++)
source("streaming/helpers/pace_mvpql_streaming.R")   # streaming solver
source("stage0/pace_core.R")                         # pace_fit_streaming()

## ---- canonical model knobs (Mel per-cell-HC, NB1, per-image frNN) -----------
## Verbatim from analysis_mel.qmd, except pql_threads = 8 (safe for ~56k cells on
## this 24GB workstation; the 1e-2 gate tolerance absorbs 4-vs-8-thread jitter).
h_bio         <- 30            # biological-proximity kernel scale (um)
kernel_radius <- 3 * h_bio     # per-image frNN truncation radius (um)
px_size       <- 0.12028       # CosMx pixel -> um
min_detection <- 0.05          # gene QC: detected in >= 5% of >= 1 cell type
min_neff      <- 30            # drop (focal, neighbour) kernel cols below this
n_iter        <- 32            # PQL outer iterations
tau_shrinkage <- "adaptive"    # empirical-Bayes pooling of variance components
pql_threads   <- 8

## 6 cell types: six Tumor sub-clusters -> "Tumour"; T-cell subtypes -> "T_Cell".
cell_types      <- c("Tumour", "Endothelial", "Fibroblast", "Macrophage",
                     "B_Cell", "T_Cell")
tumour_subtypes <- c("Tumor_a", "Tumor_b", "Tumor_c", "Tumor_d", "Tumor_e", "Tumor_f")

## Responder coding (SIMVI PD-vs-nonPD): keep ALL; non-PD = SD+PR+CR (reference),
## PD = case. resp_term "ResponderPD" -> the V_R x S interaction.
resp_ref  <- "nonPD"
resp_case <- "PD"
resp_term <- paste0("Responder", resp_case)   # "ResponderPD"

h5ad_path     <- "data/simvi_melanoma/Melanoma_5612.h5ad"
canonical_fit <- "data/simvi_melanoma/canonical_fit.rds"
sweep_out     <- "data/simvi_melanoma/sweeps/mel_radius_sweep.rds"

## Near-range radii to REFIT (R = 3 * h_tech). R = 15 (h_tech = 5) is NOT refit --
## it is taken from the canonical fit below.
radii_um_refit <- c(5, 10, 20)

## ============================================================================
## 1. Load, clean, and collapse to 6 cell types (verbatim from analysis_mel.qmd).
## ============================================================================
cat("[data] loading Melanoma_5612.h5ad ...\n")
sce <- readH5AD(h5ad_path, reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
set.seed(1994)

## Keep all valid responses; PD vs non-PD (SD + PR + CR).
keep_response <- sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR")
sce <- sce[, which(keep_response)]
sce$BEST_RESPONSE_BY_SCAN <- ifelse(
  as.character(sce$BEST_RESPONSE_BY_SCAN) == "PD", "PD", "nonPD")

## Flat cell-level data frame: metadata + counts + spatial coords (px -> um).
df_raw <- cbind(colData(sce),
                t(as.matrix(assay(sce, "counts"))),
                reducedDim(sce, "spatial")) |>
  as.data.frame() |>
  mutate(imageID = paste(SPID, fov, sep = "_"),
         x = V1 * px_size,
         y = V2 * px_size) |>
  mutate(celltype = make.names(celltype),
         BEST_RESPONSE_BY_SCAN = droplevels(as.factor(BEST_RESPONSE_BY_SCAN)),
         cell_ID = paste0(cell_ID, "_", fov))

## Harmonise label spelling, pool the T-cell subtypes, then collapse tumour.
df_raw$celltype <- recode(df_raw$celltype,
  `B.cell`       = "B_Cell",
  `T.CD8.memory` = "T_Cell",
  `T.CD8.naive`  = "T_Cell",
  `Treg`         = "T_Cell",
  `endothelial`  = "Endothelial",
  `fibroblast`   = "Fibroblast",
  `macrophage`   = "Macrophage")
t_other <- unique(grep("^T\\.", df_raw$celltype, value = TRUE))
if (length(t_other)) {
  message("folding extra T-cell labels into T_Cell: ", paste(t_other, collapse = ", "))
  df_raw$celltype[df_raw$celltype %in% t_other] <- "T_Cell"
}
df_raw$celltype <- ifelse(df_raw$celltype %in% tumour_subtypes, "Tumour", df_raw$celltype)

df_raw <- df_raw |> filter(celltype %in% cell_types)
df_raw$imageID <- as.character(df_raw$imageID)
df_raw$cell_id <- as.character(df_raw$cellID_str)
panel_genes <- make.names(rownames(sce))

## ============================================================================
## 2. Gene filter + assemble Y / df (verbatim from analysis_mel.qmd).
## ============================================================================
df_raw <- df_raw[order(df_raw$imageID), ]   # per-image kernel + RE row order

genes_in_df   <- intersect(panel_genes, names(df_raw))
celltype_chr  <- as.character(df_raw$celltype)
cells_by_type <- lapply(cell_types, function(ct) which(celltype_chr == ct))
max_detection <- numeric(length(genes_in_df)); names(max_detection) <- genes_in_df
chunk_size <- 50L
for (start in seq.int(1L, length(genes_in_df), by = chunk_size)) {
  stop_i   <- min(start + chunk_size - 1L, length(genes_in_df))
  detected <- as.matrix(df_raw[, genes_in_df[start:stop_i], drop = FALSE]) > 0
  for (j in seq_len(stop_i - start + 1L)) {
    per_type <- vapply(cells_by_type, function(idx)
      if (length(idx) == 0L) 0 else mean(detected[idx, j]), numeric(1))
    max_detection[start + j - 1L] <- max(per_type)
  }
  rm(detected)
}
invisible(gc(verbose = FALSE))
genes_in_df   <- genes_in_df[max_detection >= min_detection]
df_raw$nCount <- rowSums(as.matrix(df_raw[, genes_in_df, drop = FALSE]))
df_raw        <- df_raw[df_raw$nCount > 0, ]

df <- df_raw
df$celltype    <- factor(as.character(df$celltype), levels = cell_types)
df$imageID     <- factor(as.character(df$imageID))
df$Responder   <- forcats::fct_relevel(as.factor(df$BEST_RESPONSE_BY_SCAN), resp_ref)
df$.resp_dummy <- as.integer(df$Responder == resp_case)
Y <- as.matrix(df[, genes_in_df, drop = FALSE]); storage.mode(Y) <- "integer"

cat(sprintf("[data] %s cells, %d genes, %d images\n",
            format(nrow(df), big.mark = ","), length(genes_in_df), nlevels(df$imageID)))

## ============================================================================
## 3. Fit the PACE melanoma arm at one near-range radius.
##    Same pace_fit_streaming() call as analysis_mel.qmd; only h_tech varies.
## ============================================================================
fit_pace_at_radius <- function(radius_um) {
  h_tech <- radius_um / 3                       # E^tech truncates at 3 * h_tech
  cat(sprintf("\n[fit] near-range R = %g um  (h_tech = %.4f)  ...\n", radius_um, h_tech))
  res <- pace_fit_streaming(
    Y = Y, df = df, types = cell_types,
    celltype_col = "celltype", image_col = "imageID", coord_cols = c("x", "y"),
    h_bio = h_bio, h_tech = h_tech, eps = kernel_radius,
    contamination     = "percell_hc",   # mu = mu_bio + rho_i * a_ig
    dispersion        = "nb1",          # Mel canonical: NB1
    condition_col     = "Responder",    # PD vs non-PD disease contrast
    kernel_per_image  = TRUE,           # images = patients (no cross-FOV neighbours)
    image_re          = "intercept",    # imageID random intercept (donor effect)
    drop_sparse_neff  = min_neff,
    within_image      = TRUE,
    edge_correct      = TRUE,
    data_informed_tau = TRUE,
    n_iter = n_iter, threads = pql_threads, tau_shrinkage = tau_shrinkage,
    verbose = TRUE)
  list(radius_um = radius_um, h_tech = h_tech, U = res$fit$U, se_U = res$fit$se_U)
}

## Sequential (NOT concurrent): one fit at a time to stay within memory.
sweep <- list()
for (r in radii_um_refit) {
  sweep[[paste0("R", r)]] <- fit_pace_at_radius(r)
}

## ---- R = 15 taken directly from the canonical fit (NOT refit) --------------
canon_fit <- readRDS(canonical_fit)$fit
sweep$R15  <- list(radius_um = 15, h_tech = 5,
                   U = canon_fit$U, se_U = canon_fit$se_U)

## keep the sweep ordered by radius for tidy tables/plots
sweep <- sweep[c("R5", "R10", "R15", "R20")]

## ============================================================================
## 4. CONSISTENCY CHECK (in place of a reproduction gate).
##    (a) each refit shares the canonical gene set + neighbour-slope rownames.
##    (b) the canonical SPP1 macrophage<-tumour slopes match the manuscript.
## ============================================================================
canon_genes       <- colnames(canon_fit$U)
canon_slope_rows  <- grep("\\(Intercept\\)", rownames(canon_fit$U),
                          invert = TRUE, value = TRUE)
cat("\n===== CONSISTENCY CHECK =====\n")
for (rk in c("R5", "R10", "R20")) {
  s <- sweep[[rk]]
  genes_match <- identical(colnames(s$U), canon_genes)
  rows_match  <- all(canon_slope_rows %in% rownames(s$U))
  cat(sprintf("  %-4s gene set == canonical: %s ; canonical slope rows present: %s\n",
              rk, genes_match, rows_match))
  if (!genes_match || !rows_match)
    stop(sprintf("CONSISTENCY CHECK FAILED for %s: gene set or slope rownames ",
                 rk),
         "do not align with the canonical fit. Aborting.")
}

## canonical SPP1 sanity: baseline ~ +0.007, baseline + interaction ~ -0.06.
spp1_base_canon <- canon_fit$U["Macrophage::Tumour", "SPP1"]
spp1_int_canon  <- canon_fit$U[paste0("Macrophage::", resp_term, ":Tumour"), "SPP1"]
cat(sprintf("  canonical SPP1 Macrophage<-Tumour: baseline (non-PD) = %+.4f ; PD (baseline+interaction) = %+.4f\n",
            spp1_base_canon, spp1_base_canon + spp1_int_canon))
cat("  (manuscript: non-PD +0.007, PD -0.06)\n")
cat("CONSISTENCY CHECK PASSED.\n")

## ============================================================================
## 5. Save the sweep (list R5/R10/R15/R20; each radius_um, h_tech, U, se_U).
##    R5/R10/R20 are refits; R15 is copied from the (untouched) canonical fit.
## ============================================================================
dir.create(dirname(sweep_out), recursive = TRUE, showWarnings = FALSE)
saveRDS(sweep, sweep_out)
cat(sprintf("\nSaved sweep: %s\n", normalizePath(sweep_out)))

## ============================================================================
## 6. Cell-type markers via one-vs-rest AUC > 0.70 on CP10k (Methods:
##    "Contamination correction evaluation"). AUC is rank-based, so the per-cell
##    CP10k rescaling matters but a log transform would not; we rank CP10k directly.
##    A gene is a marker of a type if its one-vs-rest AUC for that type exceeds 0.70.
## ============================================================================
auc_threshold <- 0.70
cp10k_ranks_per_gene <- function(counts_col, libsize) {
  rank(counts_col / libsize, ties.method = "average")   # *1e4 is a monotone no-op
}
one_vs_rest_auc <- function(ranks, in_group) {
  n1 <- sum(in_group); n2 <- length(ranks) - n1
  if (n1 == 0 || n2 == 0) return(NA_real_)
  (sum(ranks[in_group]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

libsize      <- df$nCount
is_tumour    <- df$celltype == "Tumour"
is_macro     <- df$celltype == "Macrophage"
gene_names   <- colnames(Y)
auc_tumour   <- numeric(length(gene_names)); names(auc_tumour) <- gene_names
auc_macro    <- numeric(length(gene_names)); names(auc_macro)  <- gene_names
for (g in seq_along(gene_names)) {
  ranks_g          <- cp10k_ranks_per_gene(Y[, g], libsize)
  auc_tumour[g]    <- one_vs_rest_auc(ranks_g, is_tumour)
  auc_macro[g]     <- one_vs_rest_auc(ranks_g, is_macro)
}
is_tumour_marker <- auc_tumour > auc_threshold   # neighbour-marker  -> contamination
is_macro_marker  <- auc_macro  > auc_threshold   # focal-marker      -> genuine
cat(sprintf("\n[markers] Tumour markers (AUC>0.70): %d   Macrophage markers: %d   (of %d genes)\n",
            sum(is_tumour_marker), sum(is_macro_marker), length(gene_names)))

## ============================================================================
## 7. Per-radius extraction:
##    (a) SPP1 macrophage<-tumour slopes: baseline (non-PD) and PD (baseline +
##        ResponderPD interaction).
##    (b) contamination control: # significant macrophage<-tumour genes, split
##        into neighbour-marker (contamination) vs focal-marker (genuine).
##        Significance = |U/se_U| > Bonferroni z over the genes tested.
## ============================================================================
baseline_row <- "Macrophage::Tumour"
resp_row     <- paste0("Macrophage::", resp_term, ":Tumour")  # Macrophage::ResponderPD:Tumour
n_genes      <- ncol(sweep$R15$U)
z_bonf       <- qnorm(1 - 0.05 / (2 * n_genes))               # Bonferroni over genes

spp1_table <- lapply(sweep, function(s) {
  b_base <- s$U[baseline_row, "SPP1"]
  b_int  <- s$U[resp_row,     "SPP1"]
  data.frame(radius_um       = s$radius_um,
             h_tech          = s$h_tech,
             spp1_nonPD      = b_base,             # baseline slope
             spp1_interaction = b_int,             # ResponderPD interaction
             spp1_PD         = b_base + b_int)     # PD slope = baseline + interaction
}) |> bind_rows()

contam_table <- lapply(sweep, function(s) {
  z   <- s$U[baseline_row, ] / s$se_U[baseline_row, ]
  sig <- is.finite(z) & abs(z) > z_bonf
  data.frame(radius_um     = s$radius_um,
             h_tech        = s$h_tech,
             n_sig         = sum(sig, na.rm = TRUE),
             contamination = sum(sig & is_tumour_marker[names(z)], na.rm = TRUE),
             genuine       = sum(sig & is_macro_marker[names(z)],  na.rm = TRUE))
}) |> bind_rows()

cat("\n===== SPP1 macrophage<-tumour slopes by near-range radius =====\n")
print(spp1_table, row.names = FALSE)
cat(sprintf("\n(manuscript R15: non-PD +0.007, PD -0.06)\n"))

cat(sprintf("\n===== Contamination control (Bonferroni z = %.3f over %d genes) =====\n",
            z_bonf, n_genes))
print(contam_table, row.names = FALSE)
cat("\nExpectation: contamination (neighbour-marker) count falls as R grows;",
    "genuine (focal-marker) count stays ~flat.\n")

## build_wta_fit.R — PACE-MV fit for the 10x Atera WTA FFPE breast-cancer SPE.
##
## Same model as the canonical BC fit (bc_fit_joint_kernel_redesign.R): per-cell-HC
## contamination, NB2, frNN kernels, E^tech + edge correction, data-informed tau,
## celltype-only RE (single section, no Responder). Differs only in the INPUT: reads
## the WTA SpatialExperiment (18,028 genes; celltype already annotated) instead of the
## 313-panel Y_df bundle.
##
## CALIBRATION knob: R_N_GENES=<k> fits a random k-gene subset (seeded) to measure
## the per-gene time + peak memory before committing to the full ~12.5k-gene run.
## MEMORY WARNING: the solver densifies Y and builds a dense E^tech (n_cells x n_genes).
## At the full gene count that is ~25 GB -> the full fit needs gene-batching; the
## subset keeps both matrices small.
##
## Env: R_E_TECH=1 R_K_WITHIN=1 R_DROP_SPARSE_K=1 R_NEFF_MIN=30 R_DATA_INFORMED_TAU=1 \
##      R_EDGE_CORRECT=1 R_BLEED_PERCELL=1 R_DISP_MODEL=nb2 R_N_ITER=32 R_PQL_THREADS=4 \
##      [R_N_GENES=<k>]  [R_SPE=<path>]  [R_OUT=<path>]

suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(BiocParallel); library(dbscan)
  library(SpatialExperiment); library(SingleCellExperiment); library(SummarizedExperiment)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
source("scripts/zzz.R"); .load_all("scripts")
source("scripts/helpers/pace_mvpql.R")
source("scripts/helpers/pace_mvpql_joint_multi.R")

ts  <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...) cat(sprintf("[%s] %s\n", ts(), sprintf(...)))

t0     <- Sys.time()
TYPES  <- c("T_Cell","Stromal","Macrophage","Tumour","Endothelial",
            "B_Cell","Myoepithelial","Mast","Dendritic_Cell")
H_TECH <- as.numeric(Sys.getenv("R_H_TECH", unset = "5"))
H_BIO  <- as.numeric(Sys.getenv("R_H_BIO",  unset = "30"))
DET_MIN <- as.numeric(Sys.getenv("R_DET_MIN", unset = "0.05"))
N_GENES <- as.integer(Sys.getenv("R_N_GENES", unset = "0"))   # 0 = all genes
SPE_PATH <- Sys.getenv("R_SPE",
  unset = "/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/data/breast_cancer/spe_wta_ffpe_breast_cancer.rds")
OUT <- Sys.getenv("R_OUT",
  unset = if (N_GENES > 0) sprintf("data/breast_cancer/sweeps/mvpql_wta_calib%d.rds", N_GENES)
          else "data/breast_cancer/sweeps/mvpql_wta_percell_hc.rds")

## ---- 1. Load SPE; sparse gene-filter FIRST (avoid densifying 18k x 167k) ----
say("[1/7] Loading WTA SPE: %s", SPE_PATH)
spe <- readRDS(SPE_PATH)
counts <- assay(spe, "counts")                                  # genes x cells, sparse
ct_chr <- as.character(spe$celltype)

say("    detection filter (max-per-celltype >= %.2f) on sparse matrix ...", DET_MIN)
det <- vapply(TYPES, function(c) {
  idx <- which(ct_chr == c)
  if (!length(idx)) rep(0, nrow(counts)) else as.numeric(Matrix::rowMeans(counts[, idx, drop = FALSE] > 0))
}, numeric(nrow(counts)))
max_det <- apply(det, 1, max)
keep_g  <- which(max_det >= DET_MIN)
say("    %d / %d genes pass detection >= %.2f", length(keep_g), nrow(counts), DET_MIN)
counts <- counts[keep_g, , drop = FALSE]
n_full_genes <- nrow(counts)

## Library size = total over the FITTED gene set (the offset). Drop genuinely
## zero-count cells (log offset undefined). Both are computed on the full gene
## set BEFORE any calibration subset, so a gene subset never starves the offset.
libsize   <- Matrix::colSums(counts)
keep_cell <- which(libsize > 0)
if (length(keep_cell) < ncol(counts))
  say("    dropped %d zero-count cells", ncol(counts) - length(keep_cell))
counts  <- counts[, keep_cell, drop = FALSE]
libsize <- libsize[keep_cell]
ct_chr  <- ct_chr[keep_cell]
spe     <- spe[, keep_cell]

if (N_GENES > 0 && N_GENES < nrow(counts)) {
  set.seed(1)
  sub <- sort(sample(nrow(counts), N_GENES))
  counts <- counts[sub, , drop = FALSE]
  say("    CALIBRATION: random %d-gene subset (offset = full %d-gene libsize)",
      N_GENES, n_full_genes)
}

Y <- t(as.matrix(counts))                                      # cells x genes (dense)
storage.mode(Y) <- "integer"
rm(counts); invisible(gc())

df <- data.frame(
  cell_id  = as.character(spe$cell_id),
  x        = spatialCoords(spe)[, 1],
  y        = spatialCoords(spe)[, 2],
  celltype = factor(ct_chr, levels = TYPES),
  slide    = as.character(spe$slide),
  stringsAsFactors = FALSE
)
df$imageID <- factor(df$slide)
df$nCount  <- libsize                                          # full-geneset library size (valid offset)
say("    Y: %d cells x %d genes; %d slide(s); peak gc after load:", nrow(Y), ncol(Y), nlevels(df$imageID))
print(gc(reset = TRUE)[, "max used"])

## ---- 2. frNN kernel-weighted neighbour predictors (uncapped) ----
RAD <- 3 * H_BIO
say("[2/7] frNN kernels (eps=%d um, uncapped) ...", RAD)
coords <- as.matrix(df[, c("x", "y")]); n <- nrow(df); ct <- as.character(df$celltype)
fr    <- dbscan::frNN(coords, eps = RAD)
i_vec <- rep.int(seq_len(n), lengths(fr$id))
j_vec <- unlist(fr$id, use.names = FALSE); d_vec <- unlist(fr$dist, use.names = FALSE)
jc    <- match(ct[j_vec], TYPES); ok <- !is.na(jc)
i_ok  <- i_vec[ok]; d_ok <- d_vec[ok]; jc_ok <- jc[ok]
say("    %.1f neighbours/cell median; %d contributions", stats::median(lengths(fr$id)), length(i_ok))
K_tech <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok / H_TECH), dims = c(n, length(TYPES))))
K_bio  <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok^2 / H_BIO^2), dims = c(n, length(TYPES))))
colnames(K_tech) <- paste0(TYPES, "_near"); colnames(K_bio) <- TYPES
rm(fr, i_vec, j_vec, d_vec, jc, ok, i_ok, d_ok, jc_ok); invisible(gc())
for (tc in TYPES) df[[tc]] <- K_bio[, tc]

## ---- 3. drop-sparse-K + within-(slide,celltype) centering ----
if (nzchar(Sys.getenv("R_DROP_SPARSE_K"))) {
  neff_min <- as.numeric(Sys.getenv("R_NEFF_MIN", unset = "50"))
  say("[3/7] drop-sparse-K (n_eff < %.0f) ...", neff_min)
  for (f in TYPES) {
    cf <- which(df$celltype == f); if (!length(cf)) next
    for (nb in TYPES) {
      Kc <- df[[nb]][cf] - mean(df[[nb]][cf], na.rm = TRUE)
      msq <- max(Kc^2, na.rm = TRUE); neff <- if (msq > 0) sum(Kc^2, na.rm = TRUE) / msq else 0
      if (neff < neff_min) df[[nb]][cf] <- 0
    }
  }
}
if (nzchar(Sys.getenv("R_K_WITHIN"))) {
  say("    R_K_WITHIN: within-(slide, celltype) centering")
  for (tc in TYPES) df[[tc]] <- df[[tc]] - ave(df[[tc]], df$imageID, df$celltype, FUN = mean)
}

## ---- 4. X_fixed (intercept-only under E^tech; no Responder) ----
offset_vec <- log(df$nCount)
X_fixed <- matrix(1, nrow(df), 1, dimnames = list(NULL, "(Intercept)"))
re_specs <- list(list(group_col = "celltype",
                      formula = stats::as.formula(paste0("~ 1 + ", paste(TYPES, collapse = " + ")))))

## ---- 4b. per-cell-HC anchors + E^tech + edge correction ----
type_means <- t(vapply(TYPES, function(tt) colMeans(Y[df$celltype == tt, , drop = FALSE]), numeric(ncol(Y))))
rownames(type_means) <- TYPES; colnames(type_means) <- colnames(Y)

say("[4/7] homotypic-core anchors (>=50%% same-type within 30um) ...")
same_frac <- rep(NA_real_, nrow(df))
for (s in unique(df$slide)) {
  si <- which(df$slide == s); if (length(si) < 50) next
  ctl <- df$celltype[si]; frh <- dbscan::frNN(as.matrix(df[si, c("x","y")]), eps = 30)
  for (k in seq_along(si)) { idk <- frh$id[[k]]; same_frac[si[k]] <- if (length(idk)) mean(ctl[idk] == ctl[k]) else 1 }
}
core <- which(same_frac >= 0.5); core_means <- type_means
for (X in TYPES) { idx <- core[df$celltype[core] == X]; if (length(idx) >= 20) core_means[X, ] <- colMeans(Y[idx, , drop = FALSE]) }
owner_mean <- apply(core_means, 2, max); owner_t <- TYPES[apply(core_means, 2, which.max)]
percell_anchor_mask <- matrix(0, nrow(df), ncol(Y))
for (X in TYPES) {
  is_anc <- owner_t != X & owner_mean > 0.1 & (core_means[X, ] / pmax(owner_mean, 1e-9) < 0.1)
  rr <- which(df$celltype == X); if (length(rr)) percell_anchor_mask[rr, ] <- matrix(as.numeric(is_anc), length(rr), ncol(Y), byrow = TRUE)
}

say("[4b] E^tech (cross-celltype expression-weighted, per slide) ... DENSE %dx%d = %.1f GB",
    nrow(df), ncol(Y), nrow(df) * ncol(Y) * 8 / 1e9)
t_et <- Sys.time(); rad <- 3 * H_TECH
E_tech <- matrix(0, nrow(df), ncol(Y), dimnames = list(NULL, colnames(Y)))
by_slide <- split(seq_len(nrow(df)), df$imageID); ct_all <- as.character(df$celltype)
for (si in seq_along(by_slide)) {
  rows <- by_slide[[si]]; if (length(rows) < 2) next
  nn_e <- dbscan::frNN(as.matrix(df[rows, c("x","y")]), eps = rad)
  Y_im <- Y[rows, , drop = FALSE]; ct_im <- ct_all[rows]
  for (i in seq_along(rows)) {
    nl <- nn_e$id[[i]]; if (!length(nl)) next
    keep <- ct_im[nl] != ct_im[i]; if (!any(keep)) next
    E_tech[rows[i], ] <- as.numeric(crossprod(exp(-nn_e$dist[[i]][keep] / H_TECH), Y_im[nl[keep], , drop = FALSE]))
  }
}
say("    E^tech in %.1f s; range [%.2f, %.2f]", as.numeric(difftime(Sys.time(), t_et, units = "secs")), min(E_tech), max(E_tech))
if (nzchar(Sys.getenv("R_EDGE_CORRECT"))) {
  af <- numeric(nrow(df))
  for (si in seq_along(by_slide)) {
    rows <- by_slide[[si]]; cim <- as.matrix(df[rows, c("x","y")])
    af[rows] <- area_fraction(cim, rad, min(cim[,1]), max(cim[,1]), min(cim[,2]), max(cim[,2]))
  }
  E_tech <- sweep(E_tech, 1, af, "/")
  say("    edge-corrected E^tech; af median=%.3f", median(af))
}

## ---- 5. fit ----
PQL_THREADS <- as.integer(Sys.getenv("R_PQL_THREADS", unset = "4"))
data_informed_W <- NULL
if (nzchar(Sys.getenv("R_DATA_INFORMED_TAU"))) {
  re_tmp <- build_random_design_multi(df, re_specs)
  data_informed_W <- .compute_data_informed_weights(re = re_tmp, Y = Y, df = df,
    focals = re_tmp$blocks[[1]]$group_levels, TYPES = TYPES, celltype_col = "celltype", verbose = TRUE)
  rm(re_tmp)
}
say("[5/7] fit_pace_mvpql_joint_multi (%d genes, %d cells, %d threads) ...", ncol(Y), nrow(Y), PQL_THREADS)
print(gc(reset = TRUE)[, "max used"])
t_fit <- Sys.time()
fit <- fit_pace_mvpql_joint_multi(
  Y = Y, X_fixed = X_fixed, df = df, re_specs = re_specs, offset_vec = offset_vec,
  pre_offset_mat = NULL, data_informed_W = data_informed_W,
  ambient_mat = E_tech, ambient_image_idx = as.integer(df$imageID), ambient_n_images = nlevels(df$imageID),
  bleed_percell = TRUE, percell_anchor_mask = percell_anchor_mask,
  n_iter = as.integer(Sys.getenv("R_N_ITER", unset = "32")), tol = 5e-3,
  tau_shrinkage = Sys.getenv("R_TAU", unset = "adaptive"),
  BPPARAM = SerialParam(), n_threads = PQL_THREADS,
  interior_precision = 1L, chunk_size = 128L, verbose = TRUE)
fit_min <- as.numeric(difftime(Sys.time(), t_fit, units = "mins"))
say("    FIT done in %.1f min (%.3f min/gene)", fit_min, fit_min / ncol(Y))
say("    peak memory during/after fit:"); print(gc()[, "max used"])

saveRDS(list(fit = fit, gene_set = colnames(Y), K_bio = K_bio, h_tech = H_TECH, h_bio = H_BIO,
             n_genes = ncol(Y), n_cells = nrow(Y), fit_minutes = fit_min,
             cell_meta = data.frame(cell_id = df$cell_id, nCount = df$nCount,
                                    celltype = as.character(df$celltype))), OUT)
say("Saved %s", OUT)
say("TOTAL %.1f min for %d genes -> full %d genes ~ extrapolate fit phase linearly",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), ncol(Y), length(keep_g))

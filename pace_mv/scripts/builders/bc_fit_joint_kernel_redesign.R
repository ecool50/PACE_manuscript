## bc_fit_joint_kernel_redesign.R — BC canonical fit (per-cell-HC contamination).
##
## Fits the BC canonical PACE-MV model (`mvpql_percell_hc.rds`): per-cell
## homotypic-core contamination (μ = μ_bio + ρ_i·a_ig, where ρ_i is an EB-shrunk
## per-cell loading on the RAW cross-celltype E^tech ambient field a_ig, anchored
## on homotypic-core negative-control genes), isotropic edge correction, frNN
## kernel, data-informed τ, NB2 dispersion.  Supersedes the rank-1 joint-bleed
## and two-stage NB1 prefit (`mvpql_etech.rds`) paths, both removed 2026-06-15.
##
## KEY DIFFERENCE vs Mel: BC has NO Responder (single biopsy, no disease
## contrast).  Its 17 *slides* are technical FOVs of one section, so:
##   - the per-cell ambient field a_ig is built per df$slide (n_images = 17),
##   - the RE structure stays celltype-only (NO imageID RE block) — the slides
##     are technical, not biological replication, so we do not model them as a
##     variance component.
##
## Runs through the multi-block solver `fit_pace_mvpql_joint_multi`.  Single-
## celltype-block re_specs + resp_term = NULL are fully supported by the multi path.
##
## Env (BC canonical):
##   R_E_TECH=1 R_K_WITHIN=1 R_DROP_SPARSE_K=1 R_NEFF_MIN=30 \
##   R_DATA_INFORMED_TAU=1 R_EDGE_CORRECT=1 R_BLEED_PERCELL=1 \
##   R_DISP_MODEL=nb2 R_N_ITER=32 R_PQL_THREADS=4

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(forcats)
  library(SpatialExperiment); library(SingleCellExperiment); library(Matrix); library(BiocParallel)
  library(FNN); library(dbscan)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
source("scripts/zzz.R"); .load_all("scripts")
source("scripts/helpers/pace_mvpql.R")
source("scripts/helpers/pace_mvpql_joint_multi.R")

ts  <- function() format(Sys.time(),"%H:%M:%S")
say <- function(...) cat(sprintf("[%s] %s\n", ts(), sprintf(...)))

## ---- Isotropic area-fraction edge correction (copied from Mel redesign) ----
## area_fraction(i) = area(disc(i, r) ∩ image_rectangle) / (π r²)
## X_corrected = X_raw / area_fraction.  Angular quadrature, N=1000 directions.
##
## LINEAGE: this is the dependency-free re-implementation of the thesis-era
## `borderEdge()` in `lateral_spillover/scripts/model_scripts_lite/helpers/
## compute_abundance.R` (invoked there as `edge = "isotropic"`). That version
## computes the SAME quantity exactly via spatstat.geom (discs ∩ owin / πr²);
## here we approximate area(disc ∩ rectangle) by angular quadrature so the
## builder carries no spatstat dependency. Difference vs thesis: applied to the
## E^tech ambient field at r = 3·h_tech = 15 µm, not the 25 µm neighbourhood.
.compute_area_fraction <- function(coords, r,
                                    xmin = NULL, xmax = NULL,
                                    ymin = NULL, ymax = NULL,
                                    N_angles = 1000L) {
  if (is.null(xmin)) xmin <- min(coords[, 1])
  if (is.null(xmax)) xmax <- max(coords[, 1])
  if (is.null(ymin)) ymin <- min(coords[, 2])
  if (is.null(ymax)) ymax <- max(coords[, 2])
  theta <- seq(0, 2*pi, length.out = N_angles + 1L)[-1L]
  cos_t <- cos(theta); sin_t <- sin(theta)
  n <- nrow(coords)
  af <- numeric(n)
  norm_factor <- pi * r^2
  for (i in seq_len(n)) {
    x0 <- coords[i, 1]; y0 <- coords[i, 2]
    d_right  <- ifelse(cos_t > 0, (xmax - x0) / cos_t, Inf)
    d_left   <- ifelse(cos_t < 0, (xmin - x0) / cos_t, Inf)
    d_top    <- ifelse(sin_t > 0, (ymax - y0) / sin_t, Inf)
    d_bottom <- ifelse(sin_t < 0, (ymin - y0) / sin_t, Inf)
    d_max  <- pmin(d_right, d_left, d_top, d_bottom)
    r_eff  <- pmin(r, d_max)
    af[i] <- (pi * sum(r_eff^2) / N_angles) / norm_factor
  }
  af
}

t0 <- Sys.time()
TYPES <- c("T_Cell","Stromal","Macrophage","Tumour","Endothelial",
           "B_Cell","Myoepithelial","Mast","Dendritic_Cell")

H_TECH <- as.numeric(Sys.getenv("R_H_TECH", unset = "5"))
H_BIO  <- as.numeric(Sys.getenv("R_H_BIO",  unset = "30"))
OUT    <- Sys.getenv("R_OUT",
  unset = "data/breast_cancer/sweeps/mvpql_percell_hc.rds")
KNN_K  <- 50
say("BC per-cell-HC redesign: h_tech=%.1f µm, h_bio=%.1f µm  out=%s", H_TECH, H_BIO, OUT)

say("[1/7] Loading cached Y + df ...")
YDF <- Sys.getenv("R_YDF", unset = "data/breast_cancer/Y_df_for_mcsd.rds")
say("    input bundle: %s", YDF)
ymv <- readRDS(YDF)
Y <- ymv$Y; df <- ymv$df
df$celltype <- factor(as.character(df$cellType), levels = TYPES)
## image grouping = slide (17 FOVs of one biopsy).  Used by: E^tech per-slide,
## rank-1 λ_m, and K_WITHIN centering.  NO imageID RE block is added.
df$imageID <- factor(as.character(df$slide))
DET_MIN <- as.numeric(Sys.getenv("R_DET_MIN", unset = "0.05"))
say("    using max-per-celltype detection filter >= %.2f", DET_MIN)
if (DET_MIN > 0.05 + 1e-6) {
  ct_chr <- as.character(df$celltype)
  ct_idx_lst <- lapply(TYPES, function(c) which(ct_chr == c))
  max_det <- vapply(seq_len(ncol(Y)), function(gi) {
    max(vapply(ct_idx_lst, function(idx)
      if (length(idx) == 0L) 0 else mean(Y[idx, gi] > 0), numeric(1)))
  }, numeric(1))
  keep_g <- which(max_det >= DET_MIN)
  say("    %d / %d genes pass max-CT >= %.2f (dropping %d)",
      length(keep_g), ncol(Y), DET_MIN, ncol(Y) - length(keep_g))
  Y <- Y[, keep_g, drop = FALSE]
}
Y <- as.matrix(Y)
df$nCount <- rowSums(Y)
say("    Y: %d cells x %d genes; %d slides; coords x=[%.0f,%.0f] y=[%.0f,%.0f]",
    nrow(Y), ncol(Y), nlevels(df$imageID), min(df$x), max(df$x), min(df$y), max(df$y))

## ---- 2. Kernel-weighted neighbour predictors (radius frNN, uncapped) ----
## FIX (audit 2026-06-15): the locked recipe used FNN::get.knn(k=50), which caps
## each cell at its 50 nearest neighbours -> in dense tissue the Gaussian K_bio
## (exp(-(d/h_bio)^2), still ~0.2 weight at the 50th NN) is truncated, undercounting
## neighbour density density-dependently. Replace with dbscan::frNN(eps=3*h_bio):
## every cell within RAD contributes, no KNN truncation. RAD=3*h_bio=90um: both
## kernels (exp(-(d/h_bio)^2) and exp(-d/h_tech)) are <1e-3 beyond this. Self is
## excluded (frNN, like get.knn, omits the point itself).
RAD <- 3 * H_BIO
say("[2/7] Building kernel-weighted neighbour predictors (frNN, eps=%d µm, uncapped) ...", RAD)
say("    h_tech=%d µm (exp(-d/h_tech))   h_bio=%d µm (exp(-d^2/h_bio^2))", H_TECH, H_BIO)
## GLOBAL frNN is CORRECT here: BC is ONE physical section and "slide" is a (spatial) CV-fold
## label, NOT a separate FOV. Cells of different slides are spatially interleaved in a single
## coordinate frame (verified 2026-06-18: a slide-1 cell's nearest slide-2 cell is ~8-20um
## away, 0% within 2um -> no coincidences -> not overlapping sections). Cross-slide neighbours
## are GENUINE physical neighbours, so they must be included. (Mel differs: imageID = patient,
## genuinely separate samples; its builder runs frNN PER image.)
coords <- as.matrix(df[, c("x", "y")])
n <- nrow(df)
ct <- as.character(df$celltype)
fr <- dbscan::frNN(coords, eps = RAD)
i_vec <- rep.int(seq_len(n), lengths(fr$id))
j_vec <- unlist(fr$id,   use.names = FALSE)
d_vec <- unlist(fr$dist, use.names = FALSE)
jc    <- match(ct[j_vec], TYPES)
ok    <- !is.na(jc)
i_ok  <- i_vec[ok]; d_ok <- d_vec[ok]; jc_ok <- jc[ok]
say("    frNN: %.1f neighbours/cell median; %d kernel contributions",
    stats::median(lengths(fr$id)), length(i_ok))
## sparseMatrix SUMS duplicate (i, celltype) entries -> correct kernel accumulation
K_tech <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok / H_TECH),
                                         dims = c(n, length(TYPES))))
K_bio  <- as.matrix(Matrix::sparseMatrix(i = i_ok, j = jc_ok, x = exp(-d_ok^2 / H_BIO^2),
                                         dims = c(n, length(TYPES))))
colnames(K_tech) <- paste0(TYPES, "_near"); colnames(K_bio) <- TYPES
rm(fr, i_vec, j_vec, d_vec, jc, ok, i_ok, d_ok, jc_ok)
say("    K_tech mean range=[%.3f, %.3f]; K_bio mean range=[%.3f, %.3f]",
    min(colMeans(K_tech)), max(colMeans(K_tech)),
    min(colMeans(K_bio)),  max(colMeans(K_bio)))

## Inject K_bio columns into df under names TYPES so re_specs can find them
for (tc in TYPES) df[[tc]] <- K_bio[, tc]

## ---- 3. R_DROP_SPARSE_K + R_K_WITHIN (from Mel redesign) ----
DROP_SPARSE_K <- nzchar(Sys.getenv("R_DROP_SPARSE_K"))
if (DROP_SPARSE_K) {
  neff_min <- as.numeric(Sys.getenv("R_NEFF_MIN", unset = "50"))
  say("[3/7] R_DROP_SPARSE_K: dropping K_<focal,nb> with n_eff < %.0f", neff_min)
  dropped <- character()
  for (f in TYPES) {
    cells_f <- which(df$celltype == f)
    if (!length(cells_f)) next
    for (nb in TYPES) {
      K_vals <- df[[nb]][cells_f]
      Kc <- K_vals - mean(K_vals, na.rm = TRUE)
      max_sq <- max(Kc^2, na.rm = TRUE)
      n_eff <- if (max_sq > 0) sum(Kc^2, na.rm = TRUE) / max_sq else 0
      if (n_eff < neff_min) {
        df[[nb]][cells_f] <- 0
        dropped <- c(dropped, sprintf("%s::%s(n_eff=%.1f)", f, nb, n_eff))
      }
    }
  }
  say("    dropped %d (focal, nb) pairs: %s", length(dropped),
      paste(dropped, collapse = ", "))
}

K_WITHIN <- nzchar(Sys.getenv("R_K_WITHIN"))
if (K_WITHIN) {
  say("    R_K_WITHIN: replacing K_<nb> with within-(slide, celltype) centered")
  for (tc in TYPES) {
    df[[tc]] <- df[[tc]] - ave(df[[tc]], df$imageID, df$celltype, FUN = mean)
  }
  max_abs_mean <- max(vapply(TYPES, function(tc) {
    max(abs(ave(df[[tc]], df$imageID, df$celltype, FUN = mean)))
  }, numeric(1)))
  say("    K_<nb>_within max |within-group mean| = %.2e", max_abs_mean)
}

## ---- 4. X_fixed (intercept-only under E^tech; no Responder for BC) ----
offset_vec <- log(df$nCount)
E_TECH <- nzchar(Sys.getenv("R_E_TECH"))
if (E_TECH) {
  X_fixed <- matrix(1, nrow = nrow(df), ncol = 1, dimnames = list(NULL, "(Intercept)"))
  say("[4/7] X_fixed: intercept-only (E^tech replaces _near columns)")
} else {
  X_fixed <- cbind(`(Intercept)` = 1, K_tech)
  say("[4/7] X_fixed: %d x %d (intercept + tech-kernel cols)", nrow(X_fixed), ncol(X_fixed))
}

## re_specs: celltype block ONLY (no Responder interaction, no imageID block).
re_specs <- list(
  list(group_col = "celltype",
       formula   = stats::as.formula(paste0("~ 1 + ", paste(TYPES, collapse = " + ")))))
say("    re_specs: single celltype block ~ 1 + (%s)", paste(TYPES, collapse = " + "))

## ---- 4b. E^tech + log1p + edge correction + rank-1 args ----
pre_offset_mat <- NULL
ambient_E_tech_mat <- NULL
ambient_image_idx      <- NULL
ambient_n_images       <- 0L
ORTHO   <- nzchar(Sys.getenv("R_BLEED_ORTHO"))
PERCELL <- nzchar(Sys.getenv("R_BLEED_PERCELL"))
percell_anchor_mask <- NULL
if (ORTHO || PERCELL) {
  ## Per-cell-type mean expression (clean-ish identity profile), used for ORTHO neighbour
  ## centering and/or PERCELL negative-control anchor definition.
  type_means <- t(vapply(TYPES, function(tt) colMeans(Y[df$celltype == tt, , drop = FALSE]),
                         numeric(ncol(Y))))
  rownames(type_means) <- TYPES; colnames(type_means) <- colnames(Y)
  say("    built type-mean matrix (%d x %d)", nrow(type_means), ncol(type_means))
}
if (PERCELL) {
  ## HOMOTYPIC-CORE negative-control anchors (de-circularized): instead of raw type-means
  ## (computed from contaminated cells), define each type's baseline from SPATIALLY-ISOLATED
  ## cells (>=50% same-type neighbours within 30um) -> least-contaminated identity profile.
  ## A gene is a negative control for X if another type owns it AND X's CORE level < 10% of owner.
  HC_HOMO <- as.numeric(Sys.getenv("R_PERCELL_HOMOFRAC", unset = "0.5"))
  same_frac <- rep(NA_real_, nrow(df))
  for (s in unique(df$slide)) {
    si <- which(df$slide == s); if (length(si) < 50) next
    ctl <- df$celltype[si]; fr <- dbscan::frNN(as.matrix(df[si, c("x","y")]), eps = 30)
    for (k in seq_along(si)) { idk <- fr$id[[k]]; same_frac[si[k]] <- if (length(idk)) mean(ctl[idk] == ctl[k]) else 1 }
  }
  core <- which(same_frac >= HC_HOMO)
  core_means <- type_means                                   # raw fallback for sparse types
  ncore <- integer(0)
  for (X in TYPES) {
    idx <- core[df$celltype[core] == X]; ncore <- c(ncore, length(idx))
    if (length(idx) >= 20) core_means[X, ] <- colMeans(Y[idx, , drop = FALSE])
  }
  say("    PERCELL homotypic-core cells (>=%.0f%% same-type): %d/%d; per-type median=%d (raw fallback if <20)",
      100*HC_HOMO, length(core), nrow(df), as.integer(median(ncore)))
  owner_mean <- apply(core_means, 2, max)
  owner_t    <- TYPES[apply(core_means, 2, which.max)]
  percell_anchor_mask <- matrix(0, nrow(df), ncol(Y))
  ngc <- integer(0)
  for (X in TYPES) {
    is_anc <- owner_t != X & owner_mean > 0.1 & (core_means[X, ] / pmax(owner_mean, 1e-9) < 0.1)
    rr <- which(df$celltype == X)
    if (length(rr)) percell_anchor_mask[rr, ] <- matrix(as.numeric(is_anc), length(rr), ncol(Y), byrow = TRUE)
    ngc <- c(ngc, sum(is_anc))
  }
  say("    R_BLEED_PERCELL (homotypic-core anchors): genes/type median=%d range=[%d,%d]", as.integer(median(ngc)), min(ngc), max(ngc))
}
if (E_TECH) {
  say("[4b] Building E^tech (cross-celltype expression-weighted spillover) per slide ...")
  t_et <- Sys.time()
  rad  <- 3 * H_TECH
  E_tech <- matrix(0, nrow(df), ncol(Y), dimnames = list(NULL, colnames(Y)))
  by_slide <- split(seq_len(nrow(df)), df$imageID)
  ct_all <- as.character(df$celltype)
  for (si in seq_along(by_slide)) {
    rows <- by_slide[[si]]
    if (length(rows) < 2) next
    coords_im <- as.matrix(df[rows, c("x","y")])
    nn_e <- dbscan::frNN(coords_im, eps = rad)
    Y_im <- Y[rows, , drop = FALSE]
    ct_im <- ct_all[rows]
    for (i in seq_along(rows)) {
      nbr_local <- nn_e$id[[i]]
      if (!length(nbr_local)) next
      d_local   <- nn_e$dist[[i]]
      keep <- ct_im[nbr_local] != ct_im[i]
      if (!any(keep)) next
      w <- exp(-d_local[keep] / H_TECH)
      Yblk <- Y_im[nbr_local[keep], , drop = FALSE]
      if (ORTHO) Yblk <- Yblk - type_means[ct_im[nbr_local[keep]], , drop = FALSE]  # center per neighbour type
      E_tech[rows[i], ] <- as.numeric(crossprod(w, Yblk))
    }
    if (si %% 5 == 0) cat(sprintf("    [E^tech] slide %d/%d done\n", si, length(by_slide)))
  }
  say("    E^tech computed in %.1f s (range = [%.2f, %.2f])",
      as.numeric(difftime(Sys.time(), t_et, units = "secs")),
      min(E_tech), max(E_tech))

  ## R_EDGE_CORRECT: isotropic area-fraction correction at r = rad = 3*H_TECH.
  if (nzchar(Sys.getenv("R_EDGE_CORRECT"))) {
    af_etech <- numeric(nrow(df))
    for (si in seq_along(by_slide)) {
      rows <- by_slide[[si]]
      if (!length(rows)) next
      coords_im <- as.matrix(df[rows, c("x","y")])
      af_etech[rows] <- .compute_area_fraction(
        coords_im, rad,
        min(coords_im[,1]), max(coords_im[,1]),
        min(coords_im[,2]), max(coords_im[,2]))
    }
    say("    R_EDGE_CORRECT (E^tech): af median=%.3f q10=%.3f q90=%.3f (%.1f%% cells af<0.95)",
        median(af_etech), quantile(af_etech, 0.10), quantile(af_etech, 0.90),
        100 * mean(af_etech < 0.95))
    E_tech <- sweep(E_tech, 1, af_etech, "/")
    say("    Edge-corrected E^tech range = [%.2f, %.2f]", min(E_tech), max(E_tech))
  }

  if (ORTHO) {
    ## eps can be negative -> use directly (NOT log1p). This is the orthogonalized covariate.
    eps_mat <- E_tech; storage.mode(eps_mat) <- "double"
    say("    eps (density-orthogonalized spillover covariate) range = [%.2f, %.2f]", min(eps_mat), max(eps_mat))
  } else {
    log_E_tech <- log1p(E_tech)
    storage.mode(log_E_tech) <- "double"
  }

  if (nzchar(Sys.getenv("R_BLEED_PERCELL"))) {
    ## Per-cell contamination loading: mu = mu_bio + rho_i*a_ig (a = RAW E^tech ambient).
    say("    R_BLEED_PERCELL=1: per-cell shrunken contamination loading rho_i on local ambient (RAW E^tech).")
    ambient_E_tech_mat <- E_tech
    storage.mode(ambient_E_tech_mat) <- "double"
    ambient_image_idx      <- as.integer(df$imageID)
    ambient_n_images       <- nlevels(df$imageID)
    pre_offset_mat <- NULL
  } else if (nzchar(Sys.getenv("R_NO_BLEED"))) {
    ## No in-model spillover term: spillover is handled EXTERNALLY by anchor-based
    ## decontamination of the counts (R_YDF points to the decontaminated bundle).
    ## X_fixed stays intercept-only (E_TECH branch), so the model = intercept +
    ## spatial random slopes + residual, with NO technical/bleed component.
    say("    R_NO_BLEED=1: NO in-model spillover term (decontamination handles spillover externally).")
    ambient_E_tech_mat <- NULL
    ambient_n_images       <- 0L
    pre_offset_mat <- NULL
  } else {
    stop("bc_fit_joint_kernel_redesign.R: set R_BLEED_PERCELL=1 or R_NO_BLEED=1")
  }
}

## ---- 5. Joint multi-block fit ----
NO_BLEED <- nzchar(Sys.getenv("R_NO_BLEED"))
PQL_THREADS <- as.integer(Sys.getenv("R_PQL_THREADS", unset = "4"))
say("[5/7] fit_pace_mvpql_joint_multi (h_tech=%.1f, h_bio=%.1f, %d C++ threads) %s...",
    H_TECH, H_BIO, PQL_THREADS, if (ambient_n_images > 0L) "(per-cell-HC bleed)" else "")

data_informed_W <- NULL
if (nzchar(Sys.getenv("R_DATA_INFORMED_TAU"))) {
  say("    R_DATA_INFORMED_TAU: detection x normalised-K-variance weights ...")
  re_tmp <- build_random_design_multi(df, re_specs)
  data_informed_W <- .compute_data_informed_weights(
    re = re_tmp, Y = Y, df = df,
    focals = re_tmp$blocks[[1]]$group_levels, TYPES = TYPES,
    celltype_col = "celltype", verbose = TRUE)
  rm(re_tmp)
}

t_fit <- Sys.time()
fit <- fit_pace_mvpql_joint_multi(
  Y = Y, X_fixed = X_fixed, df = df, re_specs = re_specs,
  offset_vec = offset_vec,
  pre_offset_mat = pre_offset_mat,
  data_informed_W = data_informed_W,
  ambient_mat        = ambient_E_tech_mat,
  ambient_image_idx  = ambient_image_idx,
  ambient_n_images   = ambient_n_images,
  bleed_percell          = nzchar(Sys.getenv("R_BLEED_PERCELL")),
  percell_anchor_mask    = percell_anchor_mask,
  n_iter = as.integer(Sys.getenv("R_N_ITER", unset = "16")), tol = 5e-3,
  tau_shrinkage = Sys.getenv("R_TAU", unset = "adaptive"),
  BPPARAM = SerialParam(), n_threads = PQL_THREADS,
  interior_precision = 1L, chunk_size = 128L, verbose = TRUE)
say("    fit done in %.1f min", as.numeric(difftime(Sys.time(), t_fit, units = "mins")))

## ---- 6. mashr + decomposition (resp_term = NULL: BC has no Responder) ----
say("[6/7] mashr + decomposition ...")
results_mv <- mvpql_to_results_multi(fit, keep_block = "celltype")
shrunken_long <- apply_mashr_shrinkage(
  results = results_mv, focals = TYPES, neighbours = TYPES, resp_term = NULL)
say("    %s shrunk rows", format(nrow(shrunken_long), big.mark = ","))

dec <- mvpql_variance_decomposition_multi(
  fit = fit, df = df, Y = Y, vars = TYPES, X_fixed = X_fixed,
  resp_term = NULL, focal_levels = TYPES)

## Backward-compat alias: BC figure scripts read decomposition$gene_focal_4block.
## The multi decomp emits gene_focal_5block (V_state_responder = 0 for BC).
## Collapse the 5-block per-gene table to BC's 4-block % schema.
gene_focal_5block <- dec$gene_focal_5block
g5 <- gene_focal_5block
Total4 <- with(g5, celltype_offset_sq + V_state_baseline + V_state_responder +
                    V_spill + V_disp)
dec$gene_focal_4block <- tibble::tibble(
  focal             = g5$focal,
  gene              = g5$gene,
  `Cell type %`     = 100 * g5$celltype_offset_sq / pmax(Total4, 1e-12),
  `Spatial state %` = 100 * (g5$V_state_baseline + g5$V_state_responder) / pmax(Total4, 1e-12),
  `Spillover %`     = 100 * g5$V_spill / pmax(Total4, 1e-12),
  `Residual %`      = 100 * g5$V_disp  / pmax(Total4, 1e-12),
  ## raw latent-scale variance components (carried so pooled/spec-weighted
  ## aggregates can be recomputed; V_state folds the zero responder block).
  celltype_offset_sq = g5$celltype_offset_sq,
  V_state            = g5$V_state_baseline + g5$V_state_responder,
  V_spill            = g5$V_spill,
  V_disp             = g5$V_disp,
  Total              = Total4,
  spec              = g5$spec,
  focal_mean        = g5$focal_mean)

## ---- 6b. SENSITIVITY (gated, default OFF): SVCA-style intrinsic cell-state confound check ----
## NOT part of the canonical decomposition. The canonical "Spatial cell state" block is an
## associational variance share (an UPPER BOUND); this sidecar quantifies how much of that spatial
## signal is SHARED with the cell's own intrinsic state (state-spatial commonality, leave-gene-out
## top-5 PCs). Written to a SEPARATE sidecar so the canonical object is untouched. Diagnostic only
## (scale-mixed + commonality-based + co-regulation confounded). See
## project_pace_intrinsic_state_2026-06-19.
if (Sys.getenv("R_INTRINSIC_STATE", unset = "0") == "1") {
  say("[6b] SENSITIVITY R_INTRINSIC_STATE: intrinsic cell-state confound (top-5 LOO PCs) ...")
  source("scripts/helpers/pace_intrinsic_state.R")
  intrinsic_sens <- pace_intrinsic_state_decomp(
    fit = fit, Y = Y, celltype = df$celltype, nCount = df$nCount,
    focal_levels = TYPES, n_pc = 5L)
  sens_out <- sub("\\.rds$", "_intrinsic_sensitivity.rds", OUT)
  saveRDS(intrinsic_sens, sens_out)
  say("    wrote SENSITIVITY sidecar %s (%d rows)", sens_out, nrow(intrinsic_sens))
}

## ---- 7. MCSD + R^2_S (BC pairs; term == nc; lfsr<0.05) ----
say("[7/7] MCSD + R^2_S ...")
PAIRS_BC <- list(c("Stromal","Tumour"), c("Macrophage","Tumour"),
                 c("Endothelial","Tumour"), c("Myoepithelial","Tumour"))

Z_re <- fit$re_meta$Z
cells_by_ct <- lapply(TYPES, function(c) which(Z_re[, paste0(c, "::(Intercept)")] != 0))
names(cells_by_ct) <- TYPES
alpha_g <- pmax(fit$alpha, 0)

mcsd_canonical <- list()
for (p in PAIRS_BC) {
  fc <- p[1]; nc <- p[2]; pk <- paste(fc, nc, sep="_")
  s <- shrunken_long |> dplyr::filter(focal == fc, neighbour == nc, term == nc) |>
    dplyr::distinct(gene, .keep_all = TRUE)
  fm <- gene_focal_5block |> dplyr::filter(focal == fc) |>
    dplyr::select(gene, spec, focal_mean) |> dplyr::distinct(gene, .keep_all = TRUE)
  if (!nrow(s) || !nrow(fm)) { cat(sprintf("  %s: skip\n", pk)); next }

  cells_c <- cells_by_ct[[fc]]
  col_N <- paste0(fc, "::", nc)
  if (!col_N %in% colnames(Z_re)) {
    cat(sprintf("  %s: pair column %s absent (dropped by n_eff) -> skip\n", pk, col_N))
    mcsd_canonical[[pk]] <- list(scores = s[0, ], status = "dropped (n_eff)")
    next
  }
  N_t <- as.numeric(Z_re[cells_c, col_N])
  var_N <- stats::var(N_t, na.rm = TRUE)
  mu_bar_per_gene <- Matrix::colMeans(fit$mu[cells_c, , drop = FALSE])
  gene_names_fit <- colnames(fit$U)
  if (is.null(gene_names_fit)) gene_names_fit <- colnames(Y)
  names(alpha_g) <- gene_names_fit
  names(mu_bar_per_gene) <- gene_names_fit

  res_all <- s |> dplyr::inner_join(fm, by = "gene") |>
    dplyr::mutate(
      MCSD     = (estimate_shrunk^2) * (spec^2) * pmax(focal_mean, 0),
      MCSD4    = (estimate_shrunk^2) * (spec^4) * pmax(focal_mean, 0),
      mu_bar   = mu_bar_per_gene[gene],
      alpha    = alpha_g[gene],
      V_S      = estimate_shrunk^2 * var_N,
      V_resid  = log(1 + (1 + pmax(alpha, 0)) / pmax(mu_bar, 1e-6)),
      V_total  = V_S + V_resid,
      R2_S     = V_S / pmax(V_total, 1e-12)
    )
  res <- res_all |>
    dplyr::filter(lfsr < 0.05) |>
    dplyr::arrange(dplyr::desc(MCSD)) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::rename(b_clean = estimate_shrunk) |>
    dplyr::select(rank, gene, MCSD, MCSD4, b_clean, spec, focal_mean,
                  R2_S, mu_bar, alpha, V_S, V_resid, V_total, lfsr, sd_shrunk)
  status <- if (nrow(res) >= 3) "significant" else "honestly null"
  mcsd_canonical[[pk]] <- list(scores = res, status = status)
  cat(sprintf("  %s n_sig=%d top-10 (%s): %s\n",
              pk, nrow(res), status, paste(head(res$gene, 10), collapse=", ")))
}

cell_meta <- data.frame(
  cell_id  = as.character(df$cell_id),
  nCount   = df$nCount,
  slide    = as.character(df$slide),
  celltype = as.character(df$celltype),
  stringsAsFactors = FALSE)
saveRDS(list(fit = fit, shrunken_long = shrunken_long, gene_set = colnames(Y),
              decomposition = dec, mcsd_canonical = mcsd_canonical,
              h_tech = H_TECH, h_bio = H_BIO, K_tech = K_tech, K_bio = K_bio,
              cell_meta = cell_meta),
        OUT)
say("Saved %s", OUT)
say("TOTAL: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

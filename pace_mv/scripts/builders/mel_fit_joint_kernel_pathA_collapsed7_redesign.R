## mel_fit_joint_kernel_pathA_collapsed7.R — REFIT 2026-05-29.
##
## Same recipe as mel_fit_joint_kernel_pathA_expanded.R EXCEPT the celltype
## label scheme is collapsed from 12 to 7 (Tumor_a..f -> single "Tumour").
##
## Motivation: post-hoc spec-rescue test (mel_collapse_tumour_posthoc.rds)
## showed SPP1 in Macrophage jumps from spec=0.36 (fails canonical 0.5 filter)
## to spec=0.76 when Tumor_a..f are pooled. The 12-CT fragmentation was
## under-powering both the spec metric AND the per-tumour-subtype mashr slopes
## for non-tumour focals' RxS signal. This refit collapses tumour subtypes so
## the spec rescue is structural (not post-hoc weighted), and the (Mac,Tumour)
## slope is estimated from the full pooled tumour neighbour kernel rather than
## six fragmented sub-slopes.
##
## TYPES collapse:
##   Tumor_a, Tumor_b, Tumor_c, Tumor_d, Tumor_e, Tumor_f  ->  Tumour
##   Endothelial, Fibroblast, Macrophage, B_Cell,
##   T_CD8_memory, Treg                                    ->  kept as-is
##
## Output: mvpql_kernel_h05_h30_pathA_collapsed7_patdz_adapt.rds

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(forcats); library(tidyr)
  library(SingleCellExperiment); library(Matrix); library(BiocParallel)
  library(zellkonverter); library(scran); library(FNN); library(dbscan)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
source("scripts/zzz.R"); .load_all("scripts")
source("scripts/helpers/pace_mvpql.R")
source("scripts/helpers/pace_mvpql_joint_multi.R")

ts  <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...) cat(sprintf("[%s] %s\n", ts(), sprintf(...)))

t0 <- Sys.time()
H_TECH <- as.numeric(Sys.getenv("R_H_TECH", unset = "5"))
H_BIO  <- as.numeric(Sys.getenv("R_H_BIO",  unset = "30"))
KNN_K  <- 50L

## ---- Isotropic area-fraction edge correction (added 2026-06-09) ----
## Implements the manuscript spec:
##     area_fraction(i) = area(disc(i, r) ∩ image_rectangle) / (π r²)
##     X_corrected = X_raw / area_fraction
##
## Computes the disc ∩ rectangle area by angular integration:
##     area(disc ∩ rect) = (1/2) ∫_0^{2π} min(r, r_max(θ))² dθ
## where r_max(θ) is the distance from the cell to the rectangle boundary in
## direction θ.  Discretised at N_angles = 1000 directions (standard angular
## quadrature; error ~ O(1/N) so ~0.1% accuracy).  Exact for N → ∞.  No floor;
## the integral is naturally > 0 for any cell inside the rectangle (>= 0.25
## at exact corners).
##
## LINEAGE: dependency-free re-implementation of the thesis-era `borderEdge()`
## in `lateral_spillover/scripts/model_scripts_lite/helpers/compute_abundance.R`
## (invoked there as `edge = "isotropic"`), which computes the same
## area(disc ∩ owin)/πr² exactly via spatstat.geom. Here applied to the E^tech
## ambient field at r = 3·h_tech = 15 µm (thesis applied it to the 25 µm
## neighbourhood-abundance / spillover counts).
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
  norm_factor <- pi * r^2  # full disc area (denominator)
  for (i in seq_len(n)) {
    x0 <- coords[i, 1]; y0 <- coords[i, 2]
    d_right  <- ifelse(cos_t > 0, (xmax - x0) / cos_t, Inf)
    d_left   <- ifelse(cos_t < 0, (xmin - x0) / cos_t, Inf)
    d_top    <- ifelse(sin_t > 0, (ymax - y0) / sin_t, Inf)
    d_bottom <- ifelse(sin_t < 0, (ymin - y0) / sin_t, Inf)
    d_max  <- pmin(d_right, d_left, d_top, d_bottom)
    r_eff  <- pmin(r, d_max)
    ## area = (1/2) × (2π/N) × Σ r_eff² = π × Σ r_eff² / N
    af[i] <- (pi * sum(r_eff^2) / N_angles) / norm_factor
  }
  af
}
OUT    <- Sys.getenv("R_OUT",
  unset = "data/simvi_melanoma/sweeps/mvpql_percell_hc.rds")
say("Mel COLLAPSED-7CT refit: h_tech=%.1f µm, h_bio=%.1f µm  out=%s", H_TECH, H_BIO, OUT)

## COLLAPSED 7-CT scheme: pool Tumor_a..f into single "Tumour" focal.
TYPES <- c("Tumour",
           "Endothelial", "Fibroblast", "Macrophage",
           "B_Cell", "T_CD8_memory", "Treg")
TUM_RAW <- c("Tumor_a","Tumor_b","Tumor_c","Tumor_d","Tumor_e","Tumor_f")

## ---- 1. Load + clean SCE (mirror expanded) ----
say("[1/7] Loading Mel SCE ...")
sce <- readH5AD(file = "data/simvi_melanoma/Melanoma_5612.h5ad", reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
px_size <- 0.12028
set.seed(1994)
## Response coding (R_MEL_RESP):
##   pd_nonpd (CANONICAL, 2026-06-16): keep ALL (PD/SD/CR/PR); PD vs non-PD(SD+CR+PR).
##     All 26 images. Matches SIMVI's PD-vs-nonPD contrast; PD stays the case group.
##   pd_resp: drop SD; PD vs RESP(CR+PR). 18 images. (Pre-2026-06-16 canonical; now a sensitivity.)
##   binary_sd: keep SD; RECIST binary NonResp(PD+SD) vs Resp(PR+CR). 26 images. (Sensitivity.)
MEL_RESP <- Sys.getenv("R_MEL_RESP", unset = "pd_nonpd")
if (identical(MEL_RESP, "binary_sd")) {
  say("    R_MEL_RESP=binary_sd: keeping SD; NonResp=PD+SD vs Resp=PR+CR")
  sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR"))]
  sce$BEST_RESPONSE_BY_SCAN <- ifelse(
    as.character(sce$BEST_RESPONSE_BY_SCAN) %in% c("CR","PR"), "Resp", "NonResp")
  RESP_REF <- "Resp"; RESP_CASE <- "NonResp"
} else if (identical(MEL_RESP, "pd_resp")) {
  ## Path A: drop SD, keep CR+PR (collapsed to "RESP") and PD.
  say("    R_MEL_RESP=pd_resp: dropping SD; PD vs RESP(CR+PR)")
  sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "CR", "PR"))]
  sce$BEST_RESPONSE_BY_SCAN <- ifelse(
    as.character(sce$BEST_RESPONSE_BY_SCAN) %in% c("CR","PR"),
    "RESP", as.character(sce$BEST_RESPONSE_BY_SCAN))
  RESP_REF <- "RESP"; RESP_CASE <- "PD"
} else {
  ## CANONICAL pd_nonpd: keep ALL valid responses; PD vs non-PD(SD+CR+PR).
  say("    R_MEL_RESP=pd_nonpd (CANONICAL): keeping ALL; PD vs non-PD(SD+CR+PR)")
  sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR"))]
  sce$BEST_RESPONSE_BY_SCAN <- ifelse(
    as.character(sce$BEST_RESPONSE_BY_SCAN) == "PD", "PD", "nonPD")
  RESP_REF <- "nonPD"; RESP_CASE <- "PD"
}
RESP_TERM <- paste0("Responder", RESP_CASE)   # ResponderPD (pd_nonpd / pd_resp) or ResponderNonResp
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
  'B.cell'        = "B_Cell",
  'T.CD8.memory'  = "T_CD8_memory",
  'endothelial'   = "Endothelial",
  'fibroblast'    = "Fibroblast",
  'macrophage'    = "Macrophage")

## *** COLLAPSE STEP: Tumor_a..f -> "Tumour" ***
n_tum_pre <- sum(df_raw$celltype %in% TUM_RAW)
df_raw$celltype <- ifelse(df_raw$celltype %in% TUM_RAW, "Tumour", df_raw$celltype)
say("    collapsed Tumor_a..f -> Tumour (%s cells re-labelled)",
    format(n_tum_pre, big.mark=","))

df_raw <- df_raw |> dplyr::filter(celltype %in% TYPES)
df_raw$imageID <- as.character(df_raw$imageID)
df_raw$cell_id <- as.character(df_raw$cellID_str)
hvg_genes <- make.names(rownames(sce))
say("    %s cells, %d genes (full panel)", format(nrow(df_raw), big.mark=","), length(hvg_genes))
cat("    cell counts per collapsed celltype:\n")
print(table(df_raw$celltype))

## DenoIST hook (kept for compatibility; not used by default here)
.dn <- Sys.getenv("R_DENOIST_COUNTS")
if (nzchar(.dn)) {
  say("    [R_DENOIST_COUNTS] overriding counts from %s", basename(.dn))
  dn <- readRDS(.dn)
  mi <- match(df_raw$cell_id, dn$cell_id)
  if (anyNA(mi)) stop("DenoIST bundle missing ", sum(is.na(mi)), " cell_ids")
  gi <- intersect(hvg_genes, dn$genes)
  Ycr <- dn$Y[mi, match(gi, dn$genes), drop = FALSE]
  df_raw[, gi] <- Ycr
}

## ---- 2. Build kernel-weighted neighbour predictors (per imageID, radius frNN) ----
## FIX (2026-06-15): switched from per-image FNN::get.knn(k=50) to per-image
## dbscan::frNN(eps = 3*h_bio = 90um). Mel CosMx FOVs are dense (median ~2200 cells/FOV,
## within-FOV 1-NN ~5.6um), so k=50 truncated the Gaussian K_bio at the 50th neighbour
## (~31um, weight ~0.34) — 99% of cells had their 50th NN inside the kernel's ~60um range.
## frNN captures every neighbour within RAD, no KNN truncation. STILL per-image (frNN run
## within each FOV separately — no cross-FOV mixing). Mirrors the BC builder's frNN kernel.
say("[2/7] Building kernel-weighted predictors per image (frNN, eps=%d um, uncapped) ...", as.integer(3*H_BIO))
df_raw <- df_raw[order(df_raw$imageID), ]
images <- unique(df_raw$imageID); nI <- length(images)
n <- nrow(df_raw)
RAD <- 3 * H_BIO
I_acc <- integer(0); J_acc <- integer(0); WT_acc <- numeric(0); WB_acc <- numeric(0)
for (im in images) {
  rows <- which(df_raw$imageID == im)
  if (length(rows) < 2L) next
  coords_im <- as.matrix(df_raw[rows, c("x","y")])
  ct_im <- as.character(df_raw$celltype[rows])
  fr <- dbscan::frNN(coords_im, eps = RAD)
  ll <- lengths(fr$id); if (!sum(ll)) next
  i_loc <- rep.int(seq_along(rows), ll)
  j_loc <- unlist(fr$id,   use.names = FALSE)
  d     <- unlist(fr$dist, use.names = FALSE)
  jc    <- match(ct_im[j_loc], TYPES)
  ok    <- !is.na(jc)
  I_acc  <- c(I_acc,  rows[i_loc[ok]]); J_acc  <- c(J_acc,  jc[ok])
  WT_acc <- c(WT_acc, exp(-d[ok] / H_TECH)); WB_acc <- c(WB_acc, exp(-d[ok]^2 / H_BIO^2))
}
## sparseMatrix SUMS duplicate (cell, celltype) entries -> correct kernel accumulation
K_tech <- as.matrix(Matrix::sparseMatrix(i = I_acc, j = J_acc, x = WT_acc, dims = c(n, length(TYPES))))
K_bio  <- as.matrix(Matrix::sparseMatrix(i = I_acc, j = J_acc, x = WB_acc, dims = c(n, length(TYPES))))
colnames(K_tech) <- paste0(TYPES, "_near"); colnames(K_bio) <- TYPES
rm(I_acc, J_acc, WT_acc, WB_acc)
say("    K_tech mean range=[%.3f, %.3f]; K_bio mean range=[%.3f, %.3f]",
    min(colMeans(K_tech)), max(colMeans(K_tech)),
    min(colMeans(K_bio)),  max(colMeans(K_bio)))

## R_EDGE_CORRECT for K_bio / K_tech is NOT implemented in this builder.
## The manuscript spec's edge correction is defined for a radius-based
## neighbourhood at radius r.  K_bio / K_tech here are KNN-based (k=50) with
## Gaussian / exponential decay — there is no single radius r to plug into
## area_fraction = (disc ∩ rectangle) / (π r²).  Edge correction applied only
## to E^tech below, where r = rad = 3 × H_TECH is unambiguous.

## ---- 3. Detection filter (det>=5% in any focal) and df assembly ----
say("[3/7] det>=5%% gene filter + assembling df ...")
hvg_in_df <- intersect(hvg_genes, names(df_raw))
ct_chr <- as.character(df_raw$celltype)
ct_idx_lst <- lapply(TYPES, function(c) which(ct_chr == c))
chunk_size <- 50L
max_det <- numeric(length(hvg_in_df)); names(max_det) <- hvg_in_df
for (s in seq.int(1L, length(hvg_in_df), by = chunk_size)) {
  e <- min(s + chunk_size - 1L, length(hvg_in_df))
  m_chunk <- as.matrix(df_raw[, hvg_in_df[s:e], drop = FALSE]) > 0
  for (gi in seq_len(e - s + 1L)) {
    g_idx <- s + gi - 1L
    max_det[g_idx] <- max(vapply(ct_idx_lst, function(idx)
      if (length(idx) == 0L) 0 else mean(m_chunk[idx, gi]), numeric(1)))
  }; rm(m_chunk)
}; invisible(gc(verbose = FALSE))
DET_MIN <- as.numeric(Sys.getenv("R_DET_MIN", unset = "0.05"))
say("    max-per-celltype detection filter >= %.2f", DET_MIN)
hvg_in_df <- hvg_in_df[max_det >= DET_MIN]
df_raw$nCount <- rowSums(as.matrix(df_raw[, hvg_in_df, drop = FALSE]))
df_raw <- df_raw[df_raw$nCount > 0, ]
df <- df_raw
df$celltype <- factor(as.character(df$celltype), levels = TYPES)
df$imageID  <- factor(as.character(df$imageID))
df$Responder <- forcats::fct_relevel(as.factor(df$BEST_RESPONSE_BY_SCAN), RESP_REF)
df$.resp_dummy <- as.integer(df$Responder == RESP_CASE)
say("    Responder split: %s",
    paste(names(table(df$Responder)), table(df$Responder), sep="=", collapse=", "))
n_pat <- df |> dplyr::distinct(imageID, Responder) |>
  dplyr::count(Responder)
say("    Per-imageID Responder count: %s",
    paste(n_pat$Responder, n_pat$n, sep="=", collapse=", "))

## ---- Patient-level confounders ----
df$Site <- ifelse(as.character(df$SPEC_CATEGORY) %in% c("MET", "NODE"),
                  "METASTATIC", "PRIMARY")
df$Site <- forcats::fct_relevel(as.factor(df$Site), "PRIMARY")
df$Treatment <- ifelse(as.character(df$TREATMENT) == "IPI+NIVO",
                       "COMBO", "MONO")
df$Treatment <- forcats::fct_relevel(as.factor(df$Treatment), "MONO")
df$Mutation <- ifelse(as.character(df$MUTATION) == "NONE", "NONE", "DRIVER")
df$Mutation <- forcats::fct_relevel(as.factor(df$Mutation), "NONE")
say("    Site: %s  Treatment: %s  Mutation: %s",
    paste(names(table(df$Site)), table(df$Site), sep="=", collapse=","),
    paste(names(table(df$Treatment)), table(df$Treatment), sep="=", collapse=","),
    paste(names(table(df$Mutation)), table(df$Mutation), sep="=", collapse=","))

for (tc in TYPES) df[[tc]] <- K_bio[match(df$cell_id, df_raw$cell_id), tc]

## R_DROP_SPARSE_K: structurally drop K_<focal,nb> from the random-slope
## formula when the (focal, neighbour) pair has too little K variation in
## cells of the focal celltype to support a reliable slope estimate.
##
## Mechanism: zero out df[[nb]] for cells of focal f when n_eff(f, nb) <
## R_NEFF_MIN.  This makes the Z column `<f>::<nb>` (and `<f>::ResponderPD:<nb>`)
## all-zero within the focal's cells -- PQL inner solve gets no data signal,
## BLUP collapses to zero.  Non-focal cells' K values are preserved.
##
## n_eff = sum(K - mean(K))^2 / max((K - mean(K))^2)
##       ≈ "how many cells contribute meaningfully to the slope"
DROP_SPARSE_K <- nzchar(Sys.getenv("R_DROP_SPARSE_K"))
if (DROP_SPARSE_K) {
  neff_min <- as.numeric(Sys.getenv("R_NEFF_MIN", unset = "50"))
  say("    R_DROP_SPARSE_K: dropping K_<focal,nb> with n_eff < %.0f", neff_min)
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

## R_K_WITHIN: replace each K_<neighbour> column with its WITHIN-(imageID,
## focal_celltype) centered version:
##   K_<nb>_within[i] = K_<nb>[i] - mean(K_<nb> | imageID(i), celltype(i))
## By construction, K_<nb>_within has zero mean within every (image, celltype)
## cell.  The downstream random slope on this predictor therefore captures
## ONLY within-image-within-celltype gradient -- the (Responder, K_<nb>)
## cohort-level confound is structurally eliminated at the design level.
## Between-image disease effects get absorbed by the Responder fixed effect
## and the imageID random intercept where they biologically belong, not by
## the spatial-gradient random slope.
K_WITHIN <- nzchar(Sys.getenv("R_K_WITHIN"))
if (K_WITHIN) {
  say("    R_K_WITHIN: replacing K_<nb> with within-(imageID, celltype) centered")
  for (tc in TYPES) {
    df[[tc]] <- df[[tc]] - ave(df[[tc]], df$imageID, df$celltype, FUN = mean)
  }
  ## Diagnostic: confirm zero within-group means (should be ~1e-15 numerical)
  max_abs_mean <- max(vapply(TYPES, function(tc) {
    max(abs(ave(df[[tc]], df$imageID, df$celltype, FUN = mean)))
  }, numeric(1)))
  say("    K_<nb>_within max |within-group mean| = %.2e", max_abs_mean)
}

## RE_TYPES = the neighbour names used in the random-effect formula.
RE_TYPES  <- TYPES

for (tc in TYPES) df[[paste0(tc, "_imgz")]] <- as.numeric(scale(df[[tc]]))
say("    df: %s cells, %d genes, %d images", format(nrow(df), big.mark=","),
    length(hvg_in_df), nlevels(df$imageID))

## ---- 4. X_fixed (intercept + Responder + per-CT _near + Site/Treatment/Mutation x TYPES) ----
C3 <- nzchar(Sys.getenv("R_C3"))
UNIFIED <- nzchar(Sys.getenv("R_UNIFIED"))
say("[4/7] Building X_fixed + re_specs %s...",
    if (UNIFIED) "(UNIFIED)" else if (C3) "(C3)" else "(default _near per-CT)")
conf_intxn_rhs <- paste(c(
  paste0("Site:",      TYPES),
  paste0("Treatment:", TYPES),
  paste0("Mutation:",  TYPES)
), collapse = " + ")
conf_main_rhs <- "Site + Treatment + Mutation"
if (UNIFIED) {
  fixed_rhs <- paste0("1 + Responder + ", conf_main_rhs, " + ", conf_intxn_rhs)
  spill_cols <- character(0)
} else if (C3) {
  df$neighbour_total_near <- rowSums(K_tech)[match(df$cell_id, df_raw$cell_id)]
  fixed_rhs <- paste0("1 + Responder + ", conf_main_rhs,
                       " + neighbour_total_near + ", conf_intxn_rhs)
  spill_cols <- "neighbour_total_near"
} else {
  spill_cols <- paste0(TYPES, "_near")
  for (sc in spill_cols) df[[sc]] <- K_tech[match(df$cell_id, df_raw$cell_id), sc]
  ## R_E_TECH: expression-weighted technical spillover (Option B = clean cross-celltype)
  ## E^tech[i, g] = Σ_{j != i, c(j) != c(i)} y_{j,g} × exp(-d_ij / h_tech)
  ## When enabled, the per-celltype `_near` fixed effects are REMOVED from X_fixed
  ## (β_spill,g handles spillover instead, as a per-gene fixed slope on E^tech).
  E_TECH <- nzchar(Sys.getenv("R_E_TECH"))
  PER_FOCAL_NEAR <- nzchar(Sys.getenv("R_PER_FOCAL_NEAR"))
  spill_part <- if (E_TECH) {
    say("    EXPRESSION-WEIGHTED TECH SPILLOVER (R_E_TECH): dropping _near columns")
    NULL  ## drop _near entirely
  } else if (PER_FOCAL_NEAR) {
    say("    PER-FOCAL β_near: using celltype:* _near interactions")
    paste(c("celltype",
            paste0("celltype:", spill_cols)),
          collapse = " + ")
  } else {
    paste(spill_cols, collapse = " + ")
  }
  ## Site / Treatment / Mutation confounders removed (user request 2026-06-01)
  fixed_rhs <- if (is.null(spill_part)) {
    paste(c("1", "Responder"), collapse = " + ")
  } else {
    paste(c("1", "Responder", spill_part), collapse = " + ")
  }
}
MUNDLAK <- nzchar(Sys.getenv("R_MUNDLAK"))
if (MUNDLAK) {
  for (tc in TYPES) df[[paste0("bar_", tc)]] <- ave(df[[tc]], df$imageID, df$celltype, FUN = mean)
  bar_cols <- paste0("bar_", TYPES)
  fixed_rhs <- paste(c(fixed_rhs, bar_cols, paste0("Responder:", bar_cols)), collapse = " + ")
}
X_fixed <- stats::model.matrix(stats::as.formula(paste("~", fixed_rhs)), data = df)
say("    X_fixed: %d x %d", nrow(X_fixed), ncol(X_fixed))

IMG_SLOPE <- nzchar(Sys.getenv("R_IMG_SLOPE"))
PAT_DZ    <- nzchar(Sys.getenv("R_PAT_SLOPE_DZ"))
img_terms <- {tt <- Sys.getenv("R_IMG_SLOPE_TERMS"); if (nzchar(tt)) strsplit(tt, ",")[[1]] else TYPES}
if (PAT_DZ || IMG_SLOPE) for (tc in img_terms) df[[paste0(tc, "_imgz")]] <- as.numeric(scale(df[[tc]]))
img_formula <- if (PAT_DZ) {
  pz <- paste0(img_terms, "_imgz")
  stats::as.formula(paste0("~ 1 + ", paste(c(pz, paste0("Responder:", pz)), collapse = " + ")))
} else if (IMG_SLOPE) {
  stats::as.formula(paste0("~ 1 + ", paste(paste0(img_terms, "_imgz"), collapse = " + ")))
} else {
  ~ 1
}
say("    imageID RE block: %s",
    if (PAT_DZ) "1 + baseline + Responder:baseline"
    else if (IMG_SLOPE) paste("1 +", paste(img_terms, collapse="+"))
    else "intercept only")
## R_TUM_ONLY_INT: constrain Responder × K_bio interactions to Tumour only,
## not all 7 neighbour celltypes. Tests whether per-focal V_RxS Fib→Tum
## (or Mac→Tum etc.) signal survives without the joint cross-neighbour
## disease modulation structure. Drops Responder:OtherCT random slopes.
TUM_ONLY_INT <- nzchar(Sys.getenv("R_TUM_ONLY_INT"))
## RE_TYPES is the neighbour-name set used in the random-effect formula (= TYPES).
re_formula_celltype <- if (TUM_ONLY_INT) {
  say("    CELLTYPE block: Tum-only disease interaction (R_TUM_ONLY_INT)")
  stats::as.formula(
    paste0("~ 1 + Responder + (", paste(RE_TYPES, collapse = " + "),
            ") + Responder:Tumour"))
} else {
  stats::as.formula(
    paste0("~ 1 + Responder * (", paste(RE_TYPES, collapse = " + "), ")"))
}
re_specs <- list(
  list(group_col = "celltype",
       formula   = re_formula_celltype),
  list(group_col = "imageID",
       formula   = img_formula))
Y <- as.matrix(df[, hvg_in_df, drop = FALSE]); storage.mode(Y) <- "integer"
offset_vec <- log(df$nCount)

## ---- 4b. (Optional) Expression-weighted technical spillover (Option B) ----
##
## E^tech[i, g] = Σ_{j != i, c(j) != c(i)} y_{j,g} × exp(-d_ij / h_tech)
##
## Cross-celltype-only diffusion kernel (h_tech=5 µm). Pre-fit a per-gene slope
## β_spill,g via NB / log-OLS, then pass pre_offset_mat[i,g] = β_spill,g × E^tech[i,g]
## to PACE as an additive log-offset. Replaces per-celltype _near fixed effects.
pre_offset_mat <- NULL
## Ambient-field carrier args for the per-cell contamination model
## (populated when R_BLEED_PERCELL=1; else NULL/0). Variable names kept for
## historical continuity; passed to the solver's ambient_* args.
ambient_E_tech_mat <- NULL
ambient_image_idx      <- NULL
ambient_n_images       <- 0L
## ---- Per-cell rho_i spillover: homotypic-core negative-control anchors ----
percell_anchor_mask <- NULL
if (nzchar(Sys.getenv("R_BLEED_PERCELL"))) {
  HC_HOMO <- as.numeric(Sys.getenv("R_PERCELL_HOMOFRAC", unset = "0.5"))
  ct_v <- as.character(df$celltype)
  type_means <- t(vapply(TYPES, function(tt) colMeans(Y[ct_v == tt, , drop = FALSE]), numeric(ncol(Y))))
  rownames(type_means) <- TYPES; colnames(type_means) <- colnames(Y)
  same_frac <- rep(NA_real_, nrow(df))
  for (im in levels(df$imageID)) {
    si <- which(df$imageID == im); if (length(si) < 50) next
    cl <- ct_v[si]; fr <- dbscan::frNN(as.matrix(df[si, c("x","y")]), eps = 30)
    for (k in seq_along(si)) { idk <- fr$id[[k]]; same_frac[si[k]] <- if (length(idk)) mean(cl[idk] == cl[k]) else 1 }
  }
  core <- which(same_frac >= HC_HOMO); core_means <- type_means; ncore <- integer(0)
  for (X in TYPES) { idx <- core[ct_v[core] == X]; ncore <- c(ncore, length(idx)); if (length(idx) >= 20) core_means[X, ] <- colMeans(Y[idx, , drop = FALSE]) }
  owner_mean <- apply(core_means, 2, max); owner_t <- TYPES[apply(core_means, 2, which.max)]
  percell_anchor_mask <- matrix(0, nrow(df), ncol(Y)); ngc <- integer(0)
  for (X in TYPES) {
    is_anc <- owner_t != X & owner_mean > 0.1 & (core_means[X, ] / pmax(owner_mean, 1e-9) < 0.1)
    rr <- which(ct_v == X); if (length(rr)) percell_anchor_mask[rr, ] <- matrix(as.numeric(is_anc), length(rr), ncol(Y), byrow = TRUE)
    ngc <- c(ngc, sum(is_anc))
  }
  say("    R_BLEED_PERCELL (Mel homotypic-core): core %d/%d; anchors/type median=%d range=[%d,%d]",
      length(core), nrow(df), as.integer(median(ngc)), min(ngc), max(ngc))
}
if (nzchar(Sys.getenv("R_E_TECH"))) {
  say("[4b] Building E^tech (cross-celltype expression-weighted spillover) ...")
  t_et <- Sys.time()
  ## frNN per imageID
  E_tech <- matrix(0, nrow(df), length(hvg_in_df),
                    dimnames = list(NULL, hvg_in_df))
  rad <- 3 * H_TECH  ## 15 µm
  df_idx_by_image <- split(seq_len(nrow(df)), df$imageID)
  ct_all <- as.character(df$celltype)
  for (im_idx in seq_along(df_idx_by_image)) {
    rows <- df_idx_by_image[[im_idx]]
    if (length(rows) < 2) next
    coords_im <- as.matrix(df[rows, c("x","y")])
    nn <- dbscan::frNN(coords_im, eps = rad)
    Y_im <- Y[rows, , drop = FALSE]
    ct_im <- ct_all[rows]
    ## For each focal, sum y_j × exp(-d/h_tech) over j != i with c(j) != c(i)
    for (i in seq_along(rows)) {
      nbr_local <- nn$id[[i]]
      if (!length(nbr_local)) next
      d_local   <- nn$dist[[i]]
      ## Cross-celltype only
      keep <- ct_im[nbr_local] != ct_im[i]
      if (!any(keep)) next
      w <- exp(-d_local[keep] / H_TECH)
      ## Sparse-ish accumulation: w-weighted sum of neighbour Y rows
      E_tech[rows[i], ] <- as.numeric(crossprod(w, Y_im[nbr_local[keep], , drop = FALSE]))
    }
    if (im_idx %% 5 == 0) cat(sprintf("    [E^tech] image %d/%d done\n",
                                       im_idx, length(df_idx_by_image)))
  }
  say("    E^tech computed in %.1f s (range = [%.2f, %.2f])",
      as.numeric(difftime(Sys.time(), t_et, units = "secs")),
      min(E_tech), max(E_tech))

  ## R_EDGE_CORRECT: apply isotropic area-fraction correction to E^tech per
  ## the manuscript spec.  E^tech is radius-based (frNN with eps = rad), so
  ## r = rad = 3 × H_TECH = 15 µm is unambiguously the search radius.
  ##     area_fraction(i) = (disc(i, 15 µm) ∩ image_rectangle) / (π × 15²)
  ##     E^tech_corrected = E^tech_raw / area_fraction
  if (nzchar(Sys.getenv("R_EDGE_CORRECT"))) {
    ## Recompute area_fraction at the E^tech radius (= rad = 3 × H_TECH)
    af_etech <- numeric(nrow(df))
    df_idx_by_im_corr <- split(seq_len(nrow(df)), df$imageID)
    for (im_idx in seq_along(df_idx_by_im_corr)) {
      rows <- df_idx_by_im_corr[[im_idx]]
      if (!length(rows)) next
      coords_im <- as.matrix(df[rows, c("x","y")])
      xmin <- min(coords_im[, 1]); xmax <- max(coords_im[, 1])
      ymin <- min(coords_im[, 2]); ymax <- max(coords_im[, 2])
      af_etech[rows] <- .compute_area_fraction(coords_im, rad, xmin, xmax, ymin, ymax)
    }
    say("    R_EDGE_CORRECT (E^tech): af median=%.3f  q10=%.3f  q90=%.3f  (%.1f%% cells with af<0.95)",
        median(af_etech), quantile(af_etech, 0.10), quantile(af_etech, 0.90),
        100 * mean(af_etech < 0.95))
    E_tech <- sweep(E_tech, 1, af_etech, "/")
    say("    Edge-corrected E^tech range = [%.2f, %.2f]", min(E_tech), max(E_tech))
  }

  ## log_E_tech is the per-(cell, gene) covariate used by BOTH paths.
  ##   - global prefit path     : β_spill,g × log_E_tech  → pre_offset_mat
  ##   - per-image RE path      : log_E_tech IS the covariate; δ_{m, g} is fit
  ##                              jointly inside PACE
  log_E_tech <- log1p(E_tech)
  storage.mode(log_E_tech) <- "double"

  if (nzchar(Sys.getenv("R_BLEED_PERCELL"))) {
    ## Per-cell rho_i contamination loading on local ambient a_ig = RAW E^tech.
    say("    R_BLEED_PERCELL=1: per-cell shrunken contamination loading on local ambient (RAW E^tech).")
    ambient_E_tech_mat <- E_tech
    storage.mode(ambient_E_tech_mat) <- "double"
    ambient_image_idx      <- as.integer(df$imageID)
    ambient_n_images       <- length(levels(df$imageID))
    pre_offset_mat <- NULL
  } else {
    ## ============================================================
    ## Existing global gene-wise NB1 prefit path (unchanged).
    ## ============================================================
    ## R_ETECH_PREFIT: pre-fit method for β_spill,g (default = NB1 via glmmTMB,
    ## matching PACE's NB1 parameterization). "ols" = legacy log-OLS surrogate.
    ## R_ETECH_PREFIT_WORKERS: parallel workers for the per-gene loop (default 6).
    PREFIT_METHOD  <- Sys.getenv("R_ETECH_PREFIT", unset = "nb1")
    PREFIT_WORKERS <- as.integer(Sys.getenv("R_ETECH_PREFIT_WORKERS", unset = "6"))
    if (PREFIT_METHOD == "nb1" && !requireNamespace("glmmTMB", quietly = TRUE))
      stop("R_ETECH_PREFIT=nb1 requires the glmmTMB package")
    say("    Pre-fitting β_spill,g via gene-wise %s (%d worker%s) ...",
        toupper(PREFIT_METHOD), PREFIT_WORKERS, ifelse(PREFIT_WORKERS > 1, "s", ""))
    t_prefit <- Sys.time()
    PREFIT_BPPARAM <- if (PREFIT_WORKERS > 1) {
      BiocParallel::MulticoreParam(workers = PREFIT_WORKERS, RNGseed = 1L)
    } else {
      BiocParallel::SerialParam()
    }
    ## Fix B: log-linear bleed correction.
    ## Pre-fit:  y ~ log(1 + E_tech) + offset(log N)
    ## Pre-offset:  β_spill,g × log(1 + E_tech[i, g])
    ## The log(1+E) form bounds the correction: the multiplicative correction
    ## (1+E)^β scales polynomially rather than exponentially, so it saturates
    ## naturally for high-bleed regions.  No ±3 cap required.
    one_gene <- function(gi) {
      g   <- hvg_in_df[gi]
      e_g <- log_E_tech[, g]
      if (stats::sd(e_g) < 1e-10) return(list(beta = 0, ok = FALSE, fail = FALSE))
      if (PREFIT_METHOD == "nb1") {
        df_gi <- data.frame(y = Y[, g], e = e_g, off = offset_vec)
        fit_nb <- tryCatch(
          suppressWarnings(glmmTMB::glmmTMB(
            y ~ e + offset(off),
            family = glmmTMB::nbinom1(),
            data   = df_gi)),
          error = function(err) NULL)
        converged <- !is.null(fit_nb) &&
                      isTRUE(fit_nb$sdr$pdHess) &&
                      is.finite(glmmTMB::fixef(fit_nb)$cond["e"])
        if (converged)
          return(list(beta = as.numeric(glmmTMB::fixef(fit_nb)$cond["e"]),
                      ok = TRUE, fail = FALSE))
        list(beta = stats::lm.fit(cbind(1, e_g),
                                   log(Y[, g] + 0.5) - offset_vec)$coefficients[2],
             ok = FALSE, fail = TRUE)
      } else {
        list(beta = stats::lm.fit(cbind(1, e_g),
                                   log(Y[, g] + 0.5) - offset_vec)$coefficients[2],
             ok = FALSE, fail = FALSE)
      }
    }
    res_list <- BiocParallel::bplapply(seq_along(hvg_in_df), one_gene,
                                         BPPARAM = PREFIT_BPPARAM)
    beta_spill <- vapply(res_list, `[[`, numeric(1), "beta")
    names(beta_spill) <- hvg_in_df
    n_nb_ok   <- sum(vapply(res_list, `[[`, logical(1), "ok"))
    n_nb_fail <- sum(vapply(res_list, `[[`, logical(1), "fail"))
    if (PREFIT_METHOD == "nb1")
      say("    NB1(glmmTMB) log-linear converged for %d / %d genes (%.1f%%); %d fell back to OLS  (pre-fit took %.1f min)",
          n_nb_ok, length(hvg_in_df), 100 * n_nb_ok / length(hvg_in_df), n_nb_fail,
          as.numeric(difftime(Sys.time(), t_prefit, units = "mins")))
    say("    β_spill,g range = [%.3f, %.3f]  median = %.3f",
        min(beta_spill), max(beta_spill), median(beta_spill))

    pre_offset_mat <- sweep(log_E_tech, 2, beta_spill, "*")
    say("    pre_offset_mat range = [%.3f, %.3f]  (log-linear; NO cap applied)",
        min(pre_offset_mat), max(pre_offset_mat))

    ## ============================================================
    ## Option B (Phase 4 fix, 2026-06-09): per-image scalar λ_m.
    ## ============================================================
    ## Activated by R_BLEED_PER_IMAGE_LAMBDA=1.  Adds a per-image multi-
    ## plicative scaling to the GLOBAL pre_offset_mat:
    ##     pre_offset_mat[i, g]  ←  λ_{m(i)} × β_spill,g × log(1+E^tech_ig)
    ## λ_m is estimated by per-image FWL: regress observed log expression on
    ## the global bleed_offset within each image, after partialling out
    ## per-image celltype-mean structure (so λ_m captures variation NOT
    ## explainable by per-image celltype composition).  18 parameters total.
    ##
    ## Compatible with R_BLEED_PER_IMAGE_RE=0 (i.e. NOT the per-(image, gene)
    ## RE block path); the two env vars are mutually exclusive.
    if (nzchar(Sys.getenv("R_BLEED_PER_IMAGE_LAMBDA"))) {
      if (nzchar(Sys.getenv("R_BLEED_PER_IMAGE_RE")))
        stop("R_BLEED_PER_IMAGE_LAMBDA and R_BLEED_PER_IMAGE_RE are mutually exclusive")
      say("    R_BLEED_PER_IMAGE_LAMBDA=1: fitting per-image scalar λ_m by FWL.")
      images_local <- levels(df$imageID)
      lambda_m <- numeric(length(images_local))
      names(lambda_m) <- images_local
      n_genes_used_m <- integer(length(images_local))
      for (mi in seq_along(images_local)) {
        cells_m <- which(as.integer(df$imageID) == mi)
        if (length(cells_m) < 50) { lambda_m[mi] <- 1; next }
        ct_m <- factor(as.character(df$celltype[cells_m]))
        X_ctrl <- if (nlevels(ct_m) >= 2)
          stats::model.matrix(~ 1 + ct_m)
        else
          matrix(1, length(cells_m), 1)
        QR <- qr(X_ctrl)
        residfn <- function(v) v - X_ctrl %*% qr.solve(QR, v)
        cum_cov <- 0; cum_var <- 0; n_used <- 0L
        for (gi in seq_along(hvg_in_df)) {
          g <- hvg_in_df[gi]
          x_g <- pre_offset_mat[cells_m, gi]
          if (stats::sd(x_g) < 1e-10) next
          y_g <- log(Y[cells_m, gi] + 0.5) - offset_vec[cells_m]
          y_r <- residfn(y_g); x_r <- residfn(x_g)
          cum_cov <- cum_cov + sum(y_r * x_r)
          cum_var <- cum_var + sum(x_r * x_r)
          n_used  <- n_used  + 1L
        }
        lambda_m[mi]   <- if (cum_var > 1e-10) cum_cov / cum_var else 1
        n_genes_used_m[mi] <- n_used
      }
      say("    λ_m range = [%.3f, %.3f]   median %.3f",
          min(lambda_m), max(lambda_m), median(lambda_m))
      say("    λ_m per image:")
      for (mi in seq_along(images_local))
        cat(sprintf("      %-15s  λ = %+.4f   (n_genes_used = %d)\n",
                    images_local[mi], lambda_m[mi], n_genes_used_m[mi]))
      ## Apply the scaling
      cell_image_idx <- as.integer(df$imageID)
      pre_offset_mat <- pre_offset_mat * lambda_m[cell_image_idx]
      say("    Post-λ pre_offset_mat range = [%.3f, %.3f]",
          min(pre_offset_mat), max(pre_offset_mat))
    }
    rm(log_E_tech); invisible(gc(verbose = FALSE))
  }
}

## ---- 5. Joint multi-block fit ----
## R_PQL_THREADS: C++ OpenMP threads in the per-chunk solver (default 4).
## These are kernel threads sharing one R process => ZERO memory overhead per
## additional thread.
## NOTE: do NOT combine with MulticoreParam (BPPARAM stays SerialParam).
## Combining process forking with OpenMP would over-subscribe CPU AND multiply
## fork memory; pick one parallelism layer. n_threads is the right one here
## because the chunked C++ solver IS the bottleneck.
## Tested 8 threads on a 10-core Mac: zero additional speedup (PQL is
## memory-bandwidth bound, not CPU-bound at this dataset size). 4 is the
## empirical sweet spot.
NO_BLEED <- nzchar(Sys.getenv("R_NO_BLEED")) || C3
PQL_THREADS <- as.integer(Sys.getenv("R_PQL_THREADS", unset = "4"))
say("[5/7] fit_pace_mvpql_joint_multi (h_tech=%.1f, h_bio=%.1f, %d C++ threads) %s%s...",
    H_TECH, H_BIO, PQL_THREADS,
    if (NO_BLEED) "(NO_BLEED)" else "",
    if (!is.null(pre_offset_mat)) " (E^tech)" else "")

## R_DATA_INFORMED_TAU: build raw-data-informed weights for the random-slope
## prior variance.  See .compute_data_informed_weights() in helpers/pace_mvpql.R.
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
  bleed_percell = nzchar(Sys.getenv("R_BLEED_PERCELL")),
  percell_anchor_mask = percell_anchor_mask,
  offset_vec = offset_vec,
  pre_offset_mat = pre_offset_mat,
  data_informed_W = data_informed_W,
  ambient_mat        = ambient_E_tech_mat,
  ambient_image_idx  = ambient_image_idx,
  ambient_n_images   = ambient_n_images,
  n_iter = as.integer(Sys.getenv("R_N_ITER", unset = "16")), tol = 5e-3,
  tau_shrinkage = Sys.getenv("R_TAU", unset = "adaptive"),
  BPPARAM = SerialParam(), n_threads = PQL_THREADS,
  interior_precision = 1L, chunk_size = 128L, verbose = TRUE)
say("    fit done in %.1f min", as.numeric(difftime(Sys.time(), t_fit, units = "mins")))

if (nzchar(Sys.getenv("R_RAW_FIT_ONLY"))) {
  dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fit = fit, gene_set = colnames(Y)), OUT)
  say("Saved RAW FIT ONLY to %s", OUT); quit(save = "no")
}

## ---- 6. mashr + decomposition ----
say("[6/7] mashr + decomposition ...")
results_mv <- mvpql_to_results_multi(fit, keep_block = "celltype")
shrunken_long <- apply_mashr_shrinkage(
  results = results_mv, focals = TYPES, neighbours = TYPES,
  resp_term = RESP_TERM)
dec <- tryCatch(
  mvpql_variance_decomposition_multi(
    fit = fit, df = df, Y = Y,
    vars = TYPES, X_fixed = X_fixed,
    resp_term = RESP_TERM, focal_levels = TYPES),
  error = function(e) {
    say("WARNING: decomposition failed (%s); saving NULL", conditionMessage(e))
    NULL
  })

## ---- 7. MCSD per pair (collapsed-7 version) ----
say("[7/7] MCSD + R^2_RxS ...")
PAIRS_MEL <- list(c("Macrophage","Tumour"), c("Fibroblast","Tumour"),
                  c("Endothelial","Tumour"), c("T_CD8_memory","Tumour"),
                  c("Treg","Tumour"), c("B_Cell","Tumour"))
PAIRS_MEL <- Filter(function(p) all(p %in% TYPES), PAIRS_MEL)
gene_focal_5block <- if (!is.null(dec)) dec$gene_focal_5block else NULL
if (is.null(gene_focal_5block)) {
  say("    no decomposition -> skipping MCSD (fit + shrunken_long still saved)")
  mcsd_canonical <- list()
}
if (!is.null(gene_focal_5block)) {

Z_re <- fit$re_meta$Z
blk_ct <- fit$re_meta$blocks[[which(vapply(fit$re_meta$blocks, `[[`,
                                            character(1), "group_col") == "celltype")]]
cells_by_ct <- setNames(fit$re_meta$cells_by_grp_list[[
  which(vapply(fit$re_meta$blocks, `[[`, character(1), "group_col") == "celltype")]],
  blk_ct$group_levels)
alpha_g <- pmax(fit$alpha, 0)
gene_names_fit <- colnames(fit$U)
names(alpha_g) <- gene_names_fit

mcsd_canonical <- list()
for (p in PAIRS_MEL) {
  fc <- p[1]; nc <- p[2]; term <- paste0(RESP_TERM, ":", nc); pk <- paste(fc, nc, sep="_")
  s <- shrunken_long |> dplyr::filter(focal == fc, neighbour == nc, term == !!term) |>
    dplyr::distinct(gene, .keep_all = TRUE)
  fm <- gene_focal_5block |> dplyr::filter(focal == fc) |>
    dplyr::select(gene, spec, focal_mean) |> dplyr::distinct(gene, .keep_all = TRUE)
  if (!nrow(s) || !nrow(fm)) { cat(sprintf("  %s: skip\n", pk)); next }

  cells_c <- cells_by_ct[[fc]]
  col_N <- match(paste0(fc, "::", nc), colnames(Z_re))
  N_t <- as.numeric(Z_re[cells_c, col_N])
  var_N <- stats::var(N_t, na.rm = TRUE)
  col_R <- match(paste0(fc, "::", RESP_TERM), colnames(Z_re))
  R_c <- if (is.na(col_R)) df$.resp_dummy[cells_c] else as.numeric(Z_re[cells_c, col_R])
  var_RN <- stats::var(R_c * N_t, na.rm = TRUE)
  mu_bar_per_gene <- Matrix::colMeans(fit$mu[cells_c, , drop = FALSE])
  names(mu_bar_per_gene) <- gene_names_fit
  row_u <- match(paste0(fc, "::", nc), rownames(fit$U))
  u_vec <- if (is.na(row_u)) setNames(rep(0, length(gene_names_fit)), gene_names_fit)
           else setNames(as.numeric(fit$U[row_u, ]), gene_names_fit)

  res_all <- s |> dplyr::inner_join(fm, by = "gene") |>
    dplyr::mutate(
      MCSD     = (estimate_shrunk^2) * (spec^2) * pmax(focal_mean, 0),
      MCSD4    = (estimate_shrunk^2) * (spec^4) * pmax(focal_mean, 0),
      mu_bar   = mu_bar_per_gene[gene],
      alpha    = alpha_g[gene],
      u_raw    = u_vec[gene],
      V_S      = u_raw^2 * var_N,
      V_RxS    = estimate_shrunk^2 * var_RN,
      V_resid  = log(1 + (1 + pmax(alpha, 0)) / pmax(mu_bar, 1e-6)),
      V_total  = V_S + V_RxS + V_resid,
      R2_S     = V_S    / pmax(V_total, 1e-12),
      R2_RxS   = V_RxS  / pmax(V_total, 1e-12)
    )
  res <- res_all |>
    dplyr::filter(lfsr < 0.05) |>
    dplyr::arrange(dplyr::desc(MCSD)) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::rename(b_clean = estimate_shrunk) |>
    dplyr::select(rank, gene, MCSD, MCSD4, b_clean, u_raw, spec, focal_mean,
                  R2_S, R2_RxS, mu_bar, alpha,
                  V_S, V_RxS, V_resid, V_total, lfsr, sd_shrunk)
  status <- if (nrow(res) >= 3) "significant" else "honestly null"
  mcsd_canonical[[pk]] <- list(scores = res, status = status)
  cat(sprintf("  %s n_sig=%d top-10 (%s): %s\n",
              pk, nrow(res), status, paste(head(res$gene, 10), collapse=", ")))
}
} ## end if (!is.null(gene_focal_5block))

cell_meta <- data.frame(cell_id = as.character(seq_len(nrow(df))),
                        nCount = df$nCount,
                        Condition = as.character(df$Responder),
                        patientID = as.character(df$imageID),
                        imageID = as.character(df$imageID),
                        celltype = as.character(df$celltype),
                        stringsAsFactors = FALSE)
say("[+] intrinsic patient-level disease test (pace_twostage_z) ...")
source("scripts/helpers/pace_twostage_z.R")
twostage_disease <- tryCatch(
  pace_twostage_z(fit, raw_Y = Y,
    cluster_vec = as.character(df$imageID), cond_vec = as.character(df$Responder),
    case_level = RESP_CASE, libsize = df$nCount),
  error = function(e) { say("twostage failed: %s", conditionMessage(e)); NULL })
if (!is.null(twostage_disease))
  say("    %d (focal,nbr,gene) triples ; q_pair<0.10 = %d",
      nrow(twostage_disease), sum(twostage_disease$q_pair < 0.10, na.rm = TRUE))

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = fit, shrunken_long = shrunken_long, gene_set = colnames(Y),
              decomposition = dec, mcsd_canonical = mcsd_canonical,
              twostage_disease = twostage_disease,
              h_tech = H_TECH, h_bio = H_BIO,
              K_tech = K_tech, K_bio = K_bio, cell_meta = cell_meta,
              resolution = "kernel_refactor_collapsed7"),
        OUT)
say("Saved %s", OUT)
say("TOTAL: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

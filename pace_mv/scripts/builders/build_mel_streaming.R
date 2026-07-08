## build_mel_streaming.R -- MELANOMA canonical (collapsed-7CT, Responder PD-vs-nonPD)
## refit via the MEMORY-BOUNDED streaming solver (fit_pace_mvpql_streaming).
##
## SURGICAL PRINCIPLE
## ------------------
## This is a verbatim port of the DENSE Mel canonical builder
##   scripts/builders/mel_fit_joint_kernel_pathA_collapsed7_redesign.R
## with EXACTLY ONE change to the model machinery: the dense joint solver
##   fit_pace_mvpql_joint_multi(... ambient_mat = E_tech ...)
## is replaced by the streaming solver
##   fit_pace_mvpql_streaming(... ambient_W = W ...)
## where W is the SPARSE n x n E^tech weight matrix such that
##   ambient_W %*% Y  ==  dense E_tech   (exact identity, spot-checked below).
##
## EVERYTHING ELSE is identical to the dense Mel builder: h5ad load, pd_nonpd
## response coding, Tumor_a..f -> Tumour collapse, the 7 cell types, the
## per-image frNN kernel (K_tech / K_bio), the X_fixed design (Responder main
## effect + Site/Treatment/Mutation confounders + interactions under the
## canonical R_E_TECH path), the re_specs (celltype block WITH the Responder
## interaction + imageID block), the per-image E^tech ambient, homotypic-core
## per-cell anchors, the mashr/decomposition/MCSD downstream WITH resp_term, and
## the saved list schema.
##
## CANONICAL ENV (run with):
##   R_E_TECH=1 R_BLEED_PERCELL=1 R_DISP_MODEL=nb2 \
##   R_PAT_SLOPE_DZ=1 R_TAU=adaptive R_DATA_INFORMED_TAU=1 \
##   R_N_ITER=32 R_PQL_THREADS=8 R_FUSE_RHO=1 R_CHUNK_SIZE=128 \
##   R_MEL_RESP=pd_nonpd \
##   R_OUT=data/simvi_melanoma/sweeps/mvpql_mel_streaming.rds \
##   Rscript scripts/builders/build_mel_streaming.R
##
## Watch progress:  tail -f <your log>
##
## OUTPUT: data/simvi_melanoma/sweeps/mvpql_mel_streaming.rds
##   (NON-canonical path; does NOT overwrite mvpql_percell_hc.rds).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(stringr); library(forcats); library(tidyr)
  library(SingleCellExperiment); library(Matrix); library(BiocParallel)
  library(zellkonverter); library(scran); library(FNN); library(dbscan)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
source("scripts/zzz.R"); .load_all("scripts")
source("scripts/helpers/pace_mvpql.R")
source("scripts/helpers/pace_mvpql_joint_multi.R")
source("streaming/helpers/pace_mvpql_streaming.R")

## ---- Small logging helpers (timestamped so a `tail -f` watch is useful) ----
ts  <- function() format(Sys.time(), "%H:%M:%S")
say <- function(...) cat(sprintf("[%s] %s\n", ts(), sprintf(...)))
gc_peak <- function(tag) {
  g <- gc(verbose = FALSE)
  say("    [mem] %-22s peak Ncells+Vcells max used = %.0f MB",
      tag, sum(g[, "max used"] * c(56, 8)) / 1e6)
}

t0 <- Sys.time()

## ---- Kernel bandwidths (canonical defaults; ported verbatim) ----
H_TECH <- as.numeric(Sys.getenv("R_H_TECH", unset = "5"))
H_BIO  <- as.numeric(Sys.getenv("R_H_BIO",  unset = "30"))

## ---- Isotropic area-fraction edge correction (copied verbatim from the dense
## Mel builder; applied to the E^tech kernel at r = 3 * h_tech). ----
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
    af[i] <- (pi * sum(r_eff^2) / N_angles) / norm_factor
  }
  af
}

OUT <- Sys.getenv("R_OUT",
  unset = "data/simvi_melanoma/sweeps/mvpql_mel_streaming.rds")
say("Mel STREAMING (collapsed-7CT, Responder): h_tech=%.1f h_bio=%.1f  out=%s",
    H_TECH, H_BIO, OUT)

## COLLAPSED 7-CT scheme: pool Tumor_a..f into a single "Tumour" focal.
TYPES <- c("Tumour",
           "Endothelial", "Fibroblast", "Macrophage",
           "B_Cell", "T_CD8_memory", "Treg")
TUM_RAW <- c("Tumor_a","Tumor_b","Tumor_c","Tumor_d","Tumor_e","Tumor_f")

## ============================================================================
## 1. Load + clean SCE (ported verbatim from the dense Mel builder).
## ============================================================================
say("[1/7] Loading Mel SCE ...")
sce <- readH5AD(file = "data/simvi_melanoma/Melanoma_5612.h5ad", reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
px_size <- 0.12028
set.seed(1994)

## Response coding (R_MEL_RESP):
##   pd_nonpd (CANONICAL): keep ALL (PD/SD/CR/PR); PD vs non-PD(SD+CR+PR).
##   pd_resp: drop SD; PD vs RESP(CR+PR). (sensitivity)
##   binary_sd: keep SD; NonResp(PD+SD) vs Resp(PR+CR). (sensitivity)
MEL_RESP <- Sys.getenv("R_MEL_RESP", unset = "pd_nonpd")
if (identical(MEL_RESP, "binary_sd")) {
  say("    R_MEL_RESP=binary_sd: keeping SD; NonResp=PD+SD vs Resp=PR+CR")
  sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR"))]
  sce$BEST_RESPONSE_BY_SCAN <- ifelse(
    as.character(sce$BEST_RESPONSE_BY_SCAN) %in% c("CR","PR"), "Resp", "NonResp")
  RESP_REF <- "Resp"; RESP_CASE <- "NonResp"
} else if (identical(MEL_RESP, "pd_resp")) {
  say("    R_MEL_RESP=pd_resp: dropping SD; PD vs RESP(CR+PR)")
  sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "CR", "PR"))]
  sce$BEST_RESPONSE_BY_SCAN <- ifelse(
    as.character(sce$BEST_RESPONSE_BY_SCAN) %in% c("CR","PR"),
    "RESP", as.character(sce$BEST_RESPONSE_BY_SCAN))
  RESP_REF <- "RESP"; RESP_CASE <- "PD"
} else {
  say("    R_MEL_RESP=pd_nonpd (CANONICAL): keeping ALL; PD vs non-PD(SD+CR+PR)")
  sce <- sce[, which(sce$BEST_RESPONSE_BY_SCAN %in% c("PD", "SD", "CR", "PR"))]
  sce$BEST_RESPONSE_BY_SCAN <- ifelse(
    as.character(sce$BEST_RESPONSE_BY_SCAN) == "PD", "PD", "nonPD")
  RESP_REF <- "nonPD"; RESP_CASE <- "PD"
}
RESP_TERM <- paste0("Responder", RESP_CASE)   # ResponderPD (pd_nonpd / pd_resp)

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

## ============================================================================
## 2. Per-image kernel-weighted neighbour predictors (frNN, radius, uncapped).
##    Ported verbatim: PER-IMAGE frNN (no cross-FOV mixing), eps = 3 * h_bio.
## ============================================================================
say("[2/7] Building kernel-weighted predictors per image (frNN eps=%d um, uncapped) ...",
    as.integer(3*H_BIO))
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
## sparseMatrix SUMS duplicate (cell, celltype) entries -> correct accumulation.
K_tech <- as.matrix(Matrix::sparseMatrix(i = I_acc, j = J_acc, x = WT_acc, dims = c(n, length(TYPES))))
K_bio  <- as.matrix(Matrix::sparseMatrix(i = I_acc, j = J_acc, x = WB_acc, dims = c(n, length(TYPES))))
colnames(K_tech) <- paste0(TYPES, "_near"); colnames(K_bio) <- TYPES
rm(I_acc, J_acc, WT_acc, WB_acc)
say("    K_tech mean range=[%.3f, %.3f]; K_bio mean range=[%.3f, %.3f]",
    min(colMeans(K_tech)), max(colMeans(K_tech)),
    min(colMeans(K_bio)),  max(colMeans(K_bio)))

## ============================================================================
## 3. Detection filter (det>=5% in any focal) + df assembly (ported verbatim).
## ============================================================================
say("[3/7] det>=5%% gene filter + assembling df ...")
hvg_in_df <- intersect(hvg_genes, names(df_raw))
ct_chr <- as.character(df_raw$celltype)
ct_idx_lst <- lapply(TYPES, function(c) which(ct_chr == c))
chunk_size_det <- 50L
max_det <- numeric(length(hvg_in_df)); names(max_det) <- hvg_in_df
for (s in seq.int(1L, length(hvg_in_df), by = chunk_size_det)) {
  e <- min(s + chunk_size_det - 1L, length(hvg_in_df))
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
n_pat <- df |> dplyr::distinct(imageID, Responder) |> dplyr::count(Responder)
say("    Per-imageID Responder count: %s",
    paste(n_pat$Responder, n_pat$n, sep="=", collapse=", "))

## ---- Patient-level confounders (ported verbatim) ----
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

## K_bio neighbour columns onto df (matched by cell_id, as the dense builder).
for (tc in TYPES) df[[tc]] <- K_bio[match(df$cell_id, df_raw$cell_id), tc]

## ---- R_DROP_SPARSE_K (ported verbatim; optional) ----
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
  say("    dropped %d (focal, nb) pairs", length(dropped))
}

## ---- R_K_WITHIN (ported verbatim; optional) ----
K_WITHIN <- nzchar(Sys.getenv("R_K_WITHIN"))
if (K_WITHIN) {
  say("    R_K_WITHIN: replacing K_<nb> with within-(imageID, celltype) centered")
  for (tc in TYPES) {
    df[[tc]] <- df[[tc]] - ave(df[[tc]], df$imageID, df$celltype, FUN = mean)
  }
  max_abs_mean <- max(vapply(TYPES, function(tc) {
    max(abs(ave(df[[tc]], df$imageID, df$celltype, FUN = mean)))
  }, numeric(1)))
  say("    K_<nb>_within max |within-group mean| = %.2e", max_abs_mean)
}

## NOTE: R_K_PCA / R_DECORR_K_BIO_TARGET / R_MUNDLAK paths from the dense builder
## are NOT canonical (R_E_TECH default leaves them off) and are intentionally
## omitted here -- they are sensitivity-only knobs and play no role in the
## locked Mel streaming refit.

for (tc in TYPES) df[[paste0(tc, "_imgz")]] <- as.numeric(scale(df[[tc]]))
say("    df: %s cells, %d genes, %d images", format(nrow(df), big.mark=","),
    length(hvg_in_df), nlevels(df$imageID))

## ---- DIMENSION ASSERTION: must match the dense canonical (56261 x 928) ----
say("    [assert] df dims = %d cells x %d genes (dense canonical: 56261 x 928)",
    nrow(df), length(hvg_in_df))
if (nrow(df) != 56261L || length(hvg_in_df) != 928L) {
  say("    WARNING: dims (%d x %d) DIFFER from the dense canonical (56261 x 928). ",
      nrow(df), length(hvg_in_df))
  say("    WARNING: check h5ad / filters before trusting the streaming refit.")
}

## ============================================================================
## 4. X_fixed + re_specs (ported verbatim under the canonical R_E_TECH path).
## ============================================================================
say("[4/7] Building X_fixed + re_specs (canonical E^tech path) ...")
conf_intxn_rhs <- paste(c(
  paste0("Site:",      TYPES),
  paste0("Treatment:", TYPES),
  paste0("Mutation:",  TYPES)
), collapse = " + ")
conf_main_rhs <- "Site + Treatment + Mutation"

## Under the canonical R_E_TECH=1 path the per-celltype `_near` fixed effects are
## DROPPED (the per-cell rho_i contamination model handles spillover instead).
## We keep the same env gating as the dense builder so the design matches.
E_TECH <- nzchar(Sys.getenv("R_E_TECH"))
PER_FOCAL_NEAR <- nzchar(Sys.getenv("R_PER_FOCAL_NEAR"))
spill_cols <- paste0(TYPES, "_near")
for (sc in spill_cols) df[[sc]] <- K_tech[match(df$cell_id, df_raw$cell_id), sc]
spill_part <- if (E_TECH) {
  say("    EXPRESSION-WEIGHTED TECH SPILLOVER (R_E_TECH): dropping _near columns")
  NULL  ## drop _near entirely
} else if (PER_FOCAL_NEAR) {
  say("    PER-FOCAL beta_near: using celltype:* _near interactions")
  paste(c("celltype", paste0("celltype:", spill_cols)), collapse = " + ")
} else {
  paste(spill_cols, collapse = " + ")
}
fixed_rhs <- if (is.null(spill_part)) {
  paste(c("1", "Responder"), collapse = " + ")
} else {
  paste(c("1", "Responder", spill_part), collapse = " + ")
}
X_fixed <- stats::model.matrix(stats::as.formula(paste("~", fixed_rhs)), data = df)
say("    X_fixed: %d x %d", nrow(X_fixed), ncol(X_fixed))

## ---- imageID RE block (R_PAT_SLOPE_DZ canonical: 1 + baseline + Responder:baseline) ----
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

## ---- celltype RE block (Responder interaction; R_TUM_ONLY_INT optional) ----
TUM_ONLY_INT <- nzchar(Sys.getenv("R_TUM_ONLY_INT"))
RE_TYPES <- TYPES   ## neighbour names go into the RE formula (no PCA in canonical)
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
  list(group_col = "celltype", formula = re_formula_celltype),
  list(group_col = "imageID",  formula = img_formula))

## Counts matrix + offset (Y kept dense here for kernel/anchor precompute; a
## sparse copy is built later for the streaming solver).
Y <- as.matrix(df[, hvg_in_df, drop = FALSE]); storage.mode(Y) <- "integer"
offset_vec <- log(df$nCount)

## ============================================================================
## 4b. Per-cell contamination: homotypic-core anchors + SPARSE E^tech weight W.
##
## The dense Mel builder forms a DENSE n x G E^tech ambient (per-image, cross-
## celltype, weight exp(-d/h_tech), divided by the area-fraction edge correction
## af_i at r = 3*h_tech) and passes it as ambient_mat. The streaming solver wants
## the SPARSE n x n weight matrix W with W[i,j] = exp(-d_ij/h_tech)/af_i for every
## cross-type neighbour j of i within rad in the SAME image, so that
## ambient_W %*% Y == dense E_tech (identity spot-checked below).
## ============================================================================
ambient_W           <- NULL
ambient_image_idx   <- NULL
ambient_n_images    <- 0L
percell_anchor_mask <- NULL
percell_anchor_idx  <- NULL

PERCELL <- nzchar(Sys.getenv("R_BLEED_PERCELL"))
if (PERCELL) {
  ## ---- Homotypic-core anchor genes (ported verbatim from the dense Mel builder,
  ##      but stored in the MEMORY-LIGHT n_types x G + length-n index form that
  ##      the streaming solver accepts, mirroring build_bc_streaming.R) ----
  HC_HOMO <- as.numeric(Sys.getenv("R_PERCELL_HOMOFRAC", unset = "0.5"))
  ct_v <- as.character(df$celltype)
  type_means <- t(vapply(TYPES, function(tt) colMeans(Y[ct_v == tt, , drop = FALSE]),
                         numeric(ncol(Y))))
  rownames(type_means) <- TYPES; colnames(type_means) <- colnames(Y)

  ## Same-type fraction within a 30 um frNN ball, PER IMAGE (Mel imageID = SPID_fov).
  same_frac <- rep(NA_real_, nrow(df))
  for (im in levels(df$imageID)) {
    si <- which(df$imageID == im); if (length(si) < 50) next
    cl <- ct_v[si]
    fr <- dbscan::frNN(as.matrix(df[si, c("x","y")]), eps = 30)
    for (k in seq_along(si)) {
      idk <- fr$id[[k]]
      same_frac[si[k]] <- if (length(idk)) mean(cl[idk] == cl[k]) else 1
    }
  }
  core <- which(same_frac >= HC_HOMO)
  core_means <- type_means
  ncore <- integer(0)
  for (X in TYPES) {
    idx <- core[ct_v[core] == X]; ncore <- c(ncore, length(idx))
    if (length(idx) >= 20) core_means[X, ] <- colMeans(Y[idx, , drop = FALSE])
  }
  say("    R_BLEED_PERCELL (Mel homotypic-core): core %d/%d cells (>=%.0f%% same-type)",
      length(core), nrow(df), 100*HC_HOMO)
  owner_mean <- apply(core_means, 2, max)
  owner_t    <- TYPES[apply(core_means, 2, which.max)]
  ## n_types x G anchor mask + length-n celltype index (memory-light form).
  percell_anchor_mask <- matrix(0, length(TYPES), ncol(Y),
                                dimnames = list(TYPES, colnames(Y)))
  ngc <- integer(0)
  for (ti in seq_along(TYPES)) {
    X <- TYPES[ti]
    is_anc <- owner_t != X & owner_mean > 0.1 & (core_means[X, ] / pmax(owner_mean, 1e-9) < 0.1)
    percell_anchor_mask[ti, ] <- as.numeric(is_anc)
    ngc <- c(ngc, sum(is_anc))
  }
  percell_anchor_idx <- as.integer(df$celltype)
  say("    anchors per type median=%d range=[%d,%d]", as.integer(median(ngc)), min(ngc), max(ngc))
}

if (E_TECH) {
  say("[4b] Building SPARSE E^tech weight matrix W (per image) ...")
  t_et <- Sys.time()
  rad <- 3 * H_TECH   ## 15 um
  by_image <- split(seq_len(nrow(df)), df$imageID)
  ct_all   <- as.character(df$celltype)

  ## ---- edge area-fraction at r = rad, PER IMAGE, if requested ----
  EDGE <- nzchar(Sys.getenv("R_EDGE_CORRECT"))
  af_etech <- rep(1, nrow(df))
  if (EDGE) {
    for (ii_im in seq_along(by_image)) {
      rows <- by_image[[ii_im]]
      if (!length(rows)) next
      coords_im <- as.matrix(df[rows, c("x","y")])
      af_etech[rows] <- .compute_area_fraction(
        coords_im, rad,
        min(coords_im[,1]), max(coords_im[,1]),
        min(coords_im[,2]), max(coords_im[,2]))
    }
    say("    R_EDGE_CORRECT: af median=%.3f q10=%.3f q90=%.3f (%.1f%% cells af<0.95)",
        median(af_etech), quantile(af_etech, 0.10), quantile(af_etech, 0.90),
        100 * mean(af_etech < 0.95))
  }

  ## ---- build (i, j, w) triplets per image:
  ##      W[i,j] = exp(-d_ij/h_tech) / af[i] for cross-celltype neighbours j
  ##      within rad of i (self excluded; frNN excludes self) ----
  ii_all <- vector("list", length(by_image))
  jj_all <- vector("list", length(by_image))
  ww_all <- vector("list", length(by_image))
  for (ii_im in seq_along(by_image)) {
    rows <- by_image[[ii_im]]
    if (length(rows) < 2) next
    coords_im <- as.matrix(df[rows, c("x","y")])
    nn_e <- dbscan::frNN(coords_im, eps = rad)
    ct_im <- ct_all[rows]
    i_loc <- rep.int(seq_along(rows), lengths(nn_e$id))
    j_loc <- unlist(nn_e$id,   use.names = FALSE)
    d_loc <- unlist(nn_e$dist, use.names = FALSE)
    keep  <- ct_im[j_loc] != ct_im[i_loc]      ## cross-celltype only
    if (!any(keep)) next
    i_loc <- i_loc[keep]; j_loc <- j_loc[keep]; d_loc <- d_loc[keep]
    i_glb <- rows[i_loc]; j_glb <- rows[j_loc]
    ## Edge correction folds into row i (1/af_i), matching dense E_tech / af_i.
    w_glb <- exp(-d_loc / H_TECH) / af_etech[i_glb]
    ii_all[[ii_im]] <- i_glb; jj_all[[ii_im]] <- j_glb; ww_all[[ii_im]] <- w_glb
    if (ii_im %% 50 == 0) cat(sprintf("    [W] image %d/%d done\n", ii_im, length(by_image)))
  }
  ii <- unlist(ii_all, use.names = FALSE)
  jj <- unlist(jj_all, use.names = FALSE)
  ww <- unlist(ww_all, use.names = FALSE)
  rm(ii_all, jj_all, ww_all)
  ## sparseMatrix sums duplicate (i,j) entries; there is one entry per neighbour
  ## pair here, so the sum is the exact kernel weight.
  ambient_W <- Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(n, n))
  ambient_W <- methods::as(ambient_W, "CsparseMatrix")
  say("    W: %d x %d, %d nnz (%.1f nnz/row), built in %.1f s",
      nrow(ambient_W), ncol(ambient_W), length(ww), length(ww) / n,
      as.numeric(difftime(Sys.time(), t_et, units = "secs")))
  rm(ii, jj, ww)
  gc_peak("after W build")

  ## ---- VALIDATION: spot-check ambient_W %*% Y == dense E_tech on a few genes.
  ## This reproduces the dense Mel builder's E^tech construction EXACTLY
  ## (per-image frNN, cross-celltype, exp(-d/h_tech), divided by af_i) and
  ## confirms the streaming identity before the (expensive) full fit. ----
  set.seed(1)
  spot_g <- sort(sample(ncol(Y), min(5L, ncol(Y))))
  E_spot <- matrix(0, n, length(spot_g))
  for (ii_im in seq_along(by_image)) {
    rows <- by_image[[ii_im]]; if (length(rows) < 2) next
    coords_im <- as.matrix(df[rows, c("x","y")])
    nn_e <- dbscan::frNN(coords_im, eps = rad)
    Y_im <- Y[rows, spot_g, drop = FALSE]; ct_im <- ct_all[rows]
    for (i in seq_along(rows)) {
      nbr <- nn_e$id[[i]]; if (!length(nbr)) next
      dl  <- nn_e$dist[[i]]
      kk  <- ct_im[nbr] != ct_im[i]; if (!any(kk)) next
      w   <- exp(-dl[kk] / H_TECH)
      E_spot[rows[i], ] <- as.numeric(crossprod(w, Y_im[nbr[kk], , drop = FALSE]))
    }
    E_spot[rows, ] <- E_spot[rows, , drop = FALSE] / af_etech[rows]
  }
  a_spot <- as.matrix(ambient_W %*% Y[, spot_g, drop = FALSE])
  max_id <- max(abs(a_spot - E_spot))
  say("    [W identity check] max|W%%*%%Y - dense E_tech| over %d genes = %.3e",
      length(spot_g), max_id)
  if (max_id > 1e-8)
    stop("W identity check FAILED (max abs diff ", max_id, ") -- sparse ambient ",
         "does not reproduce dense E_tech. Aborting.")
  rm(E_spot, a_spot)

  ambient_image_idx <- as.integer(df$imageID)
  ambient_n_images  <- nlevels(df$imageID)
}

## ---- Sparse Y for the streaming solver (dense Y kept for downstream MCSD) ----
Y_sparse <- methods::as(methods::as(methods::as(Y, "dMatrix"), "generalMatrix"),
                        "CsparseMatrix")
if (!is.null(colnames(Y))) colnames(Y_sparse) <- colnames(Y)

## ============================================================================
## 5. data-informed tau weights (ported verbatim; optional).
## ============================================================================
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

## ============================================================================
## 5b. Streaming joint fit (the ONLY model-machinery change vs the dense builder).
## ============================================================================
PQL_THREADS <- as.integer(Sys.getenv("R_PQL_THREADS", unset = "8"))
## Guarded speed-mode knobs (defaults match build_bc_streaming.R).
ALPHA_WARMUP   <- as.numeric(Sys.getenv("R_ALPHA_WARMUP",   unset = "6"))
EARLY_STOP_TOL <- as.numeric(Sys.getenv("R_EARLY_STOP_TOL", unset = "2e-2"))
MIN_ITER       <- as.integer(Sys.getenv("R_MIN_ITER",       unset = "12"))
N_ITER         <- as.integer(Sys.getenv("R_N_ITER",         unset = "32"))
CHUNK_SIZE     <- as.integer(Sys.getenv("R_CHUNK_SIZE",     unset = "128"))
DISP_MODEL     <- Sys.getenv("R_DISP_MODEL", unset = "nb2")
TAU_SHRINK     <- Sys.getenv("R_TAU",        unset = "adaptive")
FUSE_RHO       <- nzchar(Sys.getenv("R_FUSE_RHO"))
say("[5/7] fit_pace_mvpql_streaming (%d C++ threads; n_iter=%d disp=%s tau=%s; alpha_warmup=%g early_stop_tol=%g min_iter=%d fuse_rho=%s) ...",
    PQL_THREADS, N_ITER, DISP_MODEL, TAU_SHRINK,
    ALPHA_WARMUP, EARLY_STOP_TOL, MIN_ITER, FUSE_RHO)
gc_peak("before fit")
t_fit <- Sys.time()
fit <- fit_pace_mvpql_streaming(
  Y = Y_sparse, X_fixed = X_fixed, df = df, re_specs = re_specs,
  offset_vec = offset_vec,
  data_informed_W = data_informed_W,
  ambient_W          = ambient_W,
  ambient_image_idx  = ambient_image_idx,
  ambient_n_images   = ambient_n_images,
  bleed_percell       = TRUE,
  percell_anchor_mask = percell_anchor_mask,
  percell_anchor_idx  = percell_anchor_idx,
  disp_model    = DISP_MODEL,
  tau_shrinkage = TAU_SHRINK,
  n_iter = N_ITER, tol = 5e-3,
  BPPARAM = SerialParam(), n_threads = PQL_THREADS,
  interior_precision = 1L,
  chunk_size     = CHUNK_SIZE,
  alpha_warmup   = ALPHA_WARMUP,
  early_stop_tol = EARLY_STOP_TOL,
  min_iter       = MIN_ITER,
  fuse_rho       = FUSE_RHO,
  return_mu      = TRUE,        ## need full mu + technical_offset for decomposition / MCSD
  verbose        = TRUE)
say("    fit done in %.1f min", as.numeric(difftime(Sys.time(), t_fit, units = "mins")))
gc_peak("after fit")

if (nzchar(Sys.getenv("R_RAW_FIT_ONLY"))) {
  dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(fit = fit, gene_set = colnames(Y)), OUT)
  say("Saved RAW FIT ONLY to %s", OUT); quit(save = "no")
}

## ============================================================================
## 6. mashr + decomposition (ported from the dense Mel builder; resp_term set).
## ============================================================================
say("[6/7] mashr + decomposition (resp_term=%s) ...", RESP_TERM)
results_mv <- mvpql_to_results_multi(fit, keep_block = "celltype")
shrunken_long <- apply_mashr_shrinkage(
  results = results_mv, focals = TYPES, neighbours = TYPES,
  resp_term = RESP_TERM)
say("    %s shrunk rows", format(nrow(shrunken_long), big.mark = ","))
dec <- tryCatch(
  mvpql_variance_decomposition_multi(
    fit = fit, df = df, Y = Y,
    vars = TYPES, X_fixed = X_fixed,
    resp_term = RESP_TERM, focal_levels = TYPES),
  error = function(e) {
    say("WARNING: decomposition failed (%s); saving NULL", conditionMessage(e))
    NULL
  })

## ============================================================================
## 7. MCSD per pair (collapsed-7 Mel pairs, Responder handling; ported verbatim).
## ============================================================================
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

## ---- cell_meta (same columns as the dense Mel canonical) ----
cell_meta <- data.frame(cell_id = as.character(seq_len(nrow(df))),
                        nCount = df$nCount,
                        Condition = as.character(df$Responder),
                        patientID = as.character(df$imageID),
                        imageID = as.character(df$imageID),
                        celltype = as.character(df$celltype),
                        stringsAsFactors = FALSE)

## ---- intrinsic patient-level disease test (ported verbatim) ----
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

## ============================================================================
## Save (same list schema as the dense Mel canonical; NON-canonical path).
## ============================================================================
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = fit, shrunken_long = shrunken_long, gene_set = colnames(Y),
              decomposition = dec, mcsd_canonical = mcsd_canonical,
              twostage_disease = twostage_disease,
              h_tech = H_TECH, h_bio = H_BIO,
              K_tech = K_tech, K_bio = K_bio, cell_meta = cell_meta,
              resolution = "kernel_refactor_collapsed7_streaming"),
        OUT)
say("Saved %s", OUT)
say("TOTAL: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

## pace_intrinsic_state.R
## SENSITIVITY / DIAGNOSTIC (not canonical): SVCA-style intrinsic cell-state confound check.
## The canonical "Spatial cell state" block is an associational variance share (upper bound);
## this quantifies how much of it is SHARED with the cell's own intrinsic state. It is scale-mixed
## (NB-link spatial vs OLS observed intrinsic), commonality-based (can give negative shared), and
## co-regulation confounded ("intrinsic" includes possibly-spatially-induced program genes), so it
## supports the upper-bound caveat on the canonical spatial blocks; it is NOT a headline decomposition.
##
## For each (focal cell type, gene) it splits the within-focal variance of observed
## log-CP10k expression into five buckets that sum to 100%:
##   intrinsic (unique) | spatial (unique) | shared (state x spatial) | spillover | residual
##
## - spatial   = PACE's FITTED neighbour-slope contribution (fit$U), entered as ONE fixed
##               projection so spatial stays = PACE's estimate (no re-fit of slopes).
## - spillover = PACE's FITTED technical offset (fit$technical_offset_mat).
## - intrinsic = top-`n_pc` PCs of the cell's own transcriptional state, built LEAVE-GENE-OUT
##               (gene g removed from the PCA before the eigendecomposition AND the scores),
##               so a gene never predicts itself.
## - unique / shared = commonality analysis between the intrinsic and spatial sets.
##
## The GLMM fit is NOT modified: this is a post-hoc decomposition on the existing fit, which
## avoids the self-leak that shared PCs-in-the-fit would cause for high-loading genes.
##
## Efficiency: the focal-cell gene-gene covariance is formed ONCE per focal; the leave-gene-out
## PCA for each gene is then a symmetric submatrix eigendecomposition (drop g's row/col).

## fraction of var(y) explained by an OLS fit on design X (intercept added)
.var_explained <- function(y, X) {
  if (is.null(X)) return(0)
  f <- stats::lm.fit(cbind(1, X), y)
  1 - sum(f$residuals^2) / sum((y - mean(y))^2)
}

## five-bucket commonality decomposition for one gene given response y and the three covariate
## sets (intrinsic PCs, PACE spatial projection, PACE spillover offset)
.commonality_five <- function(y, INT, SPAT, SPILL) {
  r_S   <- .var_explained(y, SPAT)
  r_I   <- .var_explained(y, INT)
  r_SI  <- .var_explained(y, cbind(SPAT, INT))
  r_all <- .var_explained(y, cbind(SPAT, INT, SPILL))
  c(intrinsic_uniq = r_SI - r_S,
    spatial_uniq   = r_SI - r_I,
    shared         = r_I + r_S - r_SI,
    spillover      = r_all - r_SI,
    residual       = 1 - r_all,
    marg_spatial   = r_S,
    marg_intrinsic = r_I)
}

## Main entry point.
##   fit        : PACE fit object (needs $U, $re_meta$Z, $technical_offset_mat)
##   Y          : cells x genes count matrix (colnames = genes; rows aligned to celltype/nCount)
##   celltype   : per-cell cell-type vector (length nrow(Y))
##   nCount     : per-cell library size (length nrow(Y))
##   focal_levels : cell types to decompose (default: all)
##   n_pc       : number of leave-gene-out state PCs (default 5)
##   min_cells  : skip focals with fewer cells than this
##   min_nonzero: a gene must be non-zero in > this many focal cells to enter the PCA basis
## Returns a tibble: focal, gene, and the five buckets (%) + marginal R2 diagnostics (%).
pace_intrinsic_state_decomp <- function(fit, Y, celltype, nCount,
                                        focal_levels = NULL, n_pc = 5L,
                                        min_cells = 50L, min_nonzero = 50L) {
  stopifnot(!is.null(fit$technical_offset_mat))
  genes <- colnames(Y)
  if (is.null(focal_levels)) focal_levels <- sort(unique(as.character(celltype)))
  Z_design <- fit$re_meta$Z
  U_blup   <- fit$U
  tech_off <- fit$technical_offset_mat

  out <- vector("list", length(focal_levels))
  for (fi in seq_along(focal_levels)) {
    focal <- focal_levels[fi]
    focal_cells <- which(as.character(celltype) == focal)
    if (length(focal_cells) < min_cells) next

    intercept_col <- paste0(focal, "::(Intercept)")
    if (!(intercept_col %in% colnames(Z_design))) next
    spatial_cols <- grep(paste0("^", focal, "::"), colnames(Z_design), value = TRUE)
    spatial_cols <- setdiff(spatial_cols, intercept_col)
    N_focal <- as.matrix(Z_design[focal_cells, spatial_cols, drop = FALSE])

    ## log-CP10k for this focal (cells x genes), centered; PCA basis = expressed genes
    lib <- nCount[focal_cells]
    Ecp <- as.matrix(log1p(Y[focal_cells, , drop = FALSE] *
                             (1e4 / lib)))                       # cells x genes
    nz  <- colSums(Ecp != 0)
    basis <- genes[nz > min_nonzero]
    M  <- scale(Ecp[, basis, drop = FALSE], center = TRUE, scale = FALSE)  # cells x basis
    Cf <- crossprod(M) / (nrow(M) - 1)                          # basis x basis covariance (once)
    basis_idx <- setNames(seq_along(basis), basis)

    rows <- vector("list", length(genes))
    for (gi in seq_along(genes)) {
      g <- genes[gi]
      ## skip genes with no within-focal variance (e.g. zero counts across all focal cells):
      ## there is nothing to decompose and var(y)=0 would give 0/0 buckets.
      if (stats::var(Ecp[, g]) < 1e-12) next
      y <- Ecp[, g] - mean(Ecp[, g])
      ## leave-gene-out top-n_pc state PCs
      if (g %in% basis) {
        j <- basis_idx[[g]]
        Csub <- Cf[-j, -j, drop = FALSE]
        V <- eigen(Csub, symmetric = TRUE)$vectors[, seq_len(n_pc), drop = FALSE]
        INT <- M[, -j, drop = FALSE] %*% V
      } else {
        V <- eigen(Cf, symmetric = TRUE)$vectors[, seq_len(n_pc), drop = FALSE]
        INT <- M %*% V
      }
      SPAT  <- as.numeric(N_focal %*% U_blup[spatial_cols, g])
      SPILL <- as.numeric(tech_off[focal_cells, g])
      b <- .commonality_five(y, INT, SPAT, SPILL)
      rows[[gi]] <- c(focal = focal, gene = g, b)
    }
    out[[fi]] <- do.call(rbind, rows)
  }

  df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
  num <- c("intrinsic_uniq","spatial_uniq","shared","spillover","residual",
           "marg_spatial","marg_intrinsic")
  for (col in num) df[[col]] <- 100 * as.numeric(df[[col]])
  tibble::tibble(
    focal = df$focal, gene = df$gene,
    `Intrinsic (uniq) %` = df$intrinsic_uniq,
    `Spatial (uniq) %`   = df$spatial_uniq,
    `Shared %`           = df$shared,
    `Spillover %`        = df$spillover,
    `Residual %`         = df$residual,
    `marg.spatial R2 %`  = df$marg_spatial,
    `marg.intrinsic R2 %`= df$marg_intrinsic,
    n_pc = as.integer(n_pc))
}

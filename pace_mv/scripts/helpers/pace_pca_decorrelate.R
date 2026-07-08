## pace_pca_decorrelate.R --- PCA pre-decorrelation of the K_<neighbour>
## random-slope design (Hodges-Reich 2010 Restricted Spatial Regression).
##
## RATIONALE.  PACE-MV's K_<neighbour> columns (Gaussian-kernel densities of
## neighbour celltypes around each cell) are highly correlated (pairwise r ≈
## 0.5–0.8 in dense regions: "shared density mode" dominates).  Under PQL,
## correlated random slopes produce coefficient blowups in opposite directions
## that sum to a small fitted slope but appear inflated individually.  Hodges &
## Reich 2010 ("Adding spatially-correlated errors can mess up the fixed effect
## you love"; _Amer. Statistician_ 64:325) call this **spatial confounding** and
## prescribe orthogonalisation of the spatial predictors before fitting (the
## Restricted Spatial Regression construction).
##
## METHOD.  We compute a single shared eigen-basis V of K'K across all cells and
## reparameterise the K block as K·V (orthogonal columns by construction).  The
## random-slope fit operates on PCs (uncorrelated → no swap inflation).  Post-
## fit, BLUPs in PC basis γ are back-transformed to the original (focal, neighbour)
## basis via β = V γ.  Pratt 1987 variance decomposition is invariant under
## orthogonal predictor transforms (V'V = I) — focal-total V_RxS is unchanged;
## per-pair shares are exact linear combinations of PC loadings.
##
## ENV HOOK: set R_K_PCA=1 in the builder env to activate.

#' Compute a shared orthonormal basis V for the K_<neighbour> columns and return
#' the PCA-rotated columns K·V plus the V matrix.
#'
#' @param df data.frame; must contain the columns named in `types`
#' @param types character vector of column names (e.g. neighbour celltype names)
#' @param center  logical; subtract column means before PCA (recommended TRUE)
#' @param scale   logical; divide by column SDs before PCA (FALSE preserves the
#'                kernel scale information; TRUE makes all neighbour densities
#'                equally weighted regardless of celltype abundance)
#' @return list with components:
#'   - K_pc: n × length(types) matrix of PCA-rotated columns
#'   - V:    length(types) × length(types) eigenvector matrix; column j is
#'           the j-th PC loading vector
#'   - center: per-column means (length(types))
#'   - scale:  per-column scale factors (length(types); 1 if !scale)
#'   - pc_names: c("PC1", ..., paste0("PC", length(types)))
pace_pca_decorrelate_K <- function(df, types,
                                      center = TRUE, scale = FALSE) {
  stopifnot(all(types %in% colnames(df)))
  K <- as.matrix(df[, types, drop = FALSE])
  storage.mode(K) <- "double"
  cm <- if (center) colMeans(K, na.rm = TRUE) else rep(0, ncol(K))
  csd <- if (scale)  apply(K, 2, stats::sd, na.rm = TRUE) else rep(1, ncol(K))
  csd[!is.finite(csd) | csd < 1e-12] <- 1
  K_cs <- sweep(K, 2L, cm,  "-")
  K_cs <- sweep(K_cs, 2L, csd, "/")
  ## Eigen of K'K / n; columns of V are loadings
  S <- crossprod(K_cs) / nrow(K_cs)
  eig <- eigen(S, symmetric = TRUE)
  V   <- eig$vectors
  pc_names <- paste0("PC", seq_along(types))
  rownames(V) <- types
  colnames(V) <- pc_names
  K_pc <- K_cs %*% V
  colnames(K_pc) <- pc_names
  list(K_pc = K_pc, V = V, center = cm, scale = csd,
       pc_names = pc_names, types = types,
       eigenvalues = eig$values)
}

#' Back-transform U BLUPs from PC basis to original K_<neighbour> basis.
#'
#' PACE's U matrix has rows like "<focal>::PC1", "<focal>::PC2", ...,
#' "<focal>::ResponderPD:PC1", ...  This function reconstructs the rows in
#' the original neighbour basis: "<focal>::<neighbour>",
#' "<focal>::ResponderPD:<neighbour>".  Only rows that match the PC pattern
#' are transformed; other rows (e.g. "<focal>::(Intercept)", imageID rows)
#' are passed through unchanged.
#'
#' @param U          q × G matrix of BLUPs in PC basis
#' @param V          length(types) × length(types) loading matrix from
#'                   pace_pca_decorrelate_K
#' @param types      original neighbour celltype names
#' @param focals     character vector of focal celltype names (rows starting
#'                   with "<focal>::" containing PCs)
#' @return list with components U_orig (q' × G; q' = q for matched rows, +
#'         passthrough for unmatched) and rownames mapped to original basis.
pace_pca_backtransform_U <- function(U, V, types, focals) {
  rn <- rownames(U)
  pc_names <- colnames(V)
  if (is.null(pc_names)) pc_names <- paste0("PC", seq_len(ncol(V)))
  G <- ncol(U)
  q <- nrow(U)
  ## Build a map from PC-basis rows to original-basis rows per (focal, term-prefix)
  ## A "term-prefix" is the optional "ResponderPD:" (or other interaction prefix)
  ## that appears between the focal-celltype name and the PC name.
  pc_regex <- paste0("^(?<focal>", paste(focals, collapse = "|"), ")::",
                     "(?<prefix>([A-Za-z0-9_]+:)*)(?<pc>", paste(pc_names, collapse = "|"), ")$")
  parsed <- regmatches(rn, regexec(pc_regex, rn, perl = TRUE))
  has_match <- sapply(parsed, length) > 0
  ## Group rows by (focal, prefix) and apply V to convert PC-cols to types
  U_orig_rows <- list()
  out_rn <- character(0)
  for (focal in focals) {
    for (prefix in unique(sapply(parsed[has_match], `[`, "prefix"))) {
      idx_block <- which(has_match &
        sapply(parsed, function(p) length(p) > 0 && p["focal"] == focal &&
                                    p["prefix"] == prefix))
      if (!length(idx_block)) next
      ## Order block by PC name
      pc_of_idx <- sapply(parsed[idx_block], `[`, "pc")
      ord <- match(pc_names, pc_of_idx)
      ok <- !is.na(ord)
      idx_ord <- idx_block[ord[ok]]
      if (length(idx_ord) < length(pc_names)) next  ## incomplete block; skip
      ## U[idx_ord, ] is length(pc_names) × G in PC basis
      U_pc <- U[idx_ord, , drop = FALSE]
      ## Back-transform: β = V γ (per gene).  V is types × pc_names.
      U_neighbour <- V %*% U_pc
      new_rn <- paste0(focal, "::", prefix, types)
      for (j in seq_along(new_rn)) {
        U_orig_rows[[new_rn[j]]] <- U_neighbour[j, ]
      }
      out_rn <- c(out_rn, new_rn)
    }
  }
  ## Combine: passthrough rows + back-transformed rows
  pass_idx <- which(!has_match)
  U_pass <- if (length(pass_idx)) U[pass_idx, , drop = FALSE] else NULL
  if (length(out_rn)) {
    U_back <- do.call(rbind, U_orig_rows[out_rn])
    if (!is.null(U_pass)) {
      U_out <- rbind(U_pass, U_back)
    } else {
      U_out <- U_back
    }
  } else {
    U_out <- U_pass
  }
  colnames(U_out) <- colnames(U)
  U_out
}

#' Build per-focal V matrices when each focal celltype gets its OWN PCA basis.
#'
#' Per-focal decorrelation can yield better identifiability than shared PCA
#' (each focal's K-correlation structure may differ), at the cost of needing
#' to rotate the Z-block columns differently for each focal.  PACE's current
#' formula-based RE construction shares ONE set of column names across focals,
#' so per-focal decorrelation requires changing the Z-build path.  For now
#' we ship the SHARED-basis version (one V across cells) as the default and
#' return the per-focal V matrices for users who want to compare.
#'
#' @param df data.frame
#' @param celltype_col name of the celltype column
#' @param types neighbour celltype names (= K column names in df)
#' @param focals focal celltype names (subset of unique(df[[celltype_col]]))
#' @return named list, one V matrix per focal
pace_pca_per_focal_V <- function(df, celltype_col, types, focals) {
  out <- list()
  for (focal in focals) {
    rows <- which(df[[celltype_col]] == focal)
    if (length(rows) < 10) {
      out[[focal]] <- NULL
      next
    }
    res <- pace_pca_decorrelate_K(df[rows, , drop = FALSE], types,
                                     center = TRUE, scale = FALSE)
    out[[focal]] <- res$V
  }
  out
}

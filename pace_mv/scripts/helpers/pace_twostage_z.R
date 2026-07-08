## pace_twostage_z.R — INTRINSIC patient-level disease significance test for PACE.
## For each (focal, neighbour, gene), compute per-patient OLS slope of normalized expression
## on PACE's focal::neighbour kernel, then Welch-t-test slopes between case vs ref patients.
## Uses log1p(y/libsize·scale) (NOT y-μ̂ residual — that would subtract off PACE's own
## shrunk estimate of the very effect we want to test). PACE provides the focal-cell mask,
## the bio-kernel proximity, and the per-cell libsize; the per-patient OLS+ t-test gives
## patient-level inference — the correct effective n for disease contrasts at low patient n.
## Replaces the post-hoc audit_disease_cohort.R diagnostic with a built-in PACE output that
## tracks the cluster-confirmed verdict (recovers liver fibrosis hits; rejects DKD/Mel null).
##
## Usage:
##   res <- pace_twostage_z(fit, raw_Y, cluster_vec, cond_vec, case_level,
##                          libsize = NULL, focal_only = NULL, min_cells_per_pat = 20,
##                          min_pat_per_arm = 3)
## Returns a data.frame: focal, neighbour, gene, b_case, b_ref, twostage_t, twostage_p,
##   n_case_pat, n_ref_pat, q_pair (BH within-pair q-value).

pace_twostage_z <- function(fit, raw_Y, cluster_vec, cond_vec, case_level,
                             libsize = NULL, focal_only = NULL,
                             min_cells_per_pat = 20, min_pat_per_arm = 3) {
  stopifnot(nrow(raw_Y) == length(cluster_vec), length(cluster_vec) == length(cond_vec))
  Z <- fit$re_meta$Z; mu <- fit$mu
  gn <- colnames(raw_Y); n_genes <- length(gn)
  TYPES <- fit$re_meta$blocks[[1]]$group_levels
  if (is.null(libsize)) libsize <- Matrix::rowSums(raw_Y)
  if (is.null(focal_only)) focal_only <- TYPES

  out <- list()
  for (fc in focal_only) {
    fc_int <- paste0(fc, "::(Intercept)")
    if (!(fc_int %in% colnames(Z))) next
    cells <- which(as.numeric(Z[, fc_int]) != 0)
    if (length(cells) < 100) next
    cl_fc <- cluster_vec[cells]; cd_fc <- cond_vec[cells]
    pats <- unique(cl_fc); n_pats <- length(pats)
    if (n_pats < 2*min_pat_per_arm) next

    ## per-cell normalized log expression: ln_ig = log1p(y/libsize · scale)
    ## scale = median(libsize) so all cells live on a comparable expression scale
    Y_fc  <- as.matrix(raw_Y[cells, , drop = FALSE])
    lib_fc <- pmax(libsize[cells], 1)
    scale_fac <- exp(median(log(lib_fc)))
    R_fc <- log1p(sweep(Y_fc, 1, scale_fac / lib_fc, "*"))

    ## majority condition per patient (cells of focal type only)
    patcond <- setNames(rep(NA_character_, n_pats), pats)
    pat_n   <- setNames(integer(n_pats),         pats)
    for (p in pats) { k <- which(cl_fc == p); tb <- table(cd_fc[k])
      patcond[p] <- names(tb)[which.max(tb)]; pat_n[p] <- length(k) }
    keep_pat <- which(pat_n >= min_cells_per_pat)
    if (length(keep_pat) < 2*min_pat_per_arm) next
    case_idx <- keep_pat[patcond[keep_pat] == case_level]
    ref_idx  <- keep_pat[patcond[keep_pat] != case_level & !is.na(patcond[keep_pat])]
    if (length(case_idx) < min_pat_per_arm || length(ref_idx) < min_pat_per_arm) next

    for (nb in setdiff(TYPES, fc)) {
      nb_col <- paste0(fc, "::", nb)
      if (!(nb_col %in% colnames(Z))) next
      prox <- as.numeric(Z[cells, nb_col])
      if (var(prox, na.rm = TRUE) == 0) next

      ## per-patient OLS slope of R on prox (centered): b_pg = Σ pc·R / Σ pc²
      B <- matrix(NA_real_, nrow = length(pats), ncol = n_genes,
                  dimnames = list(pats, gn))
      for (p in pats) {
        k <- which(cl_fc == p); if (length(k) < min_cells_per_pat) next
        pc <- prox[k] - mean(prox[k]); ss <- sum(pc * pc); if (ss <= 0) next
        B[p, ] <- as.numeric(crossprod(pc, R_fc[k, , drop = FALSE])) / ss
      }

      ## per-gene Welch t between case and ref patient slopes (vectorized)
      Bc <- B[case_idx, , drop = FALSE]; Br <- B[ref_idx, , drop = FALSE]
      mc <- colMeans(Bc, na.rm = TRUE);   mr <- colMeans(Br, na.rm = TRUE)
      vc <- apply(Bc, 2, var, na.rm = TRUE); vr <- apply(Br, 2, var, na.rm = TRUE)
      nc <- colSums(is.finite(Bc));         nr <- colSums(is.finite(Br))
      se <- sqrt(vc / nc + vr / nr)
      tt <- (mc - mr) / se
      df_w <- (vc / nc + vr / nr)^2 / ((vc / nc)^2 / pmax(nc - 1, 1) + (vr / nr)^2 / pmax(nr - 1, 1))
      pv <- 2 * pt(-abs(tt), df = pmax(df_w, 1))

      d <- data.frame(focal = fc, neighbour = nb, gene = gn,
                      b_case = mc, b_ref = mr, twostage_t = as.numeric(tt),
                      twostage_p = as.numeric(pv),
                      n_case_pat = as.integer(nc), n_ref_pat = as.integer(nr),
                      stringsAsFactors = FALSE)
      d$q_pair <- p.adjust(d$twostage_p, "BH")
      out[[paste(fc, nb, sep = "::")]] <- d
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

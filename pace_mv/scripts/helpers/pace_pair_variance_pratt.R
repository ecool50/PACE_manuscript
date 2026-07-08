## pace_pair_variance_pratt.R — Pratt 1987 signed pair-level variance decomposition.
## V_pair(c, t, g) = u_{c,t,g} · (Σ_K(c) u_{c,g})_t
## Sums exactly to V_spatial(c, g) = u' Σ_K(c) u.  O(K^2) per (focal, gene).
## Standard pair attribution under the Pratt / linear-case CGP12 / Mara-Tarantola
## framework for variance of a linear form in correlated inputs.

pace_pair_variance_pratt <- function(mv, cond_prefix = NULL, focals = NULL,
                                     cohort_label = "cohort", block_label = NULL) {
  fit <- mv$fit
  Z <- fit$re_meta$Z
  gn <- mv$gene_set
  colnames(fit$U) <- gn
  TYPES <- if (!is.null(fit$re_meta$blocks))
             fit$re_meta$blocks[[1]]$group_levels
           else fit$re_meta$group_levels
  if (is.null(focals)) focals <- TYPES
  if (is.null(block_label))
    block_label <- if (is.null(cond_prefix)) "Spatial" else "RxS"

  pair_rows <- list()
  focal_rows <- list()

  for (fc in focals) {
    fc_int <- paste0(fc, "::(Intercept)")
    if (!(fc_int %in% colnames(Z))) next
    cells <- which(as.numeric(Z[, fc_int]) != 0)
    if (length(cells) < 50) next
    term_names <- if (is.null(cond_prefix))
                    paste0(fc, "::", TYPES)
                  else paste0(fc, "::", cond_prefix, ":", TYPES)
    keep <- term_names %in% colnames(Z) & term_names %in% rownames(fit$U)
    if (!any(keep)) next
    tn <- term_names[keep]
    tt <- TYPES[keep]
    Z_fc <- as.matrix(Z[cells, tn, drop = FALSE])
    Sigma_K <- cov(Z_fc)
    U_c <- as.matrix(fit$U[tn, , drop = FALSE]); U_c[!is.finite(U_c)] <- 0
    # K × G matrix product Σ_K · U_c, then elementwise U_c * (Σ_K U_c) and sum cols → K-vector
    SU <- Sigma_K %*% U_c                                # K × G
    V_pair_gene <- U_c * SU                              # K × G : element t,g = u_{t,g}·(Σ u_g)_t
    V_pair_t <- as.numeric(rowSums(V_pair_gene))         # K
    V_pratt_total <- sum(V_pair_t)                        # = u'Σu summed over genes
    V_diag_t <- as.numeric(rowSums(U_c^2) * diag(Sigma_K))
    V_diag_total <- sum(V_diag_t)

    for (i in seq_along(tt)) {
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        cohort = cohort_label, block = block_label,
        focal = fc, neighbour = tt[i],
        V_pair_pratt = V_pair_t[i],
        V_pair_diag  = V_diag_t[i],
        within_focal_share_pct = 100 * V_pair_t[i] / V_pratt_total,
        within_focal_diag_pct  = 100 * V_diag_t[i] / V_diag_total,
        sign = sign(V_pair_t[i]),
        stringsAsFactors = FALSE
      )
    }
    focal_rows[[length(focal_rows) + 1]] <- data.frame(
      cohort = cohort_label, block = block_label, focal = fc,
      n_cells = length(cells),
      V_block_pratt = V_pratt_total,
      V_block_diag  = V_diag_total,
      cross_cov_pct = 100 * (V_pratt_total - V_diag_total) / V_pratt_total,
      n_negative_pairs = sum(V_pair_t < 0),
      stringsAsFactors = FALSE
    )
  }
  list(pair_long = do.call(rbind, pair_rows),
       focal_summary = do.call(rbind, focal_rows))
}

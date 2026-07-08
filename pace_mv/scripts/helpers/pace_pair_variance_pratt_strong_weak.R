## pace_pair_variance_pratt_strong_weak.R
## Strong/weak split of the canonical Pratt pair variance (Nan-Azriel-Schwartzman
## 2026 strong-localized vs weak-diffuse partition, applied on the Pratt basis).
##
## Mirrors scripts/helpers/pace_pair_variance_pratt.R exactly for the total
## V_pair (V_pair_gene = U_c * (Sigma_K %*% U_c), summed over genes), then splits
## each pair's per-gene Pratt contribution into:
##   - STRONG : genes that pass the significance gate (lfsr < lfsr_thresh) for
##              that (focal, neighbour) term  -> the driver-carried variance
##   - WEAK   : all remaining genes            -> the null/diffuse-carried variance
## STRONG + WEAK = V_pair_pratt (the canonical total), so no variance is created
## or lost; this only attributes it. The locked pace_pair_variance_pratt.R is
## NOT modified. lfsr comes from mv$shrunken_long; slopes come from fit$U (raw
## BLUPs), matching the canonical Pratt.

pace_pair_variance_pratt_strong_weak <- function(mv, cond_prefix = NULL,
                                                 focals = NULL, lfsr_thresh = 0.05,
                                                 cohort_label = "cohort") {
  fit <- mv$fit
  Z   <- fit$re_meta$Z
  gn  <- mv$gene_set
  colnames(fit$U) <- gn
  shr <- mv$shrunken_long

  TYPES <- if (!is.null(fit$re_meta$blocks)) fit$re_meta$blocks[[1]]$group_levels
           else fit$re_meta$group_levels
  if (is.null(focals)) focals <- TYPES

  pair_rows <- list()

  for (fc in focals) {
    fc_int <- paste0(fc, "::(Intercept)")
    if (!(fc_int %in% colnames(Z))) next
    cells <- which(as.numeric(Z[, fc_int]) != 0)
    if (length(cells) < 50) next

    term_names <- if (is.null(cond_prefix)) paste0(fc, "::", TYPES)
                  else paste0(fc, "::", cond_prefix, ":", TYPES)
    keep <- term_names %in% colnames(Z) & term_names %in% rownames(fit$U)
    if (!any(keep)) next
    tn <- term_names[keep]
    tt <- TYPES[keep]

    Sigma_K <- cov(as.matrix(Z[cells, tn, drop = FALSE]))
    U_c     <- as.matrix(fit$U[tn, , drop = FALSE])
    U_c[!is.finite(U_c)] <- 0
    V_pair_gene <- U_c * (Sigma_K %*% U_c)              # K x G, canonical per-gene Pratt

    for (i in seq_along(tt)) {
      nc  <- tt[i]
      vg  <- V_pair_gene[i, ]                            # per-gene Pratt for this pair
      shr_term <- if (is.null(cond_prefix)) nc else paste0(cond_prefix, ":", nc)
      bv <- shr[shr$focal == fc & shr$neighbour == nc & shr$term == shr_term, ]
      bv <- bv[!duplicated(bv$gene), ]
      strong_genes <- bv$gene[which(bv$lfsr < lfsr_thresh)]
      is_strong    <- gn %in% strong_genes

      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        cohort         = cohort_label,
        focal          = fc,
        neighbour      = nc,
        V_pair_pratt   = sum(vg),
        V_pair_strong  = sum(vg[is_strong]),
        V_pair_weak    = sum(vg[!is_strong]),
        n_strong       = length(strong_genes),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, pair_rows)
}

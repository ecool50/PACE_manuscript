## pace_pair_variance_pratt_prox.R — PROXIMITY-CORRECTED pairwise attribution (SENSITIVITY).
##
## Like the canonical Pratt pairwise (pace_pair_variance_pratt.R) but corrected for the variation in
## proximity (neighbour-count variance) so the attribution reflects per-unit-proximity EFFECT, not how
## abundant each neighbour happens to be. Two changes from canonical:
##   (1) use the neighbour-count CORRELATION matrix R (unit variance) instead of covariance Sigma_K,
##       which removes the Var(N_t) abundance weighting (Ellis "lever 2");
##   (2) use SE-STANDARDISED slopes u/se (Wald form) instead of raw slopes, which down-weights
##       unreliable rare-type slopes and prevents the rare-type blow-up that (1) alone causes
##       (e.g. Tumour<-Mast 65% with raw slopes -> 0.4% with u/se).
## V_pair*(c,t) = (u/se)_t * (R (u/se))_t, summed over genes. Reported as the RAW Pratt share
## V_pair / V_state_stdz, where V_state_stdz is the total over ALL neighbours INCLUDING the self type
## (focal::focal). The displayed off-diagonal columns therefore sum to <100% -- the remainder is the
## self-neighbour term, not displayed (same convention as the canonical pairwise; NOT renormalised to 1).
## This is a PRECISION-WEIGHTED effect attribution (signal-to-noise per unit proximity), NOT a variance
## share -- a sensitivity complement to the canonical variance heatmap. NOTE: slopes are already
## mashr-shrunk, so u/se layers precision weighting on top of shrinkage.

pace_pair_variance_pratt_prox <- function(mv, min_var = 1e-12, min_cells = 50) {
  if (!requireNamespace("Matrix", quietly = TRUE))               # sparse Z subsetting needs Matrix's S4 methods
    stop("pace_pair_variance_pratt_prox() requires the Matrix package")
  fit <- mv$fit; Z <- fit$re_meta$Z; U <- fit$U; seU <- fit$se_U
  TYPES <- sort(unique(as.character(mv$cell_meta$celltype)))
  out <- list()
  for (fc in TYPES) {
    ci <- paste0(fc, "::(Intercept)"); if (!(ci %in% colnames(Z))) next
    cells <- which(Z[, ci] != 0); if (length(cells) < min_cells) next
    tn <- paste0(fc, "::", TYPES)                                    # ALL neighbours INCLUDING self
    tn <- tn[tn %in% colnames(Z) & tn %in% rownames(U)]
    Zfc <- as.matrix(Z[cells, tn, drop = FALSE]); v <- apply(Zfc, 2, var); keep <- v > min_var
    tn <- tn[keep]; Zfc <- Zfc[, keep, drop = FALSE]; if (ncol(Zfc) < 2) next
    nbr <- sub(paste0("^", fc, "::"), "", tn)
    R  <- cor(Zfc)                                                   # (1) correlation -> Var(N) removed
    Uz <- U[tn, , drop = FALSE] / pmax(seU[tn, , drop = FALSE], 1e-9)# (2) SE-standardised slopes (u/se)
    Uz[!is.finite(Uz)] <- 0
    Vpair <- rowSums(Uz * (R %*% Uz))                               # summed over genes (incl self)
    tot <- sum(Vpair)                                              # V_state_stdz over ALL neighbours
    df <- data.frame(focal = fc, neighbour = nbr,
                     V_pair_stdz = Vpair,
                     raw_share_pct = if (tot != 0) 100 * Vpair / tot else NA_real_,  # V_pair/V_state, RAW
                     n_exposed = colSums(Zfc > 0), stringsAsFactors = FALSE)
    out[[fc]] <- df[df$neighbour != fc, ]                          # display off-diagonal (self in denom)
  }
  do.call(rbind, out)
}

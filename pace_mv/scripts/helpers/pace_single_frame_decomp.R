## pace_single_frame_decomp.R — CANONICAL single-frame per-(gene,focal) decomposition.
##
## One additive frame (sums to 100%) holding cell type AND the within-cell-type blocks together,
## following Elijah's formula. Magnitudes come from OBSERVED expression (log1p CP10k); the within-part
## is split by the FIT's locked within-cell-type proportions (no re-fit lm -> respects the
## "PACE spatial from the fit" rule, [[feedback_pace_model_outputs_only_2026-05-31]]).
##
## For focal cell type c and gene g, with y_i = log1p(CP10k) observed expression:
##   denominator  = sum_{i in c} (y_i - global_mean_g)^2                      [observed]
##   Cell type    = n_c * (focal_mean_{c,g} - global_mean_g)^2               [observed lineage shift]
##   within-total = sum_{i in c} (y_i - focal_mean_{c,g})^2                  [observed]
##     Spatial    = within-total * V_state / (V_state+V_spill+V_disp)        [fit proportions]
##     Spillover  = within-total * V_spill / (...)                           [fit proportions]
##     Residual   = within-total * V_disp  / (...)                           [fit proportions]
## global_mean_g = mean over ALL cells of y; focal_mean = mean over focal cells.
##
## Cell type is the observed group-mean deviation from the tissue mean (a between-vs-global SS); the
## within blocks reuse the locked gene_focal_4block decomposition. Cell-type % is scale-dependent
## (computed on log1p CP10k; raw counts and link scale differ).

single_frame_decomp_obs <- function(Y, celltype, nCount, gene_focal_block) {
  genes <- colnames(Y); ct <- as.character(celltype); TYPES <- sort(unique(ct))
  Ylog <- log1p(Y * (1e4 / nCount))                      # n x G observed log1p CP10k
  global_mean <- as.numeric(Matrix::colMeans(Ylog)); names(global_mean) <- genes
  gf <- gene_focal_block
  # Condition cohorts split the spatial component into baseline + responder
  # (gene_focal_5block has V_state_responder); no-condition cohorts have V_state.
  has_resp <- "V_state_responder" %in% names(gf)

  out <- vector("list", length(TYPES)); k <- 0L
  for (fc_type in TYPES) {
    fc <- which(ct == fc_type); n_fc <- length(fc); if (n_fc < 5) next
    Yf <- as.matrix(Ylog[fc, , drop = FALSE])
    focal_mean <- colMeans(Yf)
    SS_within  <- colSums(sweep(Yf, 2, focal_mean, "-")^2)     # per gene, observed
    SS_lineage <- n_fc * (focal_mean - global_mean)^2          # per gene, observed
    den <- SS_lineage + SS_within
    ## fit within-proportions (locked gene_focal block) for this focal
    g4 <- gf[gf$focal == fc_type, , drop = FALSE]
    idx <- match(genes, g4$gene)
    vstate <- if (has_resp) g4$V_state_baseline[idx] else g4$V_state[idx]
    vresp  <- if (has_resp) g4$V_state_responder[idx] else 0
    vspill <- g4$V_spill[idx]; vdisp <- g4$V_disp[idx]
    vt <- vstate + vresp + vspill + vdisp
    p_sp   <- ifelse(vt > 0, vstate / vt, 0)
    p_resp <- ifelse(vt > 0, vresp  / vt, 0)
    p_bl   <- ifelse(vt > 0, vspill / vt, 0)
    p_rs   <- ifelse(vt > 0, vdisp  / vt, 0)
    k <- k + 1L
    row <- data.frame(
      focal = fc_type, gene = genes,
      `Cell type %` = 100 * SS_lineage       / den,
      `Spatial %`   = 100 * SS_within * p_sp / den,
      check.names = FALSE, stringsAsFactors = FALSE)
    if (has_resp) row[["Responder spatial %"]] <- 100 * SS_within * p_resp / den
    row[["Spillover %"]] <- 100 * SS_within * p_bl / den
    row[["Residual %"]]  <- 100 * SS_within * p_rs / den
    row$SS_lineage <- SS_lineage; row$SS_within <- SS_within
    row$denom <- den; row$n_focal <- n_fc
    out[[k]] <- row
  }
  do.call(rbind, out)
}

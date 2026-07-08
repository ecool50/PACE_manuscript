## ============================================================================
##  Canonical top-driver tables for PACE-MV manuscript (LOCKED 2026-05-13)
## ============================================================================
##  Convention:
##   - Pairs ranked by pair_pct (mean R² over focal × neighbour-density × responder).
##   - Genes within pair ranked by canonical MCSD = b̂² × spec² × focal_mean.
##   - lfsr<0.05 only (no spec cutoff, no other filters).
##   - MCSD reported as SHARE OF PAIR TOTAL — sums to 1.0 across all sig genes;
##     top-3 entries are the fraction of pair-total MCSD captured by those genes.
##   - Header: pair_pct, n_sig, top3_share.
##   - Top-3 pairs × top-3 genes per cohort.
##  Source-of-truth implementation; do NOT change formulas without updating
##  feedback_mcsd_per_pair_norm.md and project_pace_mv_pipeline_final.md.
## ============================================================================
suppressPackageStartupMessages({library(dplyr); library(tibble); library(Matrix)})
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---- Cohort RDS paths (override via R_RDS_<COHORT> env if needed) ----
default_rds <- list(
  BC  = list(rds = Sys.getenv("R_RDS_BC",
              "data/breast_cancer/sweeps/mvpql_kernel_h05_h30_det10.rds"),
              resp_term = NULL),
  Mel = list(rds = Sys.getenv("R_RDS_MEL",
              "data/simvi_melanoma/sweeps/mvpql_kernel_h05_h30_det10.rds"),
              resp_term = "ResponderPD"))

## ---- per-celltype mean expression and Tan-style specificity ----
compute_per_celltype_means <- function(fit, TYPES, cells_by_ct, gn) {
  if (is.null(colnames(fit$mu))) colnames(fit$mu) <- gn
  sapply(TYPES, function(c) {
    idx <- cells_by_ct[[c]]
    if (!length(idx)) rep(0, length(gn)) else Matrix::colMeans(fit$mu[idx, , drop=FALSE])
  })
}

compute_spec <- function(mu_per_ct, focal) {
  fm <- mu_per_ct[, focal]
  others <- mu_per_ct[, setdiff(colnames(mu_per_ct), focal), drop=FALSE]
  max_other <- apply(others, 1, max)
  fm / (fm + max_other + 1e-9)
}

## ---- core compute: per-pair pair_pct + per-gene MCSD (share-of-pair) ----
compute_pair_top_drivers <- function(mv, resp_term = NULL, top_n_genes = 3L) {
  fit <- mv$fit; Z_re <- fit$re_meta$Z
  if (!is.null(fit$re_meta$blocks)) {
    blk_idx <- which(vapply(fit$re_meta$blocks, `[[`, character(1), "group_col") == "celltype")
    TYPES <- fit$re_meta$blocks[[blk_idx]]$group_levels
    cells_by_ct <- setNames(fit$re_meta$cells_by_grp_list[[blk_idx]], TYPES)
  } else {
    TYPES <- fit$re_meta$group_levels
    cells_by_ct <- lapply(TYPES, function(c) which(Z_re[, paste0(c, "::(Intercept)")] != 0))
    names(cells_by_ct) <- TYPES
  }
  gn <- colnames(fit$U)
  mu_per_ct <- compute_per_celltype_means(fit, TYPES, cells_by_ct, gn)
  alpha_g <- pmax(fit$alpha, 0); names(alpha_g) <- gn
  shr <- mv$shrunken_long
  pair_pp <- list(); pair_genes <- list()
  for (fc in TYPES) {
    cells_c <- cells_by_ct[[fc]]
    if (length(cells_c) < 50) next  ## min cells for stable RE; was 200 (2026-05-14)
    mu_bar_g <- Matrix::colMeans(fit$mu[cells_c, , drop=FALSE])
    V_resid_g <- log(1 + (1 + alpha_g) / pmax(mu_bar_g, 1e-6))
    spec_f <- compute_spec(mu_per_ct, fc)
    focal_mean_f <- mu_per_ct[, fc]
    for (nc in TYPES) {
      if (nc == fc) next
      target_term <- if (is.null(resp_term)) nc else paste0(resp_term, ":", nc)
      bv <- shr |> dplyr::filter(focal == !!fc, neighbour == !!nc, term == !!target_term) |>
        dplyr::distinct(gene, .keep_all = TRUE)
      if (!nrow(bv)) next
      b_g <- setNames(rep(0, length(gn)), gn)
      ix <- match(bv$gene, gn); ok <- !is.na(ix)
      b_g[ix[ok]] <- bv$estimate_shrunk[ok]
      col_N <- match(paste0(fc, "::", nc), colnames(Z_re))
      if (is.na(col_N)) next
      N_t <- as.numeric(Z_re[cells_c, col_N])
      vN <- if (is.null(resp_term)) stats::var(N_t, na.rm=TRUE) else {
        col_R <- match(paste0(fc, "::", resp_term), colnames(Z_re))
        R_c <- if (is.na(col_R)) rep(0, length(cells_c)) else as.numeric(Z_re[cells_c, col_R])
        stats::var(R_c * N_t, na.rm=TRUE)
      }
      R2 <- (b_g^2 * vN) / pmax(b_g^2 * vN + V_resid_g, 1e-12)
      pp <- 100 * mean(R2, na.rm=TRUE)
      mcsd_raw <- b_g^2 * spec_f^2 * focal_mean_f
      sig_idx <- which(bv$lfsr < 0.05)
      sig_genes <- bv$gene[sig_idx]
      sig_mcsd <- mcsd_raw[match(sig_genes, gn)]
      pair_sum <- if (length(sig_mcsd) && sum(sig_mcsd, na.rm=TRUE) > 0)
        sum(sig_mcsd, na.rm=TRUE) else 1
      mcsd_share <- mcsd_raw / pair_sum
      df <- tibble::tibble(gene = gn, R2_RxS = R2, b_shrunk = b_g,
                            spec = spec_f, focal_mean = focal_mean_f,
                            MCSD_raw = mcsd_raw, MCSD = mcsd_share) |>
        dplyr::inner_join(bv |> dplyr::select(gene, lfsr), by="gene") |>
        dplyr::filter(lfsr < 0.05) |>
        dplyr::arrange(dplyr::desc(MCSD)) |>
        head(top_n_genes) |>
        dplyr::mutate(sign = ifelse(b_shrunk > 0, "+", "-"),
                      R2_pct = round(R2_RxS*100, 2),
                      MCSD = round(MCSD, 3),
                      spec = round(spec, 2),
                      focal_mean = round(focal_mean, 2),
                      b_shrunk = round(b_shrunk, 3))
      pk <- paste0(fc, "->", nc)
      pair_pp[[pk]] <- tibble::tibble(pair = pk, pair_pct = pp,
                                       n_sig = length(sig_idx),
                                       top_n_share = if (nrow(df)) round(sum(df$MCSD), 3) else NA_real_)
      pair_genes[[pk]] <- df |> dplyr::mutate(pair = pk) |>
        dplyr::select(pair, gene, sign, b_shrunk, R2_pct, spec, focal_mean, MCSD, lfsr)
    }
  }
  list(pair_summary = dplyr::bind_rows(pair_pp) |> dplyr::arrange(dplyr::desc(pair_pct)),
       pair_genes   = pair_genes)
}

## ---- Print canonical top-3 pair × top-3 gene tables for each cohort ----
print_canonical_tables <- function(rds_paths = default_rds, top_n_pairs = 3L, top_n_genes = 3L) {
  for (lbl in names(rds_paths)) {
    info <- rds_paths[[lbl]]
    mv <- readRDS(info$rds)
    res <- compute_pair_top_drivers(mv, info$resp_term, top_n_genes)
    cat(sprintf("\n========= %s — top-%d pairs × top-%d genes (locked canonical) =========\n",
                lbl, top_n_pairs, top_n_genes))
    for (pk in head(res$pair_summary$pair, top_n_pairs)) {
      r <- res$pair_summary[res$pair_summary$pair == pk, ]
      cat(sprintf("\n%s  (pair_pct=%.3f, n_sig=%d, top%d_share=%.3f)\n",
                  pk, r$pair_pct, r$n_sig, top_n_genes, r$top_n_share))
      print(res$pair_genes[[pk]] |> dplyr::select(-pair), n = top_n_genes, width = 200)
    }
  }
}

## When sourced from a different working directory, expose helpers.
## When run directly (Rscript), produce the canonical tables.
if (!interactive() && identical(Sys.getenv("R_TOP_DRIVERS_SOURCE_ONLY"), "")) {
  print_canonical_tables()
}

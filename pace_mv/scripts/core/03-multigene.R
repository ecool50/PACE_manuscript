# ==============================================================================
# 03-multigene.R
# Multi-gene aggregation for variance decomposition
# Works natively with fit_light objects from 00-extract.R
# ==============================================================================

#' Multi-gene Block Summary
#'
#' Unified interface for both responder-aware and no-responder analysis.
#' When \code{resp_term} is NULL, no responder logic is applied: the per_gene
#' tibble includes \code{cellstate_SS}, the summary uses \code{ratio_to_celltype},
#' and the return list includes \code{ratio_cellstate_celltype}.
#' When \code{resp_term} is not NULL, the summary uses \code{ratio_vs_residual}
#' and \code{ratio_vs_celltype} instead.
#'
#' @param fits_light Named list of fit_light objects (one per gene)
#' @param df Original data frame
#' @param residual Residual type
#' @param resp_term Responder term name (NULL for no-responder analysis)
#' @param group_ct Cell type grouping column
#' @param focal Optional focal cell type
#' @param neighbour Optional neighbour cell type(s) to include
#' @param gene_weights Optional named vector of gene weights
#' @param weight_by_resid_var Weight genes by inverse residual variance
#' @param BPPARAM BiocParallel backend
#' @return List with summary, weights, per_gene results, and (when resp_term is
#'   NULL) ratio_cellstate_celltype
block_sum_multigene_from_fits <- function(
        fits_light, df,
        residual = "working",
        resp_term = NULL,
        group_ct = "celltype",
        focal = NULL, neighbour = NULL,
        gene_weights = NULL, weight_by_resid_var = FALSE,
        BPPARAM = BiocParallel::SerialParam()
) {
    stopifnot(is.list(fits_light), length(fits_light) > 0)
    genes <- names(fits_light)
    stopifnot(!is.null(genes), all(nzchar(genes)))
    BPPARAM <- .bp_or_serial(BPPARAM)

    has_resp <- !is.null(resp_term)

    per_gene <- BiocParallel::bplapply(genes, function(g) {
        fl <- fits_light[[g]]
        out <- coef_ss_link(fl, df, residual = residual, focal = focal,
                            group_terms = c("celltype", "imageID"),
                            resp_term = resp_term)
        blk <- block_sum_coef_ss(out, resp_term = resp_term,
                                 group_ct = group_ct, neighbour = neighbour)
        row <- tibble(
            gene = g, block = blk$summary$block, SS_link = blk$summary$SS_link,
            n_terms = blk$summary$n_terms, n_obs = out$n_obs[1],
            resid_SS_total = blk$totals$residual_SS_total,
            celltype_SS = blk$totals$celltype_SS
        )
        if (!has_resp) row$cellstate_SS <- blk$totals$cellstate_SS
        row
    }, BPPARAM = BPPARAM) |> dplyr::bind_rows()

    # Weights
    w_df <- if (!is.null(gene_weights)) {
        tibble(gene = names(gene_weights), w = as.numeric(gene_weights))
    } else if (isTRUE(weight_by_resid_var)) {
        rv <- per_gene |>
            filter(block == "Residuals") |>
            distinct(gene, resid_SS_total, n_obs) |>
            mutate(var = resid_SS_total / pmax(1, n_obs), w = 1 / pmax(var, .Machine$double.eps))
        rv[, c("gene", "w")]
    } else {
        tibble(gene = unique(per_gene$gene), w = 1)
    }
    per_gene <- left_join(per_gene, w_df, by = "gene")

    summary <- per_gene |>
        group_by(block) |>
        summarise(SS_link = sum(w * SS_link, na.rm = TRUE),
                  n_terms = sum(n_terms, na.rm = TRUE), .groups = "drop")
    total_SS    <- sum(summary$SS_link, na.rm = TRUE)
    celltype_SS <- summary$SS_link[match("Cell type", summary$block)]

    if (has_resp) {
        resid_SS <- summary$SS_link[match("Residuals", summary$block)]
        summary <- summary |>
            mutate(
                pct_total = SS_link / total_SS,
                ratio_vs_residual = if (!is.na(resid_SS) && resid_SS != 0) SS_link / resid_SS else NA_real_,
                ratio_vs_celltype = if (!is.na(celltype_SS) && celltype_SS != 0) SS_link / celltype_SS else NA_real_
            ) |>
            arrange(desc(SS_link))
        list(summary = summary, weights = w_df, per_gene = per_gene)
    } else {
        summary <- summary |>
            mutate(
                pct_total = SS_link / total_SS,
                ratio_to_celltype = if (!is.na(celltype_SS) && celltype_SS != 0) SS_link / celltype_SS else NA_real_
            ) |>
            arrange(desc(SS_link))
        list(
            summary = summary, weights = w_df, per_gene = per_gene,
            ratio_cellstate_celltype = {
                cs <- summary$SS_link[match("Spatial cell state", summary$block)]
                if (!is.na(celltype_SS) && celltype_SS != 0) cs / celltype_SS else NA_real_
            }
        )
    }
}

#' Multi-gene Block Summary (No Responder) - backward-compatible wrapper
#'
#' @inheritParams block_sum_multigene_from_fits
#' @return List with summary, weights, per_gene, and ratio_cellstate_celltype
block_sum_multigene_from_fits_no_resp <- function(
        fits_light, df, residual = "working", group_ct = "celltype",
        focal = NULL, neighbour = NULL, gene_weights = NULL,
        weight_by_resid_var = FALSE, BPPARAM = BiocParallel::SerialParam()
) {
    block_sum_multigene_from_fits(fits_light, df, residual = residual,
                                  resp_term = NULL, group_ct = group_ct,
                                  focal = focal, neighbour = neighbour,
                                  gene_weights = gene_weights,
                                  weight_by_resid_var = weight_by_resid_var,
                                  BPPARAM = BPPARAM)
}


#' Run Block Decomposition for All Focal Cell Types
#'
#' Loops over all cell types as focal, aggregating multi-gene variance
#' decomposition for each. When \code{resp_term} is NULL, uses no-responder
#' logic; otherwise uses responder-aware logic.
#'
#' @param fits_light Named list of fit_light objects (one per gene)
#' @param df Original data frame
#' @param residual Residual type
#' @param resp_term Responder term name (NULL for no-responder analysis)
#' @param group_ct Cell type grouping column
#' @param weight_by_resid_var Weight genes by inverse residual variance
#' @param BPPARAM BiocParallel backend (used across genes within each focal)
#' @return List with summary (all focals), by_focal (per-focal results)
run_block_decomp <- function(
        fits_light, df,
        residual = "working",
        resp_term = NULL,
        group_ct = "celltype",
        weight_by_resid_var = TRUE,
        BPPARAM = BiocParallel::SerialParam()
) {
    cell_types <- levels(df[[group_ct]])
    if (is.null(cell_types)) cell_types <- unique(df[[group_ct]])

    by_focal <- list()
    for (ct in cell_types) {
        message("  ", ct)
        by_focal[[ct]] <- tryCatch({
            r <- block_sum_multigene_from_fits(
                fits_light = fits_light, df = df,
                residual = residual, resp_term = resp_term,
                group_ct = group_ct,
                focal = ct, neighbour = NULL,
                weight_by_resid_var = weight_by_resid_var,
                BPPARAM = BPPARAM
            )
            n_cells <- sum(df[[group_ct]] == ct, na.rm = TRUE)
            r$focal <- ct
            r$n_cells <- n_cells
            r
        }, error = function(e) {
            message("  Skip ", ct, ": ", conditionMessage(e))
            NULL
        })
    }
    by_focal <- Filter(Negate(is.null), by_focal)

    all_summary <- lapply(names(by_focal), function(ct) {
        r <- by_focal[[ct]]
        s <- as.data.frame(r$summary, stringsAsFactors = FALSE)
        s$focal <- ct
        s$n_cells <- r$n_cells
        s
    }) |> dplyr::bind_rows()

    list(summary = all_summary, by_focal = by_focal)
}

#' Run Block Decomposition (No Responder) - backward-compatible wrapper
#'
#' @inheritParams run_block_decomp
#' @return List with summary (all focals), by_focal (per-focal results)
run_block_decomp_no_resp <- function(fits_light, df, residual = "working",
                                      group_ct = "celltype",
                                      weight_by_resid_var = TRUE,
                                      BPPARAM = BiocParallel::SerialParam()) {
    run_block_decomp(fits_light, df, residual = residual, resp_term = NULL,
                     group_ct = group_ct, weight_by_resid_var = weight_by_resid_var,
                     BPPARAM = BPPARAM)
}

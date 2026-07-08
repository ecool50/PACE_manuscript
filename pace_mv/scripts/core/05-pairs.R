# ==============================================================================
# 05-pairs.R
# Focal-neighbour pair processing for variance decomposition
# Works natively with fit_light objects from 00-extract.R
# ==============================================================================

#' Process All Focal-Neighbour Pairs
#'
#' Unified interface for both responder-aware and no-responder analysis.
#' Sequential over pairs, parallel across genes inside.
#'
#' @param fits_light Named list of fit_light objects (one per gene)
#' @param df Original data frame
#' @param focals Character vector of focal cell types
#' @param neighbours Character vector of neighbour cell types
#' @param residual Residual type
#' @param resp_term Responder term name (NULL for no-responder analysis)
#' @param group_ct Cell type grouping column
#' @param weight_by_resid_var Weight genes by inverse residual variance
#' @param BPPARAM BiocParallel backend (used within pairs, not across)
#' @return Data frame of block summaries for all pairs
process_pairs <- function(
        fits_light, df, focals, neighbours,
        residual = "working",
        resp_term = NULL,
        group_ct = "celltype",
        weight_by_resid_var = TRUE,
        BPPARAM = BiocParallel::SerialParam()
) {
    pairs_df <- expand.grid(focal = focals, neighbour = neighbours, stringsAsFactors = FALSE)
    pairs_df <- pairs_df[pairs_df$focal != pairs_df$neighbour, , drop = FALSE]
    BPPARAM <- .bp_or_serial(BPPARAM)

    res_list <- vector("list", nrow(pairs_df))
    for (i in seq_len(nrow(pairs_df))) {
        foc <- pairs_df$focal[i]; nb <- pairs_df$neighbour[i]
        message(paste(foc, "~", nb))
        sm <- tryCatch({
            r <- block_sum_multigene_from_fits(
                fits_light = fits_light, df = df, residual = residual,
                resp_term = resp_term, group_ct = group_ct,
                focal = foc, neighbour = nb,
                weight_by_resid_var = weight_by_resid_var,
                BPPARAM = BPPARAM
            )
            df_sum <- as.data.frame(r$summary, stringsAsFactors = FALSE)
            df_sum$focal <- foc; df_sum$neighbour <- nb
            if (is.null(resp_term)) {
                df_sum$ratio_cellstate_celltype <- r$ratio_cellstate_celltype
            }
            df_sum
        }, error = function(e) {
            message("Skip ", foc, "~", nb, ": ", conditionMessage(e))
            NULL
        })
        res_list[[i]] <- sm
        if (i %% 10 == 0) gc(FALSE)
    }
    res_list <- Filter(Negate(is.null), res_list)
    if (length(res_list)) dplyr::bind_rows(res_list) else NULL
}

#' Process All Focal-Neighbour Pairs (No Responder) - backward-compatible wrapper
#'
#' @inheritParams process_pairs
#' @return Data frame of block summaries for all pairs
process_pairs_no_resp <- function(fits_light, df, focals, neighbours,
                                   residual = "working", group_ct = "celltype",
                                   weight_by_resid_var = TRUE,
                                   BPPARAM = BiocParallel::SerialParam()) {
    process_pairs(fits_light, df, focals, neighbours, residual = residual,
                  resp_term = NULL, group_ct = group_ct,
                  weight_by_resid_var = weight_by_resid_var, BPPARAM = BPPARAM)
}


# ==============================================================================
# FAST VERSION: compute per-coefficient SS once, reassign blocks per pair
# ==============================================================================

#' Process All Focal-Neighbour Pairs — Fast
#'
#' Unified interface for both responder-aware and no-responder analysis.
#' Precomputes per-gene, per-focal coefficient SS once, then reassigns blocks
#' for each neighbour. For K cell types this turns K*(K-1) full decompositions
#' into K full decompositions + K*(K-1) lightweight block reassignments.
#'
#' @param fits_light Named list of fit_light objects (one per gene)
#' @param df Original data frame
#' @param focals Character vector of focal cell types
#' @param neighbours Character vector of neighbour cell types
#' @param residual Residual type
#' @param resp_term Responder term name (NULL for no-responder analysis)
#' @param group_ct Cell type grouping column
#' @param weight_by_resid_var Weight genes by inverse residual variance
#' @param weight_by_detection Downweight sparse genes by detection rate
#' @param expr_mat Expression matrix (required when weight_by_detection is TRUE)
#' @param BPPARAM BiocParallel backend (used across genes, not pairs)
#' @return Data frame of block summaries for all pairs
process_pairs_fast <- function(
        fits_light, df, focals, neighbours,
        residual = "working",
        resp_term = NULL,
        group_ct = "celltype",
        weight_by_resid_var = TRUE,
        weight_by_detection = FALSE,
        expr_mat = NULL,
        BPPARAM = BiocParallel::SerialParam()
) {
    pairs_df <- expand.grid(focal = focals, neighbour = neighbours, stringsAsFactors = FALSE)
    pairs_df <- pairs_df[pairs_df$focal != pairs_df$neighbour, , drop = FALSE]
    BPPARAM <- .bp_or_serial(BPPARAM)

    genes <- names(fits_light)
    unique_focals <- unique(pairs_df$focal)

    # ================================================================
    # PHASE 1: For each focal, compute per-gene coefficient-level SS
    #          This is the expensive step -- done once per focal
    # ================================================================
    message("=== Phase 1: computing per-coefficient SS for each focal ===")

    coef_cache <- list()
    for (foc in unique_focals) {
        message("  Computing SS for focal = ", foc)
        coef_cache[[foc]] <- BiocParallel::bplapply(genes, function(g) {
            fl <- fits_light[[g]]
            if (is.null(fl)) return(NULL)
            tryCatch(
                coef_ss_link(fl, df, focal = foc,
                             group_terms = c("celltype", "imageID"),
                             resp_term = resp_term, resp_level = "PD",
                             residual = residual),
                error = function(e) NULL
            )
        }, BPPARAM = BPPARAM) |> stats::setNames(genes)
    }

    # ================================================================
    # PHASE 2: For each pair, reassign blocks with the target neighbour
    #          This is cheap -- just filtering and summing
    # ================================================================
    message("=== Phase 2: assembling ", nrow(pairs_df), " pairs ===")

    res_list <- vector("list", nrow(pairs_df))
    for (i in seq_len(nrow(pairs_df))) {
        foc <- pairs_df$focal[i]
        nb  <- pairs_df$neighbour[i]

        cache_foc <- coef_cache[[foc]]

        sm <- tryCatch({
            # Per-gene: apply block assignment with this neighbour
            per_gene <- lapply(genes, function(g) {
                out <- cache_foc[[g]]
                if (is.null(out) || nrow(out) == 0) return(NULL)
                blk <- block_sum_coef_ss(out, resp_term = resp_term,
                                         group_ct = group_ct, neighbour = nb)
                gene_row <- tibble::tibble(
                    gene = g, block = blk$summary$block, SS_link = blk$summary$SS_link,
                    n_terms = blk$summary$n_terms, n_obs = out$n_obs[1],
                    resid_SS_total = blk$totals$residual_SS_total,
                    celltype_SS = blk$totals$celltype_SS
                )
                if (is.null(resp_term)) {
                    gene_row$cellstate_SS <- blk$totals$cellstate_SS
                }
                gene_row
            })
            per_gene <- dplyr::bind_rows(Filter(Negate(is.null), per_gene))
            if (nrow(per_gene) == 0) stop("No genes produced results")

            # Weights
            w_df <- if (isTRUE(weight_by_resid_var)) {
                rv <- per_gene |>
                    dplyr::filter(block == "Residuals") |>
                    dplyr::distinct(gene, resid_SS_total, n_obs) |>
                    dplyr::mutate(
                        var = resid_SS_total / pmax(1, n_obs),
                        w = 1 / pmax(var, .Machine$double.eps)
                    )
                rv[, c("gene", "w")]
            } else {
                tibble::tibble(gene = unique(per_gene$gene), w = 1)
            }

            # Detection rate weighting: downweight sparse genes
            if (isTRUE(weight_by_detection) && !is.null(expr_mat)) {
                focal_idx <- df[[group_ct]] == foc
                det_df <- tibble::tibble(
                    gene = unique(per_gene$gene),
                    det_rate = vapply(unique(per_gene$gene), function(g) {
                        if (g %in% colnames(expr_mat)) {
                            mean(expr_mat[focal_idx, g] > 0, na.rm = TRUE)
                        } else { 1.0 }
                    }, numeric(1))
                )
                w_df <- dplyr::left_join(w_df, det_df, by = "gene") |>
                    dplyr::mutate(w = w * det_rate)
            }

            per_gene <- dplyr::left_join(per_gene, w_df, by = "gene")

            # Aggregate
            summary <- per_gene |>
                dplyr::group_by(block) |>
                dplyr::summarise(
                    SS_link = sum(w * SS_link, na.rm = TRUE),
                    n_terms = sum(n_terms, na.rm = TRUE),
                    .groups = "drop"
                )
            total_SS <- sum(summary$SS_link, na.rm = TRUE)
            celltype_SS <- summary$SS_link[match("Cell type", summary$block)]

            if (is.null(resp_term)) {
                # No-responder mode: ratio_to_celltype + ratio_cellstate_celltype
                cellstate_SS <- summary$SS_link[match("Spatial cell state", summary$block)]

                summary <- summary |>
                    dplyr::mutate(
                        pct_total = SS_link / total_SS,
                        ratio_to_celltype = if (!is.na(celltype_SS) && celltype_SS != 0) {
                            SS_link / celltype_SS
                        } else NA_real_
                    ) |>
                    dplyr::arrange(dplyr::desc(SS_link))

                ratio_cs_ct <- if (!is.na(celltype_SS) && celltype_SS != 0 && !is.na(cellstate_SS)) {
                    cellstate_SS / celltype_SS
                } else NA_real_

                df_sum <- as.data.frame(summary, stringsAsFactors = FALSE)
                df_sum$focal <- foc
                df_sum$neighbour <- nb
                df_sum$ratio_cellstate_celltype <- ratio_cs_ct
            } else {
                # Responder mode: ratio_vs_residual + ratio_vs_celltype
                resid_SS <- summary$SS_link[match("Residuals", summary$block)]

                summary <- summary |>
                    dplyr::mutate(
                        pct_total = SS_link / total_SS,
                        ratio_vs_residual = if (!is.na(resid_SS) && resid_SS != 0) SS_link / resid_SS else NA_real_,
                        ratio_vs_celltype = if (!is.na(celltype_SS) && celltype_SS != 0) SS_link / celltype_SS else NA_real_
                    ) |>
                    dplyr::arrange(dplyr::desc(SS_link))

                df_sum <- as.data.frame(summary, stringsAsFactors = FALSE)
                df_sum$focal <- foc
                df_sum$neighbour <- nb
            }
            df_sum
        }, error = function(e) {
            message("Skip ", foc, "~", nb, ": ", conditionMessage(e))
            NULL
        })

        res_list[[i]] <- sm
    }

    res_list <- Filter(Negate(is.null), res_list)
    if (length(res_list)) dplyr::bind_rows(res_list) else NULL
}

#' Process All Focal-Neighbour Pairs (No Responder) Fast - backward-compatible wrapper
#'
#' @inheritParams process_pairs_fast
#' @return Data frame of block summaries for all pairs
process_pairs_no_resp_fast <- function(fits_light, df, focals, neighbours,
                                        residual = "working", group_ct = "celltype",
                                        weight_by_resid_var = TRUE,
                                        weight_by_detection = FALSE, expr_mat = NULL,
                                        BPPARAM = BiocParallel::SerialParam()) {
    process_pairs_fast(fits_light, df, focals, neighbours, residual = residual,
                       resp_term = NULL, group_ct = group_ct,
                       weight_by_resid_var = weight_by_resid_var,
                       weight_by_detection = weight_by_detection,
                       expr_mat = expr_mat, BPPARAM = BPPARAM)
}

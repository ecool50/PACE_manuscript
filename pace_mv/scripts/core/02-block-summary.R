# ==============================================================================
# 02-block-summary.R
# Semantic block categorization for variance decomposition
# ==============================================================================

block_sum_coef_ss <- function(out, resp_term = NULL, group_ct = "celltype", neighbour = NULL) {
    stopifnot(all(c("type","group","term","coef_id","SS_link") %in% names(out)))
    term <- as.character(out$term); term[is.na(term)] <- ""
    block <- rep("Residuals", nrow(out))
    is_ct <- out$type == "ranef" & out$group == group_ct

    has_resp <- !is.null(resp_term)

    if (any(is_ct)) {
        is_spill <- grepl("spill|_near$", term, ignore.case = TRUE)
        idx_ct   <- which(is_ct)
        idx_int  <- idx_ct[term[idx_ct] == "(Intercept)"]

        block[idx_int]  <- "Cell type"
        block[is_spill] <- "Spillover"

        if (has_resp) {
            resp_pref <- paste0("^", resp_term, ":")
            idx_resp  <- idx_ct[term[idx_ct] == resp_term]
            idx_state <- idx_ct[term[idx_ct] != "(Intercept)" & !grepl(resp_term, term[idx_ct], fixed = TRUE)]
            idx_spat  <- idx_ct[grepl(resp_pref, term[idx_ct])]
            block[idx_resp]  <- "Responder status"
            block[idx_state] <- "Spatial cell state"
            block[idx_spat]  <- "Responder spatial state"
            if (!is.null(neighbour) && length(neighbour)) {
                neighbour <- unique(as.character(neighbour))
                ok_state   <- term %in% neighbour
                ok_spatial <- term %in% paste0(resp_term, ":", neighbour)
                block[block == "Spatial cell state" & !ok_state]           <- "Residuals"
                block[block == "Responder spatial state" & !ok_spatial]    <- "Residuals"
            }
        } else {
            idx_state <- idx_ct[term[idx_ct] != "(Intercept)"]
            block[idx_state] <- "Spatial cell state"
            if (!is.null(neighbour) && length(neighbour)) {
                neighbour <- unique(as.character(neighbour))
                block[block == "Spatial cell state" & !(term %in% neighbour)] <- "Residuals"
            }
        }
    }

    # fixed-effect spillover (both modes)
    block[out$type == "fixed" & grepl("spill|_near$", term, ignore.case = TRUE)] <- "Spillover"

    # fixed-effect Responder (responder mode only). The model now includes
    # `Responder` as a fixed effect to stabilise the random Responder BLUPs;
    # capture its variance contribution under "Responder status" rather than
    # letting it fall through to "Residuals".
    if (has_resp) {
        is_fixed_resp <- out$type == "fixed" &
                         grepl(paste0("^", resp_term), term)
        block[is_fixed_resp] <- "Responder status"
    }

    annotated <- mutate(out, block = block)
    summary <- annotated |>
        group_by(block) |>
        summarise(SS_link = sum(SS_link, na.rm = TRUE), n_terms = dplyr::n(), .groups = "drop")

    total_SS    <- sum(summary$SS_link, na.rm = TRUE)
    resid_SS    <- summary$SS_link[match("Residuals", summary$block)]
    celltype_SS <- summary$SS_link[match("Cell type", summary$block)]

    if (has_resp) {
        # ---- responder mode: ratio_vs_residual + ratio_vs_celltype ----
        ratio_vs_resid    <- if (!is.na(resid_SS) && resid_SS != 0) summary$SS_link / resid_SS else NA_real_
        ratio_vs_celltype <- if (!is.na(celltype_SS) && celltype_SS != 0) summary$SS_link / celltype_SS else NA_real_
        summary <- summary |>
            mutate(pct_total        = SS_link / total_SS,
                   ratio_vs_residual  = ratio_vs_resid,
                   ratio_vs_celltype  = ratio_vs_celltype) |>
            arrange(desc(SS_link))

        totals <- list(
            total_SS              = total_SS,
            residual_SS_total     = resid_SS,
            pure_model_residual_SS = sum(annotated$SS_link[annotated$type == "residual"], na.rm = TRUE),
            folded_other_SS       = resid_SS - sum(annotated$SS_link[annotated$type == "residual"], na.rm = TRUE),
            celltype_SS           = celltype_SS,
            residual_fraction     = resid_SS / total_SS,
            explained_fraction    = (total_SS - resid_SS) / total_SS
        )
    } else {
        # ---- no-resp mode: ratio_to_celltype ----
        summary <- summary |>
            mutate(pct_total        = SS_link / total_SS,
                   ratio_to_celltype = if (!is.na(celltype_SS) && celltype_SS != 0) SS_link / celltype_SS else NA_real_) |>
            arrange(desc(SS_link))

        totals <- list(
            total_SS              = total_SS,
            residual_SS_total     = resid_SS,
            pure_model_residual_SS = sum(annotated$SS_link[annotated$type == "residual"], na.rm = TRUE),
            folded_other_SS       = resid_SS - sum(annotated$SS_link[annotated$type == "residual"], na.rm = TRUE),
            celltype_SS           = celltype_SS,
            cellstate_SS          = summary$SS_link[match("Spatial cell state", summary$block)],
            spillover_SS          = summary$SS_link[match("Spillover", summary$block)],
            ratio_cellstate_celltype = if (!is.na(celltype_SS) && celltype_SS != 0)
                summary$SS_link[match("Spatial cell state", summary$block)] / celltype_SS else NA_real_
        )
    }

    list(summary = summary, annotated = annotated, totals = totals)
}

# Backward-compatible wrapper
block_sum_coef_ss_no_resp <- function(out, group_ct = "celltype", neighbour = NULL) {
    block_sum_coef_ss(out, resp_term = NULL, group_ct = group_ct, neighbour = neighbour)
}

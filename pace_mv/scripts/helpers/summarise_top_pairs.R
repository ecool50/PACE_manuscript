#' Summarise Top Spatial Cell State Pairs with MCSD Gene Drivers
#'
#' Takes the output of process_pairs_no_resp_fast, identifies the top N pairs
#' by spatial cell state contribution, then computes MCSD for each pair to
#' identify the top gene drivers.
#'
#' @param pairs_df Data frame output from process_pairs_no_resp_fast
#' @param fits_light Named list of fit_light objects
#' @param df Original data frame
#' @param n_pairs Number of top pairs to summarise (default: 10)
#' @param n_genes Number of top MCSD genes per pair (default: 2)
#' @param spill_correct Spillover correction mode (default: "weight")
#' @param spill_source_mode Source mode for spillover (default: "any")
#' @param detection_weight Detection weighting mode (default: "post")
#' @param min_expr_prev Minimum expression prevalence (default: 0.10)
#' @param min_expr_var Minimum expression variance (default: 1.0)
#' @param BPPARAM BiocParallel backend
#' @return A tibble with columns: rank, focal, neighbour, pct_total, gene, MCSD
#' @export
summarise_top_pairs <- function(
        pairs_df, fits_light, df,
        n_pairs = 10, n_genes = 2,
        spill_correct = "weight",
        spill_source_mode = "any",
        detection_weight = "post",
        min_expr_prev = 0.10,
        min_expr_var = 1.0,
        BPPARAM = BiocParallel::SerialParam()
) {
    # Get top pairs by spatial cell state variance (SS_link)
    top_pairs <- pairs_df |>
        dplyr::filter(block == "Spatial cell state") |>
        dplyr::arrange(dplyr::desc(pct_total)) |>
        dplyr::slice_head(n = n_pairs) |>
        dplyr::select(focal, neighbour, SS_link, pct_total)
    
    # Compute MCSD for each pair
    results_list <- vector("list", nrow(top_pairs))
    for (i in seq_len(nrow(top_pairs))) {
        foc <- gsub(" ", "_", top_pairs$focal[i])
        nb  <- gsub(" ", "_", top_pairs$neighbour[i])
        message("Computing MCSD for ", foc, " ~ ", nb, " (", i, "/", nrow(top_pairs), ")")
        
        mcsd_result <- tryCatch(
            mcsd_block_from_fits_no_resp(
                fits_light, df,
                block = "Spatial cell state",
                focal = foc, neighbour = nb,
                center = TRUE,
                min_expr_prev = min_expr_prev,
                min_expr_var = min_expr_var,
                spill_correct = spill_correct,
                spill_source_mode = spill_source_mode,
                detection_weight = detection_weight,
                BPPARAM = BPPARAM
            ),
            error = function(e) {
                message("  Failed: ", conditionMessage(e))
                NULL
            }
        )
        
        if (!is.null(mcsd_result)) {
            top_genes <- mcsd_result$scores |>
                dplyr::slice_head(n = n_genes) |>
                dplyr::mutate(
                    focal = foc,
                    neighbour = nb,
                    SS_link = top_pairs$SS_link[i],
                    pct_total = top_pairs$pct_total[i],
                    pair_rank = i
                )
            results_list[[i]] <- top_genes
        }
    }
    
    out <- dplyr::bind_rows(Filter(Negate(is.null), results_list)) |>
        dplyr::select(pair_rank, focal, neighbour, SS_link, pct_total, gene, MCSD) |>
        dplyr::arrange(pair_rank)
    
    out
}


#' Summarise Top Pairs with MCSD Gene Drivers (Responder-aware)
#'
#' Takes the output of process_pairs_fast, identifies the top N pairs
#' by a specified block contribution, then computes MCSD for each pair to
#' identify the top gene drivers.
#'
#' @param pairs_df Data frame output from process_pairs_fast
#' @param fits_light Named list of fit_light objects
#' @param df Original data frame
#' @param block Block to rank pairs by and compute MCSD for
#'   (default: "Responder spatial state")
#' @param n_pairs Number of top pairs to summarise (default: 10)
#' @param n_genes Number of top MCSD genes per pair (default: 2)
#' @param resp_term Responder term name (default: "ResponderPD")
#' @param resp_level Responder level (default: "PD")
#' @param spill_correct Spillover correction mode (default: "weight")
#' @param min_expr_prev Minimum expression prevalence (default: 0.10)
#' @param min_expr_var Minimum expression variance (default: 1.0)
#' @param BPPARAM BiocParallel backend
#' @return A tibble with columns: pair_rank, focal, neighbour, pct_total, gene, MCSD
#' @export
summarise_top_pairs_resp <- function(
        pairs_df, fits_light, df,
        block = "Responder spatial state",
        n_pairs = 10, n_genes = 2,
        resp_term = "ResponderPD",
        resp_level = "PD",
        spill_correct = "weight",
        spill_source_mode = "any",
        detection_weight = "post",
        min_expr_prev = 0.10,
        min_expr_var = 1.0,
        BPPARAM = BiocParallel::SerialParam()
) {
    # Get top pairs by the specified block
    top_pairs <- pairs_df |>
        dplyr::filter(block == !!block) |>
        dplyr::arrange(dplyr::desc(pct_total)) |>
        dplyr::slice_head(n = n_pairs) |>
        dplyr::select(focal, neighbour, SS_link, pct_total)
    
    if (nrow(top_pairs) == 0) stop("No pairs found for block = '", block, "'")
    
    # Compute MCSD for each pair
    results_list <- vector("list", nrow(top_pairs))
    for (i in seq_len(nrow(top_pairs))) {
        foc <- gsub(" ", "_", top_pairs$focal[i])
        nb  <- gsub(" ", "_", top_pairs$neighbour[i])
        message("Computing MCSD for ", foc, " ~ ", nb, " (", i, "/", nrow(top_pairs), ")")
        
        mcsd_result <- tryCatch(
            mcsd_block_from_fits(
                fits_light, df,
                block = block,
                focal = foc, neighbour = nb,
                group_ct = "celltype",
                resp_term = resp_term,
                resp_level = resp_level,
                center = TRUE,
                min_expr_prev = min_expr_prev,
                min_expr_var = min_expr_var,
                spill_correct = spill_correct,
                spill_source_mode = spill_source_mode,
                detection_weight = detection_weight,
                BPPARAM = BPPARAM
            ),
            error = function(e) {
                message("  Failed: ", conditionMessage(e))
                NULL
            }
        )
        
        if (!is.null(mcsd_result)) {
            top_genes <- mcsd_result$scores |>
                dplyr::slice_head(n = n_genes) |>
                dplyr::mutate(
                    focal = foc,
                    neighbour = nb,
                    SS_link = top_pairs$SS_link[i],
                    pct_total = top_pairs$pct_total[i],
                    pair_rank = i
                )
            results_list[[i]] <- top_genes
        }
    }
    
    out <- dplyr::bind_rows(Filter(Negate(is.null), results_list))
    if (nrow(out) == 0) {
        warning("All pairs failed MCSD computation.")
        return(tibble::tibble())
    }
    
    out |>
        dplyr::select(pair_rank, focal, neighbour, SS_link, pct_total, gene, MCSD) |>
        dplyr::arrange(pair_rank)
}
#'
#' @param summary_df Output from summarise_top_pairs
#' @return Invisible; prints formatted output
#' @export
print_pair_summary <- function(summary_df) {
    pairs <- unique(summary_df$pair_rank)
    for (p in pairs) {
        rows <- summary_df[summary_df$pair_rank == p, ]
        cat(sprintf("\n#%d  %s → %s  (SS = %.1f, %.4f%% total)\n",
                    p, rows$focal[1], rows$neighbour[1], rows$SS_link[1], rows$pct_total[1] * 100))
        for (j in seq_len(nrow(rows))) {
            cat(sprintf("     %s (MCSD = %.3f)\n", rows$gene[j], rows$MCSD[j]))
        }
    }
    invisible(summary_df)
}

#' Plot top MCSD pairs with gene annotations
#'
#' Shows all genes as points but only labels the top n_label per pair.
#'
#' @param summary_df Output of summarise_top_pairs()
#' @param n_pairs Number of pairs to show (default: NULL = all)
#' @param n_label Number of genes to label per pair (default: NULL = label all)
#' @param title Plot title
#' @return ggplot object
#' @export
plot_top_pairs <- function(summary_df, n_pairs = NULL, n_label = NULL,
                           title = "Top spatial cell state pairs — MCSD gene scores",
                           pct_format = function(x) sprintf("%.2f%%", x * 100)) {
    
    library(ggplot2)
    library(ggrepel)
    
    if (!is.null(n_pairs)) {
        keep_ranks <- sort(unique(summary_df$pair_rank))[seq_len(min(n_pairs,
                                                                     length(unique(summary_df$pair_rank))))]
        summary_df <- summary_df[summary_df$pair_rank %in% keep_ranks, ]
    }
    
    summary_df$pair_label <- paste0(gsub("_", " ", summary_df$focal),
                                    " \u2192 ", gsub("_", " ", summary_df$neighbour))
    
    pair_order <- dplyr::arrange(dplyr::distinct(summary_df,
                                                 pair_label, pct_total), pct_total)
    summary_df$pair_label <- factor(summary_df$pair_label, levels = pair_order$pair_label)
    
    summary_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(summary_df,
                                                               pair_rank), gene_rank = dplyr::row_number()))
    
    # Determine which genes get text labels
    if (!is.null(n_label)) {
        summary_df$show_label <- summary_df$gene_rank <= n_label
    } else {
        summary_df$show_label <- TRUE
    }
    
    n_genes <- max(summary_df$gene_rank)
    if (n_genes == 1) {
        fill_vals <- c(`1` = "#B2182B")
        fill_labels <- c(`1` = "Top gene")
    } else if (n_genes == 2) {
        fill_vals <- c(`1` = "#B2182B", `2` = "#2166AC")
        fill_labels <- c(`1` = "Top gene", `2` = "2nd gene")
    } else {
        # Distinct palette for labelled genes, grey for the rest
        top_colours <- c("#B2182B", "#E66101", "#5E3C99", "#1B7837",
                         "#2166AC", "#D95F02", "#7570B3", "#E7298A",
                         "#66A61E", "#E6AB02")
        all_vals <- character(n_genes)
        all_labels <- character(n_genes)
        for (i in seq_len(n_genes)) {
            key <- as.character(i)
            if (i <= length(top_colours)) {
                all_vals[i] <- top_colours[i]
            } else {
                all_vals[i] <- "grey60"
            }
            if (i == 1) {
                all_labels[i] <- "Top gene"
            } else {
                all_labels[i] <- paste0(i, ordinal_suffix(i), " gene")
            }
        }
        fill_vals <- setNames(all_vals, as.character(seq_len(n_genes)))
        fill_labels <- setNames(all_labels, as.character(seq_len(n_genes)))
    }
    
    # Limit legend to only the labelled gene ranks
    if (!is.null(n_label)) {
        legend_ranks <- seq_len(n_label)
        fill_vals <- fill_vals[as.character(legend_ranks)]
        fill_labels <- fill_labels[as.character(legend_ranks)]
    }
    
    pct_labels <- dplyr::mutate(pair_order,
                                label = paste0(pair_label, "  (", pct_format(pct_total), ")"))
    label_map <- setNames(pct_labels$label, pct_labels$pair_label)
    
    # Split into background (unlabelled) and foreground (labelled) layers
    bg_data <- summary_df[!summary_df$show_label, ]
    fg_data <- summary_df[summary_df$show_label, ]

    p <- ggplot(summary_df, aes(x = MCSD, y = pair_label)) +
        geom_point(data = bg_data,
                   aes(size = MCSD, fill = factor(gene_rank)),
                   shape = 21, colour = "grey30", stroke = 0.4,
                   show.legend = FALSE) +
        geom_point(data = fg_data,
                   aes(size = MCSD, fill = factor(gene_rank)),
                   shape = 21, colour = "grey30", stroke = 0.4,
                   show.legend = TRUE) +
        ggrepel::geom_text_repel(
            data = summary_df[summary_df$show_label, ],
            aes(label = gene),
            size = 3, fontface = "italic",
            nudge_y = 0.15,
            segment.size = 0.25,
            segment.color = "grey60",
            min.segment.length = 0,
            max.overlaps = 20,
            box.padding = 0.3
        ) +
        scale_y_discrete(labels = label_map) +
        scale_size_continuous(range = c(2, 7), guide = "none") +
        scale_fill_manual(values = fill_vals, labels = fill_labels,
                          breaks = names(fill_vals), name = NULL) +
        labs(x = "MCSD (Spatial cell state)", y = NULL, title = title) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 9),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            plot.margin = margin(10, 15, 10, 10)
        )
    
    p
}

# Helper for ordinal suffixes
ordinal_suffix <- function(x) {
    ifelse(x %% 100 %in% c(11, 12, 13), "th",
           ifelse(x %% 10 == 1, "st",
                  ifelse(x %% 10 == 2, "nd",
                         ifelse(x %% 10 == 3, "rd", "th"))))
}

#' Pooled Contamination Diagnostic Scatter
#'
#' Two-panel scatter (uncorrected | corrected) pooling all genes across
#' top N pairs. Shows expression ratio vs spatial slope with regression
#' line and correlation annotation.
#'
#' @param results_corrected Corrected model results (list of per-gene results)
#' @param results_uncorrected Uncorrected model results
#' @param df Data frame with celltype and gene expression columns
#' @param expr_mat Expression matrix (cells x genes)
#' @param pairs_df Output from process_pairs_no_resp_fast
#' @param n_pairs Number of top pairs to include (default: 10)
#' @param group_ct Cell type column name
#' @param title Plot title
#' @return ggplot object
#' @export
plot_pooled_contamination <- function(
        results_corrected, results_uncorrected,
        df, expr_mat,
        pairs_df,
        n_pairs = 10,
        group_ct = "celltype",
        title = "Correction eliminates contamination confounding",
        panels = c("Uncorrected", "Corrected")
) {
    library(ggplot2)
    
    # Get top pairs
    top_pairs <- pairs_df |>
        dplyr::filter(block == "Spatial cell state") |>
        dplyr::arrange(dplyr::desc(pct_total)) |>
        dplyr::slice_head(n = n_pairs)
    top_pairs$focal <- gsub(" ", "_", top_pairs$focal)
    top_pairs$neighbour <- gsub(" ", "_", top_pairs$neighbour)
    
    genes <- intersect(names(results_corrected), names(results_uncorrected))
    genes <- intersect(genes, colnames(expr_mat))
    
    # Collect slopes across all pairs
    all_slopes <- vector("list", nrow(top_pairs))
    
    for (i in seq_len(nrow(top_pairs))) {
        foc <- top_pairs$focal[i]
        nb  <- top_pairs$neighbour[i]
        focal_idx <- df[[group_ct]] == foc
        neighbour_idx <- df[[group_ct]] == nb
        
        pair_slopes <- lapply(genes, function(g) {
            # Corrected
            res_c <- results_corrected[[g]]
            if (is.null(res_c) || isTRUE(res_c$error)) return(NULL)
            rv_c <- res_c$ran_vals
            if (is.null(rv_c)) return(NULL)
            row_c <- rv_c$level == foc & rv_c$term == nb
            if (!any(row_c)) return(NULL)
            slope_corr <- rv_c$estimate[row_c][1]
            
            # Uncorrected
            res_u <- results_uncorrected[[g]]
            if (is.null(res_u) || isTRUE(res_u$error)) return(NULL)
            rv_u <- res_u$ran_vals
            if (is.null(rv_u)) return(NULL)
            row_u <- rv_u$level == foc & rv_u$term == nb
            if (!any(row_u)) return(NULL)
            slope_uncorr <- rv_u$estimate[row_u][1]
            
            # Expression ratio
            focal_expr <- mean(expr_mat[focal_idx, g], na.rm = TRUE)
            neighbour_expr <- mean(expr_mat[neighbour_idx, g], na.rm = TRUE)
            expr_ratio <- neighbour_expr / (focal_expr + 0.01)
            
            data.frame(
                gene = g, focal = foc, neighbour = nb,
                slope_corrected = slope_corr,
                slope_uncorrected = slope_uncorr,
                expr_ratio = expr_ratio
            )
        }) |> purrr::compact() |> dplyr::bind_rows()
        
        all_slopes[[i]] <- pair_slopes
    }
    
    pooled <- dplyr::bind_rows(all_slopes)
    if (nrow(pooled) == 0) stop("No slopes extracted.")
    
    # Compute correlations
    r_uncorr <- cor(pooled$expr_ratio, pooled$slope_uncorrected, use = "complete.obs")
    r_corr   <- cor(pooled$expr_ratio, pooled$slope_corrected, use = "complete.obs")
    
    # Long format
    df_long <- pooled |>
        tidyr::pivot_longer(
            cols = c(slope_uncorrected, slope_corrected),
            names_to = "model",
            values_to = "slope"
        ) |>
        dplyr::mutate(
            model = dplyr::case_when(
                model == "slope_uncorrected" ~ "Uncorrected",
                model == "slope_corrected" ~ "Corrected"
            ),
            model = factor(model, levels = c("Uncorrected", "Corrected"))
        ) |>
        dplyr::filter(model %in% panels)

    # Annotation labels
    r_vals <- c(Uncorrected = r_uncorr, Corrected = r_corr)
    r_labels <- data.frame(
        model = factor(intersect(c("Uncorrected", "Corrected"), panels),
                       levels = c("Uncorrected", "Corrected")),
        label = paste0("italic(r) == ", sprintf("%.2f", r_vals[intersect(c("Uncorrected", "Corrected"), panels)]))
    )
    
    p <- ggplot(df_long, aes(x = expr_ratio, y = slope)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
        geom_point(alpha = 0.15, size = 0.8, colour = "grey30") +
        geom_smooth(method = "lm", se = TRUE, colour = "#2166AC", 
                    fill = "#2166AC", alpha = 0.15, linewidth = 0.7) +
        geom_text(
            data = r_labels,
            aes(label = label),
            x = Inf, y = Inf,
            hjust = 1.1, vjust = 1.5,
            size = 4, parse = TRUE, colour = "#2166AC",
            inherit.aes = FALSE
        ) +
        facet_wrap(~ model) +
        labs(
            x = "Expression ratio (neighbour / focal)",
            y = "Spatial slope (random effect)",
            title = title,
            subtitle = paste0("Pooled across all", " focal-neighbour pairs")
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 10, colour = "grey40", hjust = 0.5),
            strip.text = element_text(face = "bold", size = 11),
            panel.grid.minor = element_blank(),
            plot.margin = margin(10, 15, 10, 10)
        )
    
    p
}

#' Pooled Correction Effect and Top Spillover Candidates
#'
#' Two-panel figure:
#' Left: Change in spatial slope vs expression ratio (pooled across pairs)
#' Right: Top N genes by max expression ratio, showing uncorrected vs corrected slopes
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results
#' @param df Data frame
#' @param expr_mat Expression matrix (cells x genes)
#' @param pairs_df Output from process_pairs_no_resp_fast
#' @param n_pairs Number of pairs to pool (default: all)
#' @param n_genes Number of top spillover genes to show (default: 15)
#' @param group_ct Cell type column name
#' @return ggplot object (patchwork)
#' @export
plot_correction_effect <- function(
        results_corrected, results_uncorrected,
        df, expr_mat,
        pairs_df,
        n_pairs = NULL,
        n_genes = 15,
        group_ct = "celltype"
) {
    library(ggplot2)
    library(patchwork)
    
    # Get pairs
    all_pairs <- pairs_df |>
        dplyr::filter(block == "Spatial cell state") |>
        dplyr::arrange(dplyr::desc(pct_total))
    if (!is.null(n_pairs)) all_pairs <- dplyr::slice_head(all_pairs, n = n_pairs)
    all_pairs$focal <- gsub(" ", "_", all_pairs$focal)
    all_pairs$neighbour <- gsub(" ", "_", all_pairs$neighbour)
    
    genes <- intersect(names(results_corrected), names(results_uncorrected))
    genes <- intersect(genes, colnames(expr_mat))
    
    # Collect slopes
    all_slopes <- vector("list", nrow(all_pairs))
    for (i in seq_len(nrow(all_pairs))) {
        foc <- all_pairs$focal[i]
        nb  <- all_pairs$neighbour[i]
        focal_idx <- df[[group_ct]] == foc
        neighbour_idx <- df[[group_ct]] == nb
        
        pair_slopes <- lapply(genes, function(g) {
            res_c <- results_corrected[[g]]
            res_u <- results_uncorrected[[g]]
            if (is.null(res_c) || isTRUE(res_c$error)) return(NULL)
            if (is.null(res_u) || isTRUE(res_u$error)) return(NULL)
            
            rv_c <- res_c$ran_vals; rv_u <- res_u$ran_vals
            if (is.null(rv_c) || is.null(rv_u)) return(NULL)
            
            row_c <- rv_c$level == foc & rv_c$term == nb
            row_u <- rv_u$level == foc & rv_u$term == nb
            if (!any(row_c) || !any(row_u)) return(NULL)
            
            focal_expr <- mean(expr_mat[focal_idx, g], na.rm = TRUE)
            neighbour_expr <- mean(expr_mat[neighbour_idx, g], na.rm = TRUE)
            
            data.frame(
                gene = g, focal = foc, neighbour = nb,
                pair = paste0(gsub("_", " ", foc), " \u2192 ", gsub("_", " ", nb)),
                slope_corrected = rv_c$estimate[row_c][1],
                slope_uncorrected = rv_u$estimate[row_u][1],
                expr_ratio = neighbour_expr / (focal_expr + 0.01)
            )
        }) |> purrr::compact() |> dplyr::bind_rows()
        all_slopes[[i]] <- pair_slopes
    }
    
    pooled <- dplyr::bind_rows(all_slopes)
    pooled$delta_slope <- pooled$slope_corrected - pooled$slope_uncorrected
    pooled$is_spillover <- pooled$expr_ratio > 1.5
    
    # ---- Panel 1: Correction effect scatter ----
    p1 <- ggplot(pooled, aes(x = expr_ratio, y = delta_slope)) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
        geom_point(aes(colour = is_spillover), alpha = 0.25, size = 0.8) +
        geom_smooth(method = "lm", se = TRUE, colour = "#2166AC",
                    fill = "#2166AC", alpha = 0.15, linewidth = 0.7) +
        scale_colour_manual(
            values = c("FALSE" = "grey50", "TRUE" = "#D6604D"),
            labels = c("FALSE" = "Ratio \u2264 1.5", "TRUE" = "Ratio > 1.5"),
            name = "Spillover\ncandidate"
        ) +
        labs(
            x = "Expression ratio (neighbour / focal)",
            y = "Change in spatial slope\n(corrected \u2212 uncorrected)",
            subtitle = "Correction effect by expression ratio"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.subtitle = element_text(face = "bold", size = 11),
            panel.grid.minor = element_blank(),
            legend.position = c(0.8, 0.85),
            legend.background = element_rect(fill = "white", colour = NA),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.4, "cm")
        )
    
    # ---- Panel 2: Top spillover candidates (arrow plot) ----
    # For each gene, take the max expression ratio across all pairs
    gene_max <- pooled |>
        dplyr::group_by(gene) |>
        dplyr::summarise(
            max_ratio = max(expr_ratio, na.rm = TRUE),
            slope_uncorrected = slope_uncorrected[which.max(expr_ratio)],
            slope_corrected = slope_corrected[which.max(expr_ratio)],
            pair = pair[which.max(expr_ratio)],
            .groups = "drop"
        ) |>
        dplyr::filter(max_ratio > 1.5) |>
        dplyr::arrange(dplyr::desc(max_ratio)) |>
        dplyr::slice_head(n = n_genes)
    
    # Order genes by expression ratio (highest at top)
    gene_max <- gene_max |>
        dplyr::mutate(
            gene_label = paste0(gene, "  (", sprintf("%.0f", max_ratio), "x)"),
            gene_label = factor(gene_label, levels = rev(gene_label))
        )
    
    p2 <- ggplot(gene_max, aes(y = gene_label)) +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
        # Arrow from uncorrected to corrected
        geom_segment(
            aes(x = slope_uncorrected, xend = slope_corrected, yend = gene_label),
            arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
            colour = "grey40", linewidth = 0.45
        ) +
        # Uncorrected point (start)
        geom_point(aes(x = slope_uncorrected), 
                   shape = 21, size = 3, fill = "#D6604D",
                   colour = "grey30", stroke = 0.4) +
        # Corrected point (end)
        geom_point(aes(x = slope_corrected), 
                   shape = 21, size = 3, fill = "#4393C3",
                   colour = "grey30", stroke = 0.4) +
        labs(
            x = "Spatial slope",
            y = NULL,
            subtitle = "Top spillover candidates (ordered by expression ratio)"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.subtitle = element_text(face = "bold", size = 11),
            axis.text.y = element_text(face = "italic", size = 9),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
            plot.margin = margin(10, 15, 10, 5)
        )
    
    p1 + p2 + plot_layout(widths = c(1, 1))
}
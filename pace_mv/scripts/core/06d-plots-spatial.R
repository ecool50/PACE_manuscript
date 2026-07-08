# ==============================================================================
# 06d-plots-spatial.R
# Spatial/Statial plots and evaluation functions
# ==============================================================================

#' Statial-Style Spillover ROC
#'
#' Ranks all gene-pair relationships by significance (p-value), then
#' plots cumulative true positives vs false positives.
#'
#' Cell type markers are defined empirically via one-vs-rest AUC. A gene
#' is a marker for cell type X if AUC > auc_threshold. If that gene shows
#' a significant spatial slope in focal cell type Y (Y != X) with neighbour X,
#' it is a false positive (contamination). All other significant relationships
#' are true positives (genuine biology).
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results
#' @param df Data frame
#' @param expr_mat Expression matrix (cells x genes)
#' @param auc_threshold AUC threshold for cell type marker classification (default: 0.70)
#' @param pval_cutoff P-value threshold; only relationships with pval < cutoff are
#'   included. Default NULL uses all relationships (original behaviour).
#' @param group_ct Cell type column name
#' @param top_n Zoom into top N relationships (default: 500)
#' @return List with plots, marker table, and data
#' @export
plot_statial_roc <- function(
        results_corrected, results_uncorrected,
        df, expr_mat,
        auc_threshold = 0.70,
        pval_cutoff = NULL,
        group_ct = "celltype",
        top_n = 500
) {
    library(ggplot2)
    cell_types <- sort(unique(df[[group_ct]]))
    genes_use <- intersect(names(results_corrected), names(results_uncorrected))
    genes_use <- intersect(genes_use, colnames(expr_mat))

    # ========== 1. Compute cell type markers via one-vs-rest AUC ==========
    message("Computing cell type markers (one-vs-rest AUC)...")

    compute_auc <- function(scores, labels) {
        n1 <- sum(labels)
        n0 <- sum(!labels)
        if (n1 == 0 || n0 == 0) return(0.5)
        r <- rank(scores)
        auc <- (sum(r[labels]) - n1 * (n1 + 1) / 2) / (n1 * n0)
        auc
    }

    marker_list <- dplyr::bind_rows(lapply(cell_types, function(ct) {
        is_ct <- df[[group_ct]] == ct
        aucs <- vapply(genes_use, function(g) {
            compute_auc(expr_mat[, g], is_ct)
        }, numeric(1))
        data.frame(gene = genes_use, cell_type = ct, auc = aucs,
                   is_marker = aucs >= auc_threshold, stringsAsFactors = FALSE)
    }))

    markers <- dplyr::filter(marker_list, is_marker)
    message(sprintf("  Found %d marker-cell type assignments (%d unique genes) at AUC >= %.2f",
                    nrow(markers), length(unique(markers$gene)), auc_threshold))

    # ========== 2. Extract all spatial slope relationships ==========
    extract_all_relationships <- function(results, label) {
        dplyr::bind_rows(purrr::compact(lapply(names(results), function(g) {
            res <- results[[g]]
            if (is.null(res) || isTRUE(res$error) || is.null(res$ran_vals)) return(NULL)
            rv <- res$ran_vals
            rv_slopes <- rv[rv$term != "(Intercept)", , drop = FALSE]
            if (nrow(rv_slopes) == 0) return(NULL)

            out <- data.frame(
                gene = rv_slopes$gene,
                focal = rv_slopes$level,
                neighbour = rv_slopes$term,
                estimate = rv_slopes$estimate,
                std.error = rv_slopes$std.error,
                model = label,
                stringsAsFactors = FALSE
            )
            if ("pval" %in% names(rv_slopes)) {
                out$pval <- rv_slopes$pval
            } else {
                z <- out$estimate / pmax(out$std.error, .Machine$double.eps)
                out$pval <- 2 * pnorm(-abs(z))
            }
            out
        })))
    }

    rel_corr   <- extract_all_relationships(results_corrected, "Corrected")
    rel_uncorr <- extract_all_relationships(results_uncorrected, "Uncorrected")

    # ========== 3. Label FP/TP ==========
    marker_lookup <- dplyr::transmute(markers, gene, neighbour = cell_type,
                                      is_neighbour_marker = TRUE)

    label_fp <- function(rel_df) {
        rel_df <- dplyr::filter(rel_df, gene %in% genes_use)
        rel_df <- dplyr::left_join(rel_df, marker_lookup, by = c("gene", "neighbour"))
        rel_df$is_neighbour_marker[is.na(rel_df$is_neighbour_marker)] <- FALSE
        rel_df$is_fp <- rel_df$is_neighbour_marker
        rel_df
    }

    message("Labelling relationships...")
    rel_corr   <- label_fp(rel_corr)
    rel_uncorr <- label_fp(rel_uncorr)

    # ========== 4. Apply p-value cutoff if specified ==========
    if (!is.null(pval_cutoff)) {
        n_before_corr   <- nrow(rel_corr)
        n_before_uncorr <- nrow(rel_uncorr)
        rel_corr   <- dplyr::filter(rel_corr,   !is.na(pval), pval < pval_cutoff)
        rel_uncorr <- dplyr::filter(rel_uncorr, !is.na(pval), pval < pval_cutoff)
        message(sprintf("  P-value cutoff %.2g: corrected %d -> %d, uncorrected %d -> %d",
                        pval_cutoff,
                        n_before_corr, nrow(rel_corr),
                        n_before_uncorr, nrow(rel_uncorr)))
    }

    n_fp_corr <- sum(rel_corr$is_fp)
    n_tp_corr <- sum(!rel_corr$is_fp)
    n_fp_uncorr <- sum(rel_uncorr$is_fp)
    n_tp_uncorr <- sum(!rel_uncorr$is_fp)
    message(sprintf("  Corrected:   %d TP, %d FP out of %d relationships",
                    n_tp_corr, n_fp_corr, nrow(rel_corr)))
    message(sprintf("  Uncorrected: %d TP, %d FP out of %d relationships",
                    n_tp_uncorr, n_fp_uncorr, nrow(rel_uncorr)))

    # ========== 5. Build cumulative curves ==========
    build_curve <- function(rel_df) {
        rel_df <- dplyr::arrange(dplyr::filter(rel_df, !is.na(pval), !is.na(is_fp)), pval)
        rel_df$cum_tp <- cumsum(!rel_df$is_fp)
        rel_df$cum_fp <- cumsum(rel_df$is_fp)
        rel_df
    }

    curve_corr   <- build_curve(rel_corr)
    curve_uncorr <- build_curve(rel_uncorr)

    # ========== 6. Plot ==========
    plot_data <- dplyr::mutate(
        dplyr::bind_rows(
            dplyr::select(curve_corr, cum_fp, cum_tp, model),
            dplyr::select(curve_uncorr, cum_fp, cum_tp, model)
        ),
        model = factor(model, levels = c("Uncorrected", "Corrected"))
    )

    base_theme <- theme_minimal(base_size = 14) +
        theme(
            plot.title    = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 11, colour = "#666666"),
            legend.position   = c(0.7, 0.2),
            legend.background = element_rect(fill = "white", colour = NA),
            legend.text       = element_text(size = 12),
            panel.grid.minor  = element_blank()
        )

    expected_slope <- n_tp_corr / max(n_fp_corr, 1)

    subtitle_text <- sprintf(
        "FP: neighbour's cell type marker (one-vs-rest AUC \u2265 %.2f) significant in focal cell type",
        auc_threshold
    )
    if (!is.null(pval_cutoff)) {
        subtitle_text <- paste0(subtitle_text, sprintf(" | p < %s", format(pval_cutoff, scientific = TRUE)))
    }

    p_full <- ggplot(plot_data, aes(x = cum_fp, y = cum_tp, colour = model)) +
        geom_abline(slope = expected_slope, intercept = 0,
                    linetype = "dashed", colour = "grey60", linewidth = 0.3) +
        geom_line(linewidth = 0.8) +
        scale_colour_manual(values = c(Uncorrected = "#4393C3", Corrected = "#D6604D"),
                            name = NULL) +
        labs(
            x = "Cumulative FP (neighbour cell type marker relationships)",
            y = "Cumulative TP (non-marker relationships)",
            title = "Spillover correction evaluation",
            subtitle = subtitle_text
        ) +
        base_theme

    p_zoom <- NULL
    if (!is.null(top_n)) {
        zoom_data <- dplyr::mutate(
            dplyr::bind_rows(
                dplyr::select(dplyr::slice_head(curve_corr, n = top_n), cum_fp, cum_tp, model),
                dplyr::select(dplyr::slice_head(curve_uncorr, n = top_n), cum_fp, cum_tp, model)
            ),
            model = factor(model, levels = c("Uncorrected", "Corrected"))
        )

        p_zoom <- ggplot(zoom_data, aes(x = cum_fp, y = cum_tp, colour = model)) +
            geom_line(linewidth = 0.8) +
            scale_colour_manual(values = c(Uncorrected = "#4393C3", Corrected = "#D6604D"),
                                name = NULL) +
            labs(
                x = "Cumulative FP (neighbour cell type marker relationships)",
                y = "Cumulative TP (non-marker relationships)",
                title = paste0("Top ", top_n, " most significant relationships"),
                subtitle = "Corrected model enriches for true spatial biology over contamination"
            ) +
            base_theme
    }

    list(
        plot_full         = p_full,
        plot_zoom         = p_zoom,
        markers           = markers,
        marker_summary    = marker_list,
        curve_corrected   = curve_corr,
        curve_uncorrected = curve_uncorr,
        n_relationships   = nrow(rel_corr),
        n_tp              = n_tp_corr,
        n_fp              = n_fp_corr,
        n_tp_uncorrected  = n_tp_uncorr,
        n_fp_uncorrected  = n_fp_uncorr
    )
}

#' FP Proportion Plot for Spillover Correction Evaluation
#'
#' Shows the false positive proportion (neighbour marker relationships)
#' as a function of the top K most significant relationships. Lower
#' FP proportion after correction = contamination removed.
#'
#' @param roc Output from plot_statial_roc
#' @param k_values Vector of K values to evaluate (default: seq)
#' @param k_max Maximum K to plot (default: 1000)
#' @param title Plot title
#' @return List with line plot, bar plot, and summary table
#' @export
plot_fp_proportion <- function(
        roc,
        k_max = 1000,
        k_step = 10,
        bar_k = c(100, 200, 500),
        title = "False positive rate in top spatial relationships"
) {
    library(ggplot2)
    library(patchwork)

    curve_corr <- roc$curve_corrected
    curve_uncorr <- roc$curve_uncorrected

    # ========== 1. Running FP proportion ==========
    k_values <- seq(k_step, min(k_max, nrow(curve_corr), nrow(curve_uncorr)), by = k_step)

    running <- lapply(k_values, function(k) {
        fp_corr <- sum(curve_corr$is_fp[1:k])
        fp_uncorr <- sum(curve_uncorr$is_fp[1:k])
        data.frame(
            k = k,
            fp_pct_corrected = 100 * fp_corr / k,
            fp_pct_uncorrected = 100 * fp_uncorr / k
        )
    }) |> dplyr::bind_rows()

    running_long <- running |>
        tidyr::pivot_longer(
            cols = c(fp_pct_corrected, fp_pct_uncorrected),
            names_to = "model",
            values_to = "fp_pct"
        ) |>
        dplyr::mutate(
            model = dplyr::case_when(
                model == "fp_pct_corrected" ~ "Corrected",
                model == "fp_pct_uncorrected" ~ "Uncorrected"
            ),
            model = factor(model, levels = c("Uncorrected", "Corrected"))
        )

    p_line <- ggplot(running_long, aes(x = k, y = fp_pct, colour = model)) +
        geom_line(linewidth = 0.8) +
        scale_colour_manual(
            values = c("Uncorrected" = "#4393C3", "Corrected" = "#D6604D"),
            name = NULL
        ) +
        labs(
            x = "Top K most significant relationships",
            y = "FP proportion (%)",
            title = title,
            subtitle = "FP = neighbour's cell type marker with significant spatial slope in focal"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 11, colour = "#666666"),
            legend.position = c(0.80, 0.85),
            legend.background = element_rect(fill = "white", colour = NA),
            legend.text = element_text(size = 12),
            panel.grid.minor = element_blank()
        )

    # ========== 2. Bar chart at specific K values ==========
    bar_data <- lapply(bar_k, function(k) {
        k_use <- min(k, nrow(curve_corr), nrow(curve_uncorr))
        fp_corr <- sum(curve_corr$is_fp[1:k_use])
        fp_uncorr <- sum(curve_uncorr$is_fp[1:k_use])
        data.frame(
            k = paste0("Top ", k),
            Corrected = 100 * fp_corr / k_use,
            Uncorrected = 100 * fp_uncorr / k_use
        )
    }) |> dplyr::bind_rows()

    bar_data$k <- factor(bar_data$k, levels = paste0("Top ", bar_k))

    bar_long <- bar_data |>
        tidyr::pivot_longer(
            cols = c(Corrected, Uncorrected),
            names_to = "model",
            values_to = "fp_pct"
        ) |>
        dplyr::mutate(model = factor(model, levels = c("Uncorrected", "Corrected")))

    p_bar <- ggplot(bar_long, aes(x = k, y = fp_pct, fill = model)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_text(
            aes(label = sprintf("%.1f%%", fp_pct)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4
        ) +
        scale_fill_manual(
            values = c("Uncorrected" = "#4393C3", "Corrected" = "#D6604D"),
            name = NULL
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(
            x = NULL,
            y = "FP proportion (%)",
            title = "False positive rate by significance tier"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(face = "bold", size = 16),
            legend.position = "top",
            legend.text = element_text(size = 12),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
        )

    list(
        plot_line = p_line,
        plot_bar = p_bar,
        combined = p_line + p_bar + plot_layout(widths = c(2, 1)),
        summary = bar_data,
        running = running
    )
}

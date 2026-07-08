# ==============================================================================
# 06c-plots-pairs.R
# Pair-level comparison and waterfall plots
# ==============================================================================

#' Spillover Volcano Plot
#'
#' For a given focal-neighbour pair, plots spatial slope vs expression ratio
#' for each gene. Highlights spillover candidates and shows correction effect.
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results (optional, shows arrows if provided)
#' @param df Data frame
#' @param expr_mat Expression matrix (cells x genes)
#' @param focal Focal cell type (underscore format)
#' @param neighbour Neighbour cell type (underscore format)
#' @param spill_threshold Expression ratio threshold for spillover (default: 1.5)
#' @param n_label Number of top genes to label (default: 10)
#' @param group_ct Cell type column name
#' @param source_mode "neighbour" or "any" for expression ratio computation
#' @return ggplot object
#' @export
plot_spillover_volcano <- function(
        results_corrected,
        results_uncorrected = NULL,
        df, expr_mat,
        focal, neighbour,
        spill_threshold = 1.5,
        n_label = 10,
        group_ct = "celltype",
        source_mode = c("neighbour", "any")
) {
    library(ggplot2)
    library(ggrepel)

    source_mode <- match.arg(source_mode)

    genes <- names(results_corrected)
    genes <- intersect(genes, colnames(expr_mat))

    focal_idx <- df[[group_ct]] == focal
    neighbour_idx <- df[[group_ct]] == neighbour

    # Compute expression ratios
    if (source_mode == "any") {
        other_types <- setdiff(unique(df[[group_ct]]), focal)
        other_means <- sapply(other_types, function(ct) {
            idx <- df[[group_ct]] == ct
            sapply(genes, function(g) mean(expr_mat[idx, g], na.rm = TRUE))
        })
    }

    # Build data frame
    plot_data <- lapply(genes, function(g) {
        res_c <- results_corrected[[g]]
        if (is.null(res_c) || isTRUE(res_c$error)) return(NULL)
        rv_c <- res_c$ran_vals
        if (is.null(rv_c)) return(NULL)
        row_c <- rv_c$level == focal & rv_c$term == neighbour
        if (!any(row_c)) return(NULL)

        slope_corr <- rv_c$estimate[row_c][1]
        se_corr <- rv_c$std.error[row_c][1]

        # Uncorrected slope if available
        slope_uncorr <- NA
        if (!is.null(results_uncorrected)) {
            res_u <- results_uncorrected[[g]]
            if (!is.null(res_u) && !isTRUE(res_u$error) && !is.null(res_u$ran_vals)) {
                rv_u <- res_u$ran_vals
                row_u <- rv_u$level == focal & rv_u$term == neighbour
                if (any(row_u)) slope_uncorr <- rv_u$estimate[row_u][1]
            }
        }

        focal_expr <- mean(expr_mat[focal_idx, g], na.rm = TRUE)

        if (source_mode == "neighbour") {
            neighbour_expr <- mean(expr_mat[neighbour_idx, g], na.rm = TRUE)
            expr_ratio <- neighbour_expr / (focal_expr + 0.01)
        } else {
            max_other <- max(other_means[g, ], na.rm = TRUE)
            expr_ratio <- max_other / (focal_expr + 0.01)
        }

        data.frame(
            gene = g,
            slope_corrected = slope_corr,
            slope_uncorrected = slope_uncorr,
            se = se_corr,
            expr_ratio = expr_ratio,
            is_spillover = expr_ratio > spill_threshold
        )
    }) |> purrr::compact() |> dplyr::bind_rows()

    if (nrow(plot_data) == 0) stop("No data to plot.")

    # Identify genes to label
    # Top spillover candidates (high ratio + positive corrected slope)
    spill_candidates <- plot_data |>
        dplyr::filter(is_spillover) |>
        dplyr::arrange(dplyr::desc(expr_ratio)) |>
        dplyr::slice_head(n = ceiling(n_label / 2))

    # Top biological genes (low ratio + large absolute slope)
    bio_candidates <- plot_data |>
        dplyr::filter(!is_spillover) |>
        dplyr::arrange(dplyr::desc(abs(slope_corrected))) |>
        dplyr::slice_head(n = ceiling(n_label / 2))

    label_genes <- unique(c(spill_candidates$gene, bio_candidates$gene))
    plot_data$label <- ifelse(plot_data$gene %in% label_genes, plot_data$gene, NA)

    # Build plot
    pair_label <- paste0(gsub("_", " ", focal), " \u2192 ", gsub("_", " ", neighbour))

    p <- ggplot(plot_data, aes(x = slope_corrected, y = expr_ratio))

    # If uncorrected available, show arrows
    if (!is.null(results_uncorrected) && any(!is.na(plot_data$slope_uncorrected))) {
        arrow_data <- plot_data |> dplyr::filter(!is.na(slope_uncorrected))
        p <- p +
            geom_segment(
                data = arrow_data,
                aes(x = slope_uncorrected, xend = slope_corrected,
                    y = expr_ratio, yend = expr_ratio),
                arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                colour = "grey70", linewidth = 0.3, alpha = 0.4
            )
    }

    p <- p +
        geom_hline(yintercept = spill_threshold, linetype = "dotted",
                   colour = "#B2182B", linewidth = 0.4, alpha = 0.6) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   colour = "grey50", linewidth = 0.3) +
        geom_point(aes(colour = is_spillover), size = 1.5, alpha = 0.6) +
        ggrepel::geom_text_repel(
            aes(label = label),
            size = 3, fontface = "italic",
            max.overlaps = 20, segment.size = 0.2,
            segment.color = "grey50", box.padding = 0.3,
            min.segment.length = 0
        ) +
        scale_colour_manual(
            values = c("FALSE" = "grey40", "TRUE" = "#D6604D"),
            labels = c("FALSE" = paste0("Ratio \u2264 ", spill_threshold),
                       "TRUE" = paste0("Ratio > ", spill_threshold)),
            name = "Spillover candidate"
        ) +
        scale_y_continuous(trans = "log1p",
                           breaks = c(0, 1, 2, 5, 10, 20, 50, 100))

    y_label <- if (source_mode == "any") {
        "Expression ratio (max other / focal)"
    } else {
        paste0("Expression ratio (", gsub("_", " ", neighbour), " / ", gsub("_", " ", focal), ")")
    }

    p <- p +
        labs(
            x = "Spatial slope (corrected)",
            y = y_label,
            title = paste0("Spillover landscape: ", pair_label),
            caption = if (!is.null(results_uncorrected))
                "Arrows show correction: uncorrected \u2192 corrected" else NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0),
            panel.grid.minor = element_blank(),
            legend.position = c(0.15, 0.9),
            legend.background = element_rect(fill = "white", colour = NA, linewidth = 0),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            plot.caption = element_text(size = 8, colour = "grey50"),
            plot.margin = margin(10, 15, 10, 10)
        )

    p
}

#' Pair-Specific Correction Heatmap
#'
#' For a given focal-neighbour pair, shows the spatial slope for each gene
#' before and after correction. Genes selected by largest absolute change.
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type
#' @param n_genes Number of top genes to show (default: 20)
#' @param title Plot title (auto-generated if NULL)
#' @return ggplot object
#' @export
plot_pair_correction <- function(
        results_corrected, results_uncorrected,
        focal, neighbour,
        n_genes = 20,
        title = NULL
) {
    library(ggplot2)

    genes <- intersect(names(results_corrected), names(results_uncorrected))

    # Extract slopes for this pair
    extract_pair_slopes <- function(results) {
        lapply(genes, function(g) {
            res <- results[[g]]
            if (is.null(res) || isTRUE(res$error) || is.null(res$ran_vals)) return(NULL)
            rv <- res$ran_vals
            row <- rv$level == focal & rv$term == neighbour
            if (!any(row)) return(NULL)
            data.frame(gene = g, slope = rv$estimate[row][1])
        }) |> purrr::compact() |> dplyr::bind_rows()
    }

    corr <- extract_pair_slopes(results_corrected) |>
        dplyr::rename(slope_corr = slope)
    uncorr <- extract_pair_slopes(results_uncorrected) |>
        dplyr::rename(slope_uncorr = slope)

    combined <- dplyr::inner_join(corr, uncorr, by = "gene") |>
        dplyr::mutate(delta = slope_corr - slope_uncorr)

    # Top genes by absolute delta
    top <- combined |>
        dplyr::arrange(dplyr::desc(abs(delta))) |>
        dplyr::slice_head(n = n_genes)

    # Long format
    plot_data <- top |>
        tidyr::pivot_longer(
            cols = c(slope_uncorr, slope_corr),
            names_to = "model",
            values_to = "slope"
        ) |>
        dplyr::mutate(
            model = dplyr::case_when(
                model == "slope_uncorr" ~ "Uncorrected",
                model == "slope_corr" ~ "Corrected"
            ),
            model = factor(model, levels = c("Uncorrected", "Corrected"))
        )

    # Order genes by uncorrected slope
    gene_order <- top |> dplyr::arrange(slope_uncorr) |> dplyr::pull(gene)
    plot_data$gene <- factor(plot_data$gene, levels = gene_order)

    pair_label <- paste0(gsub("_", " ", focal), " \u2192 ", gsub("_", " ", neighbour))
    if (is.null(title)) title <- paste0("Slope correction: ", pair_label)

    # Arrow data
    arrow_data <- top |>
        dplyr::mutate(gene = factor(gene, levels = gene_order))

    p <- ggplot() +
        # Arrows from uncorrected to corrected
        geom_segment(
            data = arrow_data,
            aes(x = slope_uncorr, xend = slope_corr, y = gene, yend = gene),
            arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
            colour = "grey50", linewidth = 0.4
        ) +
        # Uncorrected dots
        geom_point(
            data = plot_data |> dplyr::filter(model == "Uncorrected"),
            aes(x = slope, y = gene),
            shape = 21, size = 3.5, fill = "#D6604D", colour = "grey30", stroke = 0.4
        ) +
        # Corrected dots
        geom_point(
            data = plot_data |> dplyr::filter(model == "Corrected"),
            aes(x = slope, y = gene),
            shape = 21, size = 3.5, fill = "#4393C3", colour = "grey30", stroke = 0.4
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
        labs(
            x = "Spatial slope",
            y = NULL,
            title = title
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0),
            axis.text.y = element_text(face = "italic", size = 10),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
            plot.margin = margin(10, 15, 10, 10)
        )

    p
}


#' Multi-Pair Correction Figure
#'
#' Combines pair-specific correction plots for multiple pairs.
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results
#' @param pairs List of c(focal, neighbour) vectors
#' @param n_genes Number of genes per pair (default: 15)
#' @param ncol Number of columns in layout (default: length(pairs))
#' @return patchwork object
#' @export
plot_multi_pair_correction <- function(
        results_corrected, results_uncorrected,
        pairs,
        n_genes = 15,
        ncol = NULL
) {
    library(patchwork)

    plots <- lapply(pairs, function(pr) {
        plot_pair_correction(
            results_corrected, results_uncorrected,
            focal = pr[1], neighbour = pr[2],
            n_genes = n_genes
        )
    })

    if (is.null(ncol)) ncol <- length(plots)

    wrap_plots(plots, ncol = ncol) +
        plot_annotation(
            caption = expression(
                paste(phantom() %<-% phantom(),
                      "  arrow: uncorrected \u2192 corrected;   ",
                      scriptstyle("\u25CF"), " uncorrected   ",
                      scriptstyle("\u25CF"), " corrected")
            ),
            theme = theme(plot.caption = element_text(size = 9, colour = "grey50"))
        )
}


#' Pair-Specific Correction Waterfall
#'
#' Horizontal bars showing delta slope (corrected - uncorrected) for top genes.
#' Bars left = contamination removed, bars right = biology revealed.
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type
#' @param n_genes Number of top genes (default: 20)
#' @param r_uncorr Optional r value to annotate (uncorrected)
#' @param r_corr Optional r value to annotate (corrected)
#' @param title Plot title (auto-generated if NULL)
#' @return ggplot object
#' @export
plot_pair_waterfall <- function(
        results_corrected, results_uncorrected,
        focal, neighbour,
        n_genes = 20,
        r_uncorr = NULL, r_corr = NULL,
        title = NULL
) {
    library(ggplot2)

    genes <- intersect(names(results_corrected), names(results_uncorrected))

    extract_pair_slopes <- function(results) {
        lapply(genes, function(g) {
            res <- results[[g]]
            if (is.null(res) || isTRUE(res$error) || is.null(res$ran_vals)) return(NULL)
            rv <- res$ran_vals
            row <- rv$level == focal & rv$term == neighbour
            if (!any(row)) return(NULL)
            data.frame(gene = g, slope = rv$estimate[row][1])
        }) |> purrr::compact() |> dplyr::bind_rows()
    }

    corr <- extract_pair_slopes(results_corrected) |> dplyr::rename(slope_corr = slope)
    uncorr <- extract_pair_slopes(results_uncorrected) |> dplyr::rename(slope_uncorr = slope)

    combined <- dplyr::inner_join(corr, uncorr, by = "gene") |>
        dplyr::mutate(delta = slope_corr - slope_uncorr)

    top <- combined |>
        dplyr::arrange(dplyr::desc(abs(delta))) |>
        dplyr::slice_head(n = n_genes)

    # Order by delta value
    top <- top |>
        dplyr::arrange(delta) |>
        dplyr::mutate(
            gene = factor(gene, levels = gene),
            direction = ifelse(delta < 0, "Contamination removed", "Biology revealed")
        )

    pair_label <- paste0(gsub("_", " ", focal), " \u2192 ", gsub("_", " ", neighbour))
    if (is.null(title)) title <- pair_label

    # Subtitle with r values
    subtitle <- NULL
    if (!is.null(r_uncorr) && !is.null(r_corr)) {
        subtitle <- bquote(italic(r) ~ "=" ~ .(sprintf("%.2f", r_uncorr)) ~
                               "\u2192" ~ .(sprintf("%.2f", r_corr)))
    }

    p <- ggplot(top, aes(x = delta, y = gene, fill = direction)) +
        geom_col(width = 0.7) +
        geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
        scale_fill_manual(
            values = c("Contamination removed" = "#4393C3",
                       "Biology revealed" = "#D6604D"),
            name = NULL
        ) +
        labs(
            x = expression(Delta ~ "slope (corrected \u2212 uncorrected)"),
            y = NULL,
            title = title,
            subtitle = subtitle
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0),
            plot.subtitle = element_text(size = 11, colour = "#555555"),
            axis.text.y = element_text(face = "italic", size = 10),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size = 10),
            plot.margin = margin(10, 15, 10, 10)
        )

    p
}


#' Multi-Pair Waterfall Figure
#'
#' @param results_corrected Corrected model results
#' @param results_uncorrected Uncorrected model results
#' @param pairs List of c(focal, neighbour) vectors
#' @param r_values Optional list of c(r_uncorr, r_corr) per pair
#' @param n_genes Number of genes per pair (default: 15)
#' @param ncol Number of columns (default: length(pairs))
#' @return patchwork object
#' @export
plot_multi_pair_waterfall <- function(
        results_corrected, results_uncorrected,
        pairs,
        r_values = NULL,
        n_genes = 15,
        ncol = NULL
) {
    library(patchwork)

    plots <- lapply(seq_along(pairs), function(i) {
        pr <- pairs[[i]]
        rv <- if (!is.null(r_values)) r_values[[i]] else NULL
        plot_pair_waterfall(
            results_corrected, results_uncorrected,
            focal = pr[1], neighbour = pr[2],
            n_genes = n_genes,
            r_uncorr = rv[1], r_corr = rv[2]
        )
    })

    if (is.null(ncol)) ncol <- length(plots)
    wrap_plots(plots, ncol = ncol)
}

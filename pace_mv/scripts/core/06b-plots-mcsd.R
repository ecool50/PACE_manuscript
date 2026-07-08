# ==============================================================================
# 06b-plots-mcsd.R
# MCSD gene scoring and gene-level plots
# ==============================================================================

plot_mcsd_top <- function(mcsd_results,
                          top_n = 5,
                          title = NULL,
                          fill = "darkred",
                          label_accuracy = 0.001,
                          y_expand = 1.1,
                          ylab_prefix = "MCSD",
                          block_label = NULL) {
    stopifnot(is.list(mcsd_results), "scores" %in% names(mcsd_results))
    scores <- mcsd_results$scores
    if (!all(c("gene","MCSD") %in% names(scores))) {
        stop("mcsd_results$scores must contain columns 'gene' and 'MCSD'.")
    }

    # Pick a default title from meta if available
    if (is.null(title) && "meta" %in% names(mcsd_results)) {
        meta <- mcsd_results$meta
        foc  <- if (!is.null(meta$focal) && nzchar(meta$focal)) meta$focal else NULL
        nb   <- if (!is.null(meta$neighbour) && length(meta$neighbour)) {
            paste(meta$neighbour, collapse = ", ")
        } else NULL
        if (!is.null(foc) && !is.null(nb)) {
            title <- paste(foc, "-", nb)
        } else if (!is.null(foc)) {
            title <- foc
        }
    }
    if (is.null(title)) title <- ""

    # Y-axis label (optionally includes block)
    if (is.null(block_label) && "meta" %in% names(mcsd_results) && !is.null(mcsd_results$meta$block)) {
        block_label <- mcsd_results$meta$block
    }
    y_lab <- if (!is.null(block_label) && nzchar(block_label)) {
        sprintf("%s (%s)", ylab_prefix, block_label)
    } else {
        ylab_prefix
    }

    # Top N genes
    df_top <- scores %>%
        dplyr::as_tibble() %>%
        dplyr::slice_max(order_by = .data$MCSD, n = top_n, with_ties = FALSE)

    # Order bars by MCSD descending (top at top after coord_flip)
    df_top <- df_top %>%
        dplyr::mutate(gene = factor(.data$gene,
                                    levels = rev(.data$gene[order(.data$MCSD, decreasing = TRUE)])))

    ymax <- max(df_top$MCSD, na.rm = TRUE)
    if (!is.finite(ymax)) ymax <- 1

    ggplot2::ggplot(df_top, ggplot2::aes(x = .data$gene, y = .data$MCSD)) +
        ggplot2::geom_col(fill = fill, width = 0.8) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = label_accuracy)) +
        ggplot2::labs(x = NULL, y = y_lab, title = title) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::geom_text(
            ggplot2::aes(label = scales::number(.data$MCSD, accuracy = label_accuracy)),
            hjust = -0.1, size = 3
        ) +
        ggplot2::expand_limits(y = ymax * y_expand)
}

plot_raw_boxlines_gene <- function(
        df, gene,
        focal = "Fibroblast", neighbour = "Tumour",
        cov_name = "Tumour",            # numeric covariate to bin
        nbins = 5,
        fill_color = "forestgreen",
        line_color = "forestgreen",
        used_rows = NULL,               # optional: rownames/indices to use
        out_by_gene = NULL,             # optional: only to derive used_rows
        fit_light = NULL,               # fit_light object for slope overlay
        ran_vals = NULL,                # tidy ran_vals for annotation
        show_slope_line = FALSE,        # overlay fitted slope on link scale
        show_annotation = FALSE,        # annotate with BLUP estimate & p-value
        group_ct = "celltype"
) {
    # derive used_rows if not provided
    if (is.null(used_rows)) {
        if (!is.null(fit_light)) {
            used_rows <- fit_light$used_rows
        } else if (!is.null(out_by_gene)) {
            og <- out_by_gene[[gene]]
            if (is.null(og) || is.null(rownames(og$out)))
                stop("out_by_gene[[gene]]$out must exist and have rownames, or pass used_rows.")
            used_rows <- rownames(og$out)
        } else {
            used_rows <- rownames(df)
            if (is.null(used_rows)) used_rows <- seq_len(nrow(df))
        }
    }

    # slice the data
    df_loc <- df[used_rows, , drop = FALSE]

    # checks
    if (!gene %in% colnames(df_loc)) stop("Gene '", gene, "' not found in df.")
    if (!cov_name %in% colnames(df_loc)) stop("cov_name '", cov_name, "' not in df.")
    idx <- which(df_loc[[group_ct]] == focal)
    if (!length(idx)) stop("No rows for focal = '", focal, "'.")

    # raw expression & covariate
    expr <- as.numeric(df_loc[idx, gene, drop = TRUE])
    x    <- df_loc[[cov_name]][idx]

    keep <- is.finite(expr) & is.finite(x)
    expr <- expr[keep]; x <- x[keep]
    if (!length(expr)) stop("No finite values after filtering NA/Inf.")

    # bin covariate (adaptive quantiles if many unique values)
    if (length(unique(x)) > nbins) {
        br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
        if (length(br) <= 2) br <- pretty(x, nbins)
    } else {
        br <- pretty(x, nbins)
    }
    bin <- cut(x, breaks = br, include.lowest = TRUE, ordered_result = TRUE)

    df_plot <- tibble::tibble(bin = bin, expr = expr, x_raw = x)
    mean_tbl <- df_plot |>
        dplyr::group_by(bin) |>
        dplyr::summarise(mu = mean(expr, na.rm = TRUE),
                         x_mid = mean(x_raw, na.rm = TRUE),
                         .groups = "drop")

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = bin, y = expr)) +
        ggplot2::geom_boxplot(width = 0.65, outlier.alpha = 0.2, fatten = 1.1, color = fill_color) +
        ggplot2::geom_point(
            data = mean_tbl,
            ggplot2::aes(x = bin, y = mu),
            color = line_color, size = 3, inherit.aes = FALSE
        ) +
        ggplot2::labs(
            title = paste0(gene, " expression"),
            x = paste0("Number of ", neighbour, " neighbours"),
            y = paste0(focal, " ", gene, " expression")
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )

    # Overlay fitted slope line from model
    if (show_slope_line && !is.null(fit_light)) {
        re_list <- fit_light$ranef
        beta0 <- NULL; b_slope <- NULL

        if (group_ct %in% names(re_list)) {
            tab <- re_list[[group_ct]]
            levs <- rownames(tab); terms <- colnames(tab)

            if (focal %in% levs) {
                # Intercept BLUP
                if ("(Intercept)" %in% terms) {
                    beta0 <- fit_light$fixef["(Intercept)"] + tab[focal, "(Intercept)"]
                }
                # Neighbour slope BLUP
                if (cov_name %in% terms) {
                    b_slope <- tab[focal, cov_name]
                }
            }
        }

        if (!is.null(beta0) && !is.null(b_slope)) {
            # Generate fitted line on response scale across bin midpoints
            x_seq <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100)
            eta <- as.numeric(beta0) + as.numeric(b_slope) * x_seq
            mu_fitted <- exp(eta)  # log link -> response

            # Map continuous x to bin positions for overlay
            bin_seq <- cut(x_seq, breaks = br, include.lowest = TRUE, ordered_result = TRUE)
            bin_num <- as.numeric(bin_seq)
            # Fractional position within bin for smooth line
            bin_frac <- bin_num + (x_seq - tapply(x, bin, min, na.rm = TRUE)[bin_seq]) /
                pmax(tapply(x, bin, function(v) diff(range(v, na.rm = TRUE)))[bin_seq], 1e-6) - 0.5

            # Simpler: use numeric position based on bin index
            df_line <- tibble::tibble(bin = bin_seq, y = mu_fitted) |>
                dplyr::filter(!is.na(bin))

            p <- p + ggplot2::geom_line(
                data = df_line,
                ggplot2::aes(x = bin, y = y, group = 1),
                color = "red", linewidth = 1.2, inherit.aes = FALSE
            )
        }
    }

    # Annotation with BLUP estimate and p-value
    if (show_annotation && !is.null(ran_vals)) {
        rv_nb <- ran_vals |>
            dplyr::filter(
                group == group_ct,
                level == focal,
                term == cov_name
            )

        if (nrow(rv_nb) > 0) {
            est <- rv_nb$estimate[1]
            pv  <- rv_nb$pval[1]
            pv_str <- if (pv < 0.001) "< 0.001" else sprintf("%.3f", pv)
            label <- sprintf("b = %.3f, p %s", est, pv_str)

            p <- p + ggplot2::annotate(
                "text",
                x = nlevels(bin) * 0.7,
                y = max(expr, na.rm = TRUE) * 0.95,
                label = label,
                size = 3, hjust = 0.5, color = "red"
            )
        }
    }

    p
}

#' Plot Raw Boxlines for Top MCSD Genes of a Focal-Neighbour Pair
#'
#' Computes MCSD for a given pair, then generates boxline plots for
#' the top N genes using plot_raw_boxlines_gene.
#'
#' @param fits_light Named list of fit_light objects
#' @param results Model results (list of per-gene results with ran_vals)
#' @param df Data frame
#' @param focal Focal cell type (underscore format, e.g. "B_Cell")
#' @param neighbour Neighbour cell type (underscore format, e.g. "T_Cell")
#' @param n_genes Number of top MCSD genes to plot (default: 5)
#' @param nbins Number of bins for boxline plot (default: 4)
#' @param spill_correct Spillover correction mode (default: "weight")
#' @param spill_source_mode Source mode (default: "any")
#' @param detection_weight Detection weighting (default: "post")
#' @param min_expr_prev Minimum prevalence (default: 0.10)
#' @param min_expr_var Minimum variance (default: 1.0)
#' @param ncol Number of columns in grid layout (default: NULL, auto)
#' @param BPPARAM BiocParallel backend
#' @return patchwork object with arranged plots
#' @export
plot_top_mcsd_boxlines <- function(
        fits_light, results, df,
        focal, neighbour,
        n_genes = 5,
        nbins = 4,
        spill_correct = "weight",
        spill_source_mode = "any",
        detection_weight = "post",
        min_expr_prev = 0.10,
        min_expr_var = 1.0,
        show_slope_line = FALSE,
        show_annotation = TRUE,
        ncol = NULL,
        BPPARAM = BiocParallel::SerialParam()
) {
    library(patchwork)

    # Compute MCSD for this pair
    message("Computing MCSD for ", focal, " ~ ", neighbour)
    mcsd_result <- mcsd_block_from_fits_no_resp(
        fits_light, df,
        block = "Spatial cell state",
        focal = focal, neighbour = neighbour,
        center = TRUE,
        min_expr_prev = min_expr_prev,
        min_expr_var = min_expr_var,
        spill_correct = spill_correct,
        spill_source_mode = spill_source_mode,
        detection_weight = detection_weight,
        BPPARAM = BPPARAM
    )

    top_genes <- mcsd_result$scores |>
        dplyr::slice_head(n = n_genes)

    message("Top ", n_genes, " genes: ", paste(top_genes$gene, collapse = ", "))

    # Extract ran_vals from results
    all_rv <- lapply(names(results), function(g) {
        res <- results[[g]]
        if (is.null(res) || isTRUE(res$error) || is.null(res$ran_vals)) return(NULL)
        res$ran_vals
    }) |> purrr::compact() |> dplyr::bind_rows()

    # Generate plots
    plots <- lapply(seq_len(nrow(top_genes)), function(i) {
        g <- top_genes$gene[i]
        mcsd_val <- top_genes$MCSD[i]

        fl <- fits_light[[g]]
        rv <- all_rv |> dplyr::filter(gene == g)

        if (is.null(fl) || nrow(rv) == 0) {
            message("  Skipping ", g, " (no fit or ran_vals)")
            return(NULL)
        }

        tryCatch({
            p <- plot_raw_boxlines_gene(
                df = df,
                gene = g,
                focal = focal,
                neighbour = neighbour,
                cov_name = neighbour,
                nbins = nbins,
                fit_light = fl,
                ran_vals = rv,
                show_slope_line = show_slope_line,
                show_annotation = show_annotation
            ) +
                ggplot2::labs(
                    subtitle = paste0(g, "  (MCSD = ", sprintf("%.3f", mcsd_val), ")")
                ) +
                ggplot2::theme(
                    plot.subtitle = ggplot2::element_text(face = "bold.italic", size = 10)
                )
            p
        }, error = function(e) {
            message("  Error plotting ", g, ": ", conditionMessage(e))
            NULL
        })
    })

    plots <- Filter(Negate(is.null), plots)

    if (length(plots) == 0) stop("No plots generated.")

    if (is.null(ncol)) ncol <- min(length(plots), 3)

    wrap_plots(plots, ncol = ncol) +
        plot_annotation(
            title = paste0(gsub("_", " ", focal), " \u2192 ", gsub("_", " ", neighbour)),
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 14)
            )
        )
}

# ==============================================================================
# Responder Interaction Boxlines
# ==============================================================================

#' Raw Boxlines for a Single Gene, Split by Responder
#'
#' Bins neighbour count, shows boxplots coloured by responder status
#' with group mean lines. Averages within imageID first.
#'
#' @param df Data frame
#' @param gene Gene name
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type (covariate column name)
#' @param resp_var Responder column name
#' @param sample_id Sample ID column
#' @param nbins Number of bins
#' @param fit_light Optional fit_light for used_rows
#' @param ran_vals Optional ran_vals for annotation
#' @param resp_term Responder term for annotation lookup
#' @param show_annotation Show BLUP estimates
#' @param group_ct Cell type column
#' @param palette Colour palette
#' @return ggplot object
plot_raw_boxlines_gene_resp <- function(
        df, gene,
        focal = "Fibroblast", neighbour = "Tumour",
        resp_var = "Responder",
        sample_id = "imageID",
        nbins = 4,
        breaks = NULL,
        fit_light = NULL,
        ran_vals = NULL,
        resp_term = "ResponderPD",
        show_annotation = TRUE,
        show_lines = TRUE,
        group_ct = "celltype",
        palette = c("SD" = "#4393C3", "PD" = "#D6604D")
) {
    # Derive used_rows
    used_rows <- if (!is.null(fit_light)) fit_light$used_rows else {
        rn <- rownames(df); if (is.null(rn)) seq_len(nrow(df)) else rn
    }

    df_loc <- df[used_rows, , drop = FALSE]

    if (!gene %in% colnames(df_loc)) stop("Gene '", gene, "' not found in df.")
    if (!neighbour %in% colnames(df_loc)) stop("Neighbour '", neighbour, "' not in df.")

    idx <- which(df_loc[[group_ct]] == focal)
    if (!length(idx)) stop("No rows for focal = '", focal, "'.")

    df_focal <- df_loc[idx, , drop = FALSE]

    expr <- as.numeric(df_focal[[gene]])
    x    <- as.numeric(df_focal[[neighbour]])
    resp <- df_focal[[resp_var]]
    sid  <- df_focal[[sample_id]]

    keep <- is.finite(expr) & is.finite(x) & !is.na(resp) & !is.na(sid)
    df_work <- data.frame(
        expr = expr[keep], x = x[keep],
        resp = resp[keep], sid = sid[keep]
    )
    if (nrow(df_work) == 0) stop("No finite values after filtering.")

    # Bin neighbour count
    if (!is.null(breaks)) {
        br <- breaks
    } else if (length(unique(df_work$x)) > nbins) {
        br <- unique(stats::quantile(df_work$x, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
        if (length(br) <= 2) br <- pretty(df_work$x, nbins)
    } else {
        br <- pretty(df_work$x, nbins)
    }
    df_work$bin <- cut(df_work$x, breaks = br, include.lowest = TRUE, ordered_result = TRUE)
    df_work <- df_work[!is.na(df_work$bin), ]

    # Average within sample x bin x responder
    df_agg <- df_work |>
        dplyr::group_by(bin, resp, sid) |>
        dplyr::summarise(expr_mean = mean(expr, na.rm = TRUE), .groups = "drop")

    # Group means for lines
    df_mean <- df_agg |>
        dplyr::group_by(bin, resp) |>
        dplyr::summarise(mu = mean(expr_mean, na.rm = TRUE), .groups = "drop")

    pd <- ggplot2::position_dodge(width = 0.75)

    p <- ggplot2::ggplot(df_agg, ggplot2::aes(x = bin, y = expr_mean, fill = resp)) +
        ggplot2::geom_boxplot(
            position = pd, width = 0.6,
            outlier.shape = NA, alpha = 0.5, colour = "#555555"
        ) +
        ggplot2::geom_point(
            ggplot2::aes(colour = resp),
            position = ggplot2::position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
            size = 1.2, alpha = 0.4, show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(values = palette, name = resp_var) +
        ggplot2::scale_colour_manual(values = palette, name = resp_var) +
        ggplot2::labs(
            title = paste0(gene, " expression"),
            x = paste0("Number of ", gsub("_", " ", neighbour), " neighbours"),
            y = paste0(gsub("_", " ", focal), " ", gene, " expression")
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
            legend.position = "right"
        )

    # Trend lines
    if (show_lines) {
        p <- p +
            ggplot2::geom_line(
                data = df_mean,
                ggplot2::aes(x = bin, y = mu, colour = resp, group = resp),
                position = pd, linewidth = 0.9
            ) +
            ggplot2::geom_point(
                data = df_mean,
                ggplot2::aes(x = bin, y = mu, colour = resp),
                position = pd, shape = 23, size = 2.5, fill = "white", stroke = 0.6
            )
    }

    # Annotation: base slope + responder interaction slope
    if (show_annotation && !is.null(ran_vals)) {
        # Base slope (shared)
        rv_base <- ran_vals |>
            dplyr::filter(group == group_ct, level == focal, term == neighbour)
        # Responder interaction slope
        resp_nb <- paste0(resp_term, ":", neighbour)
        rv_resp <- ran_vals |>
            dplyr::filter(group == group_ct, level == focal, term == resp_nb)

        labels <- c()
        if (nrow(rv_base) > 0) {
            pv <- rv_base$pval[1]
            pv_str <- if (pv < 0.001) "< 0.001" else sprintf("%.3f", pv)
            labels <- c(labels, sprintf("base: b = %.3f, p %s", rv_base$estimate[1], pv_str))
        }
        if (nrow(rv_resp) > 0) {
            pv <- rv_resp$pval[1]
            pv_str <- if (pv < 0.001) "< 0.001" else sprintf("%.3f", pv)
            labels <- c(labels, sprintf("PD\u00D7nb: b = %.3f, p %s", rv_resp$estimate[1], pv_str))
        }

        if (length(labels)) {
            p <- p + ggplot2::annotate(
                "text",
                x = nlevels(df_work$bin) * 0.65,
                y = max(df_agg$expr_mean, na.rm = TRUE) * 0.97,
                label = paste(labels, collapse = "\n"),
                size = 3.5, hjust = 0.5, color = "#333333"
            )
        }
    }

    p
}


#' Top MCSD Boxlines with Responder Split
#'
#' Computes MCSD for a given pair and block, then generates
#' responder-split boxline plots for the top N genes.
#'
#' @param fits_light Named list of fit_light objects
#' @param results Model results (with ran_vals)
#' @param df Data frame
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type
#' @param block Block for MCSD computation
#' @param n_genes Number of top genes to plot
#' @param nbins Number of bins
#' @param resp_var Responder column name
#' @param resp_term Responder term
#' @param resp_level Responder level
#' @param sample_id Sample ID column
#' @param spill_correct Spillover correction
#' @param spill_source_mode Source mode
#' @param detection_weight Detection weighting
#' @param palette Colour palette
#' @param ncol Number of columns
#' @param BPPARAM BiocParallel backend
#' @return patchwork object
plot_top_mcsd_boxlines_resp <- function(
        fits_light, results, df,
        focal, neighbour,
        block = "Responder spatial state",
        n_genes = 5,
        nbins = 4,
        breaks = NULL,
        resp_var = "Responder",
        resp_term = "ResponderPD",
        resp_level = "PD",
        sample_id = "imageID",
        spill_correct = "weight",
        spill_source_mode = "any",
        detection_weight = "post",
        min_expr_prev = 0.10,
        min_expr_var = 1.0,
        palette = c("SD" = "#4393C3", "PD" = "#D6604D"),
        show_annotation = TRUE,
        show_lines = TRUE,
        ncol = NULL,
        BPPARAM = BiocParallel::SerialParam()
) {
    library(patchwork)

    message("Computing MCSD for ", focal, " ~ ", neighbour, " (", block, ")")
    mcsd_result <- mcsd_block_from_fits(
        fits_light, df,
        block = block,
        focal = focal, neighbour = neighbour,
        resp_term = resp_term, resp_level = resp_level,
        center = TRUE,
        min_expr_prev = min_expr_prev,
        min_expr_var = min_expr_var,
        spill_correct = spill_correct,
        spill_source_mode = spill_source_mode,
        detection_weight = detection_weight,
        BPPARAM = BPPARAM
    )

    top_genes <- mcsd_result$scores |>
        dplyr::slice_head(n = n_genes)

    message("Top ", n_genes, " genes: ", paste(top_genes$gene, collapse = ", "))

    # Extract ran_vals
    all_rv <- lapply(names(results), function(g) {
        res <- results[[g]]
        if (is.null(res) || isTRUE(res$error) || is.null(res$ran_vals)) return(NULL)
        res$ran_vals
    }) |> purrr::compact() |> dplyr::bind_rows()

    plots <- lapply(seq_len(nrow(top_genes)), function(i) {
        g <- top_genes$gene[i]
        mcsd_val <- top_genes$MCSD[i]

        fl <- fits_light[[g]]
        rv <- all_rv |> dplyr::filter(gene == g)

        if (is.null(fl) || nrow(rv) == 0) {
            message("  Skipping ", g, " (no fit or ran_vals)")
            return(NULL)
        }

        tryCatch({
            p <- plot_raw_boxlines_gene_resp(
                df = df, gene = g,
                focal = focal, neighbour = neighbour,
                resp_var = resp_var,
                sample_id = sample_id,
                nbins = nbins,
                breaks = breaks,
                fit_light = fl,
                ran_vals = rv,
                resp_term = resp_term,
                show_annotation = show_annotation,
                show_lines = show_lines,
                palette = palette
            ) +
                ggplot2::labs(
                    subtitle = paste0(g, "  (MCSD = ", sprintf("%.3f", mcsd_val), ")")
                ) +
                ggplot2::theme(
                    plot.subtitle = ggplot2::element_text(face = "bold.italic", size = 10)
                )
            p
        }, error = function(e) {
            message("  Error plotting ", g, ": ", conditionMessage(e))
            NULL
        })
    })

    plots <- Filter(Negate(is.null), plots)
    if (length(plots) == 0) stop("No plots generated.")
    if (is.null(ncol)) ncol <- min(length(plots), 3)

    wrap_plots(plots, ncol = ncol) +
        plot_annotation(
            title = paste0(gsub("_", " ", focal), " \u2192 ", gsub("_", " ", neighbour),
                           "  (", block, ")"),
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", size = 14)
            )
        )
}

#' Plot per-gene block decomposition for top MCSD genes
#'
#' @param bsm Output of block_sum_multigene_from_fits() for a specific focal/neighbour
#' @param mcsd_scores Output of mcsd_block_from_fits()$scores (or mcsd_block_from_fits_no_resp()$scores)
#' @param top_n Number of top MCSD genes to show (default: 20)
#' @param title Optional title override
#' @return ggplot object
plot_block_decomp_top_mcsd <- function(bsm, mcsd_scores, top_n = 20, title = NULL) {

    library(ggplot2)
    library(dplyr)

    # Get top genes by MCSD
    top_genes <- head(mcsd_scores$gene, top_n)

    # Filter per-gene decomposition to top genes
    # Compute proportion of total (including residuals) but don't display residuals
    pg <- bsm$per_gene |>
        filter(gene %in% top_genes) |>
        group_by(gene) |>
        mutate(prop = SS_link / sum(SS_link)) |>
        ungroup() |>
        filter(block != "Residuals")

    # Order genes by MCSD rank (left = highest)
    pg$gene <- factor(pg$gene, levels = top_genes)

    # Clean up block names and set order
    block_order <- c("Cell type", "Spatial cell state",
                     "Responder spatial state", "Responder status", "Spillover")
    pg$block <- factor(pg$block, levels = block_order)

    block_cols <- c(
        "Cell type"                = "#003049",
        "Spatial cell state"       = "#780000",
        "Responder spatial state"  = "#c1121f",
        "Responder status"         = "#669bbc",
        "Spillover"                = "#E5A100"
    )

    if (is.null(title)) {
        title <- paste("Block decomposition \u2014 top", top_n, "MCSD genes")
    }

    keep_levels <- c("Cell type", "Spatial cell state",
                     "Responder spatial state", "Responder status", "Spillover")

    p <- ggplot(pg, aes(x = gene, y = prop, fill = block)) +
        geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_fill_manual(values = block_cols, breaks = keep_levels,
                          limits = keep_levels, drop = FALSE) +
        labs(
            x = NULL,
            y = "Proportion of total",
            title = title,
            fill = "Block"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.title = element_text(hjust = 0.5)
        )

    p
}


#' Convenience wrapper: compute MCSD + block decomp and plot
#'
#' @param fits_light Named list of fit_light objects
#' @param df Data frame with cell metadata (must include _near columns)
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type
#' @param block MCSD block to rank by (default: "Responder spatial state")
#' @param top_n Number of top genes (default: 20)
#' @param resp_term Responder term name
#' @param group_ct Cell type grouping column
#' @param spill_correct Spillover correction method for MCSD
#' @param spill_source_mode Spillover source mode for MCSD
#' @param detection_weight Detection weight mode for MCSD
#' @param BPPARAM BiocParallel backend
#' @return List with plot, bsm, and mcsd_scores
plot_top_mcsd_block_decomp <- function(fits_light, df,
                                       focal, neighbour,
                                       block = "Responder spatial state",
                                       top_n = 20,
                                       resp_term = "ResponderPD",
                                       group_ct = "celltype",
                                       spill_correct = "weight",
                                       spill_source_mode = "any",
                                       detection_weight = "post",
                                       BPPARAM = BiocParallel::SerialParam()) {

    # 1. Compute MCSD scores
    message("Computing MCSD (", block, ") for ", focal, " -> ", neighbour, "...")
    mcsd_out <- mcsd_block_from_fits(
        fits_light = fits_light,
        df = df,
        block = block,
        focal = focal,
        neighbour = neighbour,
        group_ct = group_ct,
        resp_term = resp_term,
        spill_correct = spill_correct,
        spill_source_mode = spill_source_mode,
        detection_weight = detection_weight,
        BPPARAM = BPPARAM
    )

    # 2. Per-gene block decomposition
    message("Computing per-gene block decomposition...")
    bsm <- block_sum_multigene_from_fits(
        fits_light = fits_light,
        df = df,
        focal = focal,
        neighbour = neighbour,
        BPPARAM = BPPARAM
    )

    # 3. Plot
    title <- paste0(focal, " \u2192 ", neighbour, ": top ", top_n,
                    " genes (", block, ")")
    p <- plot_block_decomp_top_mcsd(bsm, mcsd_out$scores, top_n = top_n,
                                    title = title)

    list(plot = p, bsm = bsm, mcsd_scores = mcsd_out$scores)
}


# ==============================================================================
# Non-responder versions
# ==============================================================================

#' Plot per-gene block decomposition for top MCSD genes (non-responder)
#'
#' @param bsm Output of block_sum_multigene_from_fits() for a specific focal/neighbour
#' @param mcsd_scores Output of mcsd_block_from_fits_no_resp()$scores
#' @param top_n Number of top MCSD genes to show (default: 20)
#' @param title Optional title override
#' @return ggplot object
plot_block_decomp_top_mcsd_no_resp <- function(bsm, mcsd_scores, top_n = 20, title = NULL) {

    library(ggplot2)
    library(dplyr)

    top_genes <- head(mcsd_scores$gene, top_n)

    pg <- bsm$per_gene |>
        filter(gene %in% top_genes) |>
        group_by(gene) |>
        mutate(prop = SS_link / sum(SS_link)) |>
        ungroup() |>
        filter(block != "Residuals")

    pg$gene <- factor(pg$gene, levels = top_genes)

    block_order <- c("Cell type", "Spatial cell state", "Spillover")
    pg$block <- factor(pg$block, levels = block_order)

    block_cols <- c(
        "Cell type"          = "#003049",
        "Spatial cell state" = "#780000",
        "Spillover"          = "#E5A100"
    )

    if (is.null(title)) {
        title <- paste("Block decomposition \u2014 top", top_n, "MCSD genes")
    }

    keep_levels <- c("Cell type", "Spatial cell state", "Spillover")

    p <- ggplot(pg, aes(x = gene, y = prop, fill = block)) +
        geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_fill_manual(values = block_cols, breaks = keep_levels,
                          limits = keep_levels, drop = FALSE) +
        labs(
            x = NULL,
            y = "Proportion of total",
            title = title,
            fill = "Block"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.title = element_text(hjust = 0.5)
        )

    p
}


#' Convenience wrapper: compute MCSD + block decomp and plot (non-responder)
#'
#' @param fits_light Named list of fit_light objects
#' @param df Data frame with cell metadata (must include _near columns)
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type
#' @param block MCSD block to rank by (default: "Spatial cell state")
#' @param top_n Number of top genes (default: 20)
#' @param group_ct Cell type grouping column
#' @param spill_correct Spillover correction method for MCSD
#' @param spill_source_mode Spillover source mode for MCSD
#' @param detection_weight Detection weight mode for MCSD
#' @param BPPARAM BiocParallel backend
#' @return List with plot, bsm, and mcsd_scores
plot_top_mcsd_block_decomp_no_resp <- function(fits_light, df,
                                               focal, neighbour,
                                               block = "Spatial cell state",
                                               top_n = 20,
                                               group_ct = "celltype",
                                               spill_correct = "weight",
                                               spill_source_mode = "any",
                                               detection_weight = "post",
                                               BPPARAM = BiocParallel::SerialParam()) {

    message("Computing MCSD (", block, ") for ", focal, " -> ", neighbour, "...")
    mcsd_out <- mcsd_block_from_fits_no_resp(
        fits_light = fits_light,
        df = df,
        block = block,
        focal = focal,
        neighbour = neighbour,
        group_ct = group_ct,
        spill_correct = spill_correct,
        spill_source_mode = spill_source_mode,
        detection_weight = detection_weight,
        BPPARAM = BPPARAM
    )

    message("Computing per-gene block decomposition...")
    bsm <- block_sum_multigene_from_fits(
        fits_light = fits_light,
        df = df,
        focal = focal,
        neighbour = neighbour,
        BPPARAM = BPPARAM
    )

    title <- paste0(focal, " \u2192 ", neighbour, ": top ", top_n,
                    " genes (", block, ")")
    p <- plot_block_decomp_top_mcsd_no_resp(bsm, mcsd_out$scores, top_n = top_n,
                                            title = title)

    list(plot = p, bsm = bsm, mcsd_scores = mcsd_out$scores)
}


#' Volcano Plot: Z-score vs MCSD
#'
#' Plots spatial slope z-score (x-axis) against MCSD importance score (y-axis)
#' for each gene. Highlights genes that are both strong spatial effects and
#' important contributors to the coordinated spatial program.
#'
#' @param drivers Output from \code{mcsd_top_drivers()} or any data frame with
#'   columns: gene, MCSD, z_score, direction
#' @param n_label Number of top genes to label (default: 10)
#' @param z_threshold Z-score threshold lines (default: 2)
#' @param mcsd_threshold MCSD threshold line (default: NULL, no line)
#' @param focal Focal cell type (for title)
#' @param neighbour Neighbour cell type (for title)
#' @param title Custom title (overrides focal/neighbour auto-title)
#' @param cols Named vector of colours for direction (default: up=red, down=blue)
#' @param point_size Point size (default: 2)
#' @param point_alpha Point alpha (default: 0.7)
#' @param base_size Base font size (default: 14)
#' @return A ggplot object
#' @export
plot_mcsd_volcano <- function(
        drivers,
        n_label = 10,
        z_threshold = 2,
        mcsd_threshold = NULL,
        focal = NULL,
        neighbour = NULL,
        title = NULL,
        cols = c("down" = "#003049", "up" = "#780000"),
        point_size = 2,
        point_alpha = 0.7,
        base_size = 14
) {
    stopifnot(all(c("gene", "MCSD", "z_score") %in% names(drivers)))

    # Add direction if not present
    if (!"direction" %in% names(drivers)) {
        drivers <- drivers |>
            dplyr::mutate(direction = ifelse(z_score > 0, "up", "down"))
    }

    # Auto-title
    if (is.null(title) && !is.null(focal) && !is.null(neighbour)) {
        title <- paste0(focal, " \u2192 ", neighbour, ": spatial drivers")
    }

    # Label top genes by MCSD
    top_genes <- drivers |>
        dplyr::arrange(dplyr::desc(MCSD)) |>
        dplyr::slice_head(n = n_label)

    drivers <- drivers |>
        dplyr::mutate(label = ifelse(gene %in% top_genes$gene, gene, NA_character_))

    p <- ggplot2::ggplot(drivers, ggplot2::aes(x = z_score, y = MCSD,
                                                colour = direction)) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggrepel::geom_text_repel(
            ggplot2::aes(label = label),
            colour = "black", size = 3.5,
            max.overlaps = 20, seed = 42,
            min.segment.length = 0.2
        ) +
        ggplot2::geom_vline(xintercept = c(-z_threshold, z_threshold),
                            linetype = "dashed", colour = "grey50", linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
        ggplot2::scale_colour_manual(
            values = cols,
            labels = c("down" = "Downregulated", "up" = "Upregulated"),
            name = "Direction"
        ) +
        ggplot2::labs(
            x = "Spatial slope z-score",
            y = "MCSD score",
            title = title
        ) +
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(hjust = 0.5)
        )

    # Optional MCSD threshold line
    if (!is.null(mcsd_threshold)) {
        p <- p + ggplot2::geom_hline(yintercept = mcsd_threshold,
                                      linetype = "dashed", colour = "grey50",
                                      linewidth = 0.4)
    }

    p
}

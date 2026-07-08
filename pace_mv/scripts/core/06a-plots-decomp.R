# ==============================================================================
# 06a-plots-decomp.R
# Variance decomposition visualization
# ==============================================================================

plot_stacked_bar_focal <- function(
        blks_all,
        block_levels = c(
            "Cell type",
            "Spatial cell state",
            "Responder spatial state",
            "Responder status",
            "Spillover",
            "Residuals"             # kept for completeness; excluded from stack
        ),
        title = "",
        zoom_blocks = c("Spatial cell state", "Responder spatial state", "Spillover")
) {
    library(ggplot2)

    # keep blocks (exclude Residuals from the stack)
    keep_levels <- setdiff(block_levels, "Residuals")

    # residual share per focal (useful for annotations if needed)
    resid_by_focal <- blks_all %>%
        dplyr::filter(.data$block == "Residuals") %>%
        dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " ")) %>%
        dplyr::transmute(focal, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

    # data for stacked bars (no renormalisation)
    blks_plot <- blks_all %>%
        dplyr::filter(.data$block != "Residuals") %>%
        dplyr::mutate(
            focal = stringr::str_replace(.data$focal, "_", " "),
            block = factor(.data$block, levels = keep_levels)
        ) %>%
        tidyr::complete(focal, block, fill = list(pct_total = 0)) %>%
        dplyr::left_join(resid_by_focal, by = "focal")

    # NPG colors
    cols <- c(
        "Cell type"               = "#003049",
        "Spatial cell state"      = "#780000",
        "Responder spatial state" = "#c1121f",
        "Responder status"        = "#669bbc",
        "Spillover"               = "#E5A100"
    )
    cols <- cols[keep_levels]

    # Shared focal ordering: by Spatial cell state descending
    focal_order <- blks_plot %>%
        dplyr::filter(block == "Spatial cell state") %>%
        dplyr::arrange(dplyr::desc(pct_total)) %>%
        dplyr::pull(focal)
    blks_plot$focal <- factor(blks_plot$focal, levels = focal_order)

    # --- Left panel: full stacked bar ---
    p_full <- ggplot(blks_plot, aes(x = focal, y = pct_total, fill = block)) +
        geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
        labs(
            x = "Focal cell type",
            y = "Proportion of total variance",
            fill = "Block",
            title = title
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )

    # --- Right panel: zoomed grouped bar (small blocks only) ---
    zoom_blocks <- intersect(zoom_blocks, keep_levels)

    zoom_data <- blks_plot %>%
        dplyr::filter(block %in% zoom_blocks)

    p_zoom <- ggplot(zoom_data, aes(x = focal, y = pct_total, fill = block)) +
        geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
        geom_text(aes(label = sprintf("%.2f%%", pct_total * 100)),
                  position = position_dodge(width = 0.75),
                  vjust = -0.4, size = 2.5) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                           expand = expansion(mult = c(0, 0.15))) +
        scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
        labs(
            x = "Focal cell type",
            y = NULL,
            fill = "Block",
            title = paste(zoom_blocks, collapse = " & ")
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 13),
            legend.position = "none"
        )

    p_full <- p_full + theme(legend.position = "bottom")

    p <- p_full + p_zoom +
        patchwork::plot_layout(widths = c(1, 1))

    p
}

plot_stacked_bar_focal_no_resp <- function(
        blks_all,
        block_levels = c(
            "Cell type",
            "Spatial cell state",
            "Spillover",
            "Residuals"             # kept for completeness; excluded from stack
        ),
        title = "",
        show_inset = TRUE
) {
    library(ggplot2)

    # keep blocks (exclude Residuals from the stack)
    keep_levels <- setdiff(block_levels, "Residuals")

    # residual share per focal (useful for annotations if needed)
    resid_by_focal <- blks_all %>%
        dplyr::filter(.data$block == "Residuals") %>%
        dplyr::mutate(focal = stringr::str_replace(.data$focal, "_", " ")) %>%
        dplyr::transmute(focal, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

    # data for stacked bars (no renormalisation)
    blks_plot <- blks_all %>%
        dplyr::filter(.data$block != "Residuals") %>%
        dplyr::mutate(
            focal = stringr::str_replace(.data$focal, "_", " "),
            block = factor(.data$block, levels = keep_levels)
        ) %>%
        tidyr::complete(focal, block, fill = list(pct_total = 0)) %>%
        dplyr::left_join(resid_by_focal, by = "focal")

    # NPG colors
    cols <- c(
        "Cell type"          = "#003049",
        "Spatial cell state" = "#780000",
        "Spillover"          = "#E76F51"
    )
    cols <- cols[keep_levels]

    # Shared focal ordering: by Spatial cell state descending
    focal_order <- blks_plot %>%
        dplyr::filter(block == "Spatial cell state") %>%
        dplyr::arrange(dplyr::desc(pct_total)) %>%
        dplyr::pull(focal)
    blks_plot$focal <- factor(blks_plot$focal, levels = focal_order)

    # --- Left panel: full stacked bar ---
    p_full <- ggplot(blks_plot, aes(x = focal, y = pct_total, fill = block)) +
        geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
        labs(
            x = "Focal cell type",
            y = "Proportion of total variance",
            fill = "Block",
            title = title
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )

    # --- Right panel: zoomed grouped bar (Spatial cell state + Spillover only) ---
    zoom_cols <- c("Spatial cell state" = cols[["Spatial cell state"]],
                   "Spillover"          = cols[["Spillover"]])

    zoom_data <- blks_plot %>%
        dplyr::filter(block %in% c("Spatial cell state", "Spillover"))

    p_zoom <- ggplot(zoom_data, aes(x = focal, y = pct_total, fill = block)) +
        geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
        geom_text(aes(label = sprintf("%.2f%%", pct_total * 100)),
                  position = position_dodge(width = 0.75),
                  vjust = -0.4, size = 2.5) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                           expand = expansion(mult = c(0, 0.15))) +
        scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
        labs(
            x = "Focal cell type",
            y = NULL,
            fill = "Block",
            title = "Spatial cell state & Spillover"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 13),
            legend.position = "none"
        )

    p_full <- p_full + theme(legend.position = "bottom")

    p <- p_full + p_zoom +
        patchwork::plot_layout(widths = c(1, 1))

    p
}

plot_stacked_bar_gene <- function(
        blks_all,
        block_levels = c(
            "Cell type",
            "Spatial cell state",
            "Responder spatial state",
            "Responder status",
            "Spillover",
            "Residuals"             # kept for completeness; excluded from stack
        ),
        title = ""
) {
    # keep blocks (exclude Residuals from the stack)
    keep_levels <- setdiff(block_levels, "Residuals")

    # residual share per gene (useful for annotations if needed)
    resid_by_gene <- blks_all %>%
        dplyr::filter(.data$block == "Residuals") %>%
        dplyr::mutate(gene = stringr::str_replace(.data$gene, "_", " ")) %>%
        dplyr::transmute(gene, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

    # data for stacked bars (no renormalisation)
    blks_plot <- blks_all %>%
        dplyr::filter(.data$block != "Residuals") %>%
        dplyr::mutate(
            gene = stringr::str_replace(.data$gene, "_", " "),
            block = factor(.data$block, levels = keep_levels)
        ) %>%
        tidyr::complete(gene, block, fill = list(pct_total = 0)) %>%
        dplyr::left_join(resid_by_gene, by = "gene")

    # NPG colors (Responder spatial state = red)
    cols <- c(
        "Cell type"               = "#003049", # teal
        "Spatial cell state"      = "#780000", # purple-gray
        "Responder spatial state" = "#c1121f", # red
        "Responder status"        = "#669bbc", # blue
        "Spillover"               = "#E5A100"
    )
    cols <- cols[keep_levels]  # enforce requested order

    ggplot2::ggplot(blks_plot, ggplot2::aes(x = gene, y = pct_total, fill = block)) +
        ggplot2::geom_col(width = 0.8, position = ggplot2::position_stack(reverse = TRUE)) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
        ggplot2::labs(
            x = "Gene",
            y = "Proportion of total",
            fill = "Block",
            title = title
        ) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )
}

# ==============================================================================
# make_pair_composite()
# Per-pair composite for the BC analysis:
#   Panel C = top-N genes ranked by MCSD;
#   Panel D = their per-gene single-frame decomposition (% of total variance).
#   mv     - PACE bundle carrying $mcsd_canonical[[focal_neighbour]]$scores
#   sf     - single-frame decomposition (only needed for the gene panel)
#   panels - "both" (default) | "mcsd" | "gene"
# ==============================================================================
make_pair_composite <- function(focus_focal, focus_neighbour, mv, sf = NULL, n_top = 5,
                                panels = c("both", "mcsd", "gene")) {
    panels      <- match.arg(panels)
    focus_key   <- paste(focus_focal, focus_neighbour, sep = "_")
    top_drivers <- head(mv$mcsd_canonical[[focus_key]]$scores, n_top)

    # Panel C: top-N genes ranked by MCSD.
    p_mcsd <- ggplot2::ggplot(top_drivers, ggplot2::aes(MCSD, stats::reorder(gene, MCSD))) +
        ggplot2::geom_col(fill = "#F5A623", width = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", MCSD)), hjust = -0.15, size = 3) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
        ggplot2::labs(x = "MCSD (spatial cell state)", y = NULL,
                      title = sprintf("%s <- %s", focus_focal, focus_neighbour)) +
        ggplot2::theme_classic(base_size = 11)
    if (panels == "mcsd") return(p_mcsd)

    # Panel D: per-gene single-frame decomposition (% of total variance), top genes.
    if (is.null(sf)) stop("`sf` (single-frame decomposition) is required for the gene panel.")
    # Condition cohorts carry an extra "Responder spatial %" block (else 3 blocks).
    has_resp    <- "Responder spatial %" %in% names(sf)
    gene_blocks <- c("Cell type", "Spatial cell state",
                     if (has_resp) "Responder spatial state", "Spillover")
    block_cols  <- c("Cell type" = "#1F3B57", "Spatial cell state" = "#8B1A1A",
                     "Responder spatial state" = "#E07B39", "Spillover" = "#F6C9C9")

    gd <- dplyr::filter(sf, focal == focus_focal, gene %in% top_drivers$gene)
    gene_var <- data.frame(gene = gd$gene,
                           `Cell type`          = gd[["Cell type %"]],
                           `Spatial cell state` = gd[["Spatial %"]],
                           `Spillover`          = gd[["Spillover %"]],
                           check.names = FALSE)
    if (has_resp) gene_var[["Responder spatial state"]] <- gd[["Responder spatial %"]]
    gene_var <- gene_var |>
        tidyr::pivot_longer(-gene, names_to = "Block", values_to = "pct") |>
        dplyr::mutate(Block = factor(Block, levels = gene_blocks),
                      gene  = factor(gene, levels = top_drivers$gene))

    p_gene <- ggplot2::ggplot(gene_var, ggplot2::aes(gene, pct / 100, fill = Block)) +
        ggplot2::geom_col(width = 0.78, position = ggplot2::position_stack(reverse = TRUE)) +
        ggplot2::scale_fill_manual(values = block_cols[gene_blocks]) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::labs(x = "Gene", y = "Proportion of total variance",
                      title = "Per-gene decomposition (single-frame)") +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    if (panels == "gene") return(p_gene)

    patchwork::wrap_plots(p_mcsd, p_gene)
}

plot_stacked_bar_gene_no_resp <- function(
        blks_all,
        block_levels = c(
            "Cell type",
            "Spatial cell state",
            "Spillover",
            "Residuals"             # kept for completeness; excluded from stack
        ),
        title = ""
) {
    # keep blocks (exclude Residuals from the stack)
    keep_levels <- setdiff(block_levels, "Residuals")

    # residual share per gene (useful for annotations if needed)
    resid_by_gene <- blks_all %>%
        dplyr::filter(.data$block == "Residuals") %>%
        dplyr::mutate(gene = stringr::str_replace(.data$gene, "_", " ")) %>%
        dplyr::transmute(gene, resid_pct = .data$pct_total, explained = 1 - .data$pct_total)

    # data for stacked bars (no renormalisation)
    blks_plot <- blks_all %>%
        dplyr::filter(.data$block != "Residuals") %>%
        dplyr::mutate(
            gene = stringr::str_replace(.data$gene, "_", " "),
            block = factor(.data$block, levels = keep_levels)
        ) %>%
        tidyr::complete(gene, block, fill = list(pct_total = 0)) %>%
        dplyr::left_join(resid_by_gene, by = "gene")

    # NPG colors (Responder spatial state = red)
    cols <- c(
        "Cell type"               = "#003049", # teal
        "Spatial cell state"      = "#780000", # purple-gray
        "Responder spatial state" = "#c1121f", # red
        "Responder status"        = "#669bbc", # blue
        "Spillover"               = "#E5A100"
    )
    cols <- cols[keep_levels]  # enforce requested order

    ggplot2::ggplot(blks_plot, ggplot2::aes(x = gene, y = pct_total, fill = block)) +
        ggplot2::geom_col(width = 0.8, position = ggplot2::position_stack(reverse = TRUE)) +
        ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        ggplot2::scale_fill_manual(values = cols, breaks = keep_levels, limits = keep_levels, drop = FALSE) +
        ggplot2::labs(
            x = "Gene",
            y = "Proportion of total",
            fill = "Block",
            title = title
        ) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )
}

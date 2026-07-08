# ==============================================================================
# compose_inset_figure.R
# Compose a spatial cell type plot with zoomed gene-expression insets
# ==============================================================================

#' Build a zoomed inset from a plotCellState g1 plot
#'
#' Subsets all layer data to the ROI so density and points are recomputed
#' on the zoomed region only, avoiding the washed-out density problem.
#'
#' @param g1_plot A ggplot object (e.g., from plotCellState(...)$g1)
#' @param xlim Numeric vector of length 2: x range for zoom
#' @param ylim Numeric vector of length 2: y range for zoom
#' @param title Optional title above the inset (e.g., gene name)
#' @param border_col Colour of the inset border (should match ROI box)
#' @param base_size Base font size (default: 12)
#' @return A ggplot object for the zoomed inset
#' Build a zoomed inset in the plotCellState g1 style
#'
#' Rebuilds the plot from data rather than trying to mutate ggplot layers
#' (which breaks with ggrastr). Subsets data to the ROI first, then builds
#' density raster + neighbour points + expression-coloured focal points.
#'
#' @param df_from Data frame of focal cells (with x, y, marker columns)
#' @param df_to Data frame of neighbour cells (with x, y columns)
#' @param xlim,ylim Numeric vectors of length 2 defining the ROI
#' @param marker Column name for expression colouring
#' @param x_col,y_col Column names for coordinates (default: "x", "y")
#' @param title Optional title (e.g., gene name)
#' @param border_col Border colour for the inset panel
#' @param focal_size Point size for focal cells (default: 2)
#' @param neighbour_size Point size for neighbour cells (default: 1)
#' @param base_size Base font size (default: 12)
#' @return A ggplot object
make_zoom_inset <- function(df_from, df_to, xlim, ylim,
                            marker,
                            x_col = "x", y_col = "y",
                            title = NULL, border_col = "black",
                            focal_size = 2, neighbour_size = 1,
                            neighbour_colour = "#4E79A7",
                            show_density = FALSE,
                            show_legend = TRUE,
                            colour_limits = NULL,
                            rasterise = TRUE, raster_dpi = 300,
                            base_size = 12) {
    # Subset both data frames to ROI
    roi_from <- df_from[df_from[[x_col]] >= xlim[1] & df_from[[x_col]] <= xlim[2] &
                        df_from[[y_col]] >= ylim[1] & df_from[[y_col]] <= ylim[2], ]
    roi_to   <- df_to[df_to[[x_col]] >= xlim[1] & df_to[[x_col]] <= xlim[2] &
                      df_to[[y_col]] >= ylim[1] & df_to[[y_col]] <= ylim[2], ]

    p <- ggplot2::ggplot()

    # Layer 1: Neighbour density raster (optional)
    if (show_density && nrow(roi_to) > 10) {
        p <- p + ggplot2::stat_density_2d(
            data = roi_to,
            ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]],
                         fill = ggplot2::after_stat(density)),
            geom = "raster", contour = FALSE
        ) +
            ggplot2::scale_fill_distiller(palette = "Blues", direction = 1,
                                          guide = "none")
    }

    # Layer 2: Neighbour points (uses palette colour to match main plot)
    if (nrow(roi_to) > 0) {
        nb_layer <- ggplot2::geom_point(
            data = roi_to,
            ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
            size = neighbour_size / 2, colour = neighbour_colour
        )
        p <- p + if (rasterise) ggrastr::rasterise(nb_layer, dpi = raster_dpi) else nb_layer
    }

    # Layer 3: Focal cells coloured by expression
    if (nrow(roi_from) > 0) {
        legend_guide <- if (show_legend) ggplot2::guide_colourbar() else "none"
        focal_layer <- ggplot2::geom_point(
            data = roi_from,
            ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]],
                         colour = .data[[marker]]),
            size = focal_size, shape = 19
        )
        p <- p + (if (rasterise) ggrastr::rasterise(focal_layer, dpi = raster_dpi) else focal_layer) +
            ggplot2::scale_colour_gradientn(
                colours = c("black", "darkred", "red", "orange", "yellow"),
                name = marker,
                limits = colour_limits,
                guide = legend_guide
            )
    }

    legend_pos <- if (show_legend) "bottom" else "none"
    p + ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
        ggplot2::labs(title = title, x = NULL, y = NULL) +
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            legend.position = legend_pos,
            legend.key.width = grid::unit(0.8, "cm"),
            legend.key.height = grid::unit(0.25, "cm"),
            legend.title = ggplot2::element_text(size = base_size - 3),
            legend.text = ggplot2::element_text(size = base_size - 4),
            plot.title = ggplot2::element_text(size = base_size - 2,
                                               hjust = 0.5, face = "bold"),
            panel.border = ggplot2::element_rect(colour = border_col,
                                                  fill = NA, linewidth = 1.2),
            axis.title = ggplot2::element_blank(),
            axis.text  = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank()
        )
}


#' Compose a spatial cell type plot with gene-expression zoom insets
#'
#' Main plot shows all cells coloured by cell type. Insets show zoomed ROIs
#' rebuilt from data in plotCellState g1 style (density raster +
#' expression-coloured focal points).
#' ROI rectangles on the main plot match the inset border colours.
#'
#' @param df Data frame with x, y, celltype columns and gene expression columns
#' @param genes Character vector of gene/marker names (used as column names
#'   in df and as inset titles)
#' @param focal Focal cell type (e.g., "Myoepithelial")
#' @param neighbour Neighbour cell type for density (e.g., "Tumour")
#' @param xlim1 Numeric vector of length 2: x range for ROI 1
#' @param ylim1 Numeric vector of length 2: y range for ROI 1
#' @param xlim2 Numeric vector of length 2: x range for ROI 2
#' @param ylim2 Numeric vector of length 2: y range for ROI 2
#' @param ct_cols Named vector of cell type colours. If NULL, uses Tableau 10.
#' @param roi_col1 Colour for ROI 1 box and inset borders (default: "steelblue")
#' @param roi_col2 Colour for ROI 2 box and inset borders (default: "darkorange")
#' @param roi_linewidth Linewidth for ROI rectangles (default: 0.9)
#' @param x_col Column name for x coordinates (default: "x")
#' @param y_col Column name for y coordinates (default: "y")
#' @param celltype_col Column name for cell type (default: "celltype")
#' @param point_size Point size for main plot (default: 0.8)
#' @param point_alpha Point alpha for main plot (default: 0.7)
#' @param focal_size Point size for focal cells in insets (default: 2)
#' @param neighbour_size Point size for neighbour cells in insets (default: 1)
#' @param width_ratio Ratio of main plot to inset grid width (default: c(3, 2))
#' @param base_size Base font size (default: 12)
#' @return A patchwork plot object
compose_inset_figure <- function(
        df,
        genes,
        focal,
        neighbour,
        xlim1, ylim1,
        xlim2, ylim2,
        ct_cols = NULL,
        roi_col1 = "black",       # deliberately OUTSIDE the cell-type palette
        roi_col2 = "#E6007E",     # magenta: distinct from Stromal-orange / Tumour-blue etc.
        roi_linewidth = 0.9,
        x_col = "x",
        y_col = "y",
        celltype_col = "celltype",
        point_size = 0.8,
        point_alpha = 0.7,
        focal_size = 2,
        neighbour_size = 1,
        show_density = FALSE,
        draw_arrows = TRUE,
        arrow_col1 = NULL,
        arrow_col2 = NULL,
        arrow_linewidth = 0.8,
        line_y_nudge = -0.03,
        rasterise = TRUE,
        raster_dpi = 300,
        subsample_main = 25000,
        width_ratio = c(3, 2),
        base_size = 12
) {
    library(ggplot2)
    library(patchwork)

    stopifnot(is.character(genes), length(genes) >= 1)

    # Default Tableau 10 palette
    if (is.null(ct_cols)) {
        ct_cols <- c(
            Tumour         = "#4E79A7",
            Myoepithelial  = "#E15759",
            Stromal        = "#F28E2B",
            Macrophage     = "#59A14F",
            T_Cell         = "#76B7B2",
            B_Cell         = "#EDC948",
            Endothelial    = "#AF7AA1",
            Dendritic_Cell = "#FF9DA7",
            Mast           = "#9C755F",
            Fibroblast     = "#BAB0AC"
        )
    }

    # Pre-filter focal and neighbour data frames
    df_from <- df[df[[celltype_col]] == focal, ]
    df_to   <- df[df[[celltype_col]] == neighbour, ]

    # Get neighbour colour from palette so insets match main plot
    nb_colour <- if (neighbour %in% names(ct_cols)) ct_cols[[neighbour]] else "#4E79A7"

    # --- Subsample for main plot (insets use full data) ----------------------
    if (!is.null(subsample_main) && nrow(df) > subsample_main) {
        set.seed(42)
        df_main <- df[sample(nrow(df), subsample_main), ]
    } else {
        df_main <- df
    }

    # --- Main plot: cell type spatial map + ROI boxes -------------------------
    # Clean up cell type labels for the legend (replace _ with space)
    ct_labels <- setNames(gsub("_", " ", names(ct_cols)), names(ct_cols))

    main_points <- if (rasterise) {
        ggrastr::rasterise(geom_point(size = point_size, alpha = point_alpha),
                           dpi = raster_dpi)
    } else {
        geom_point(size = point_size, alpha = point_alpha)
    }

    p_main <- ggplot(df_main, aes(x = .data[[x_col]], y = .data[[y_col]],
                              colour = .data[[celltype_col]])) +
        main_points +
        scale_colour_manual(values = ct_cols, labels = ct_labels,
                            name = "Cell type",
                            guide = guide_legend(override.aes = list(size = 5, alpha = 1))) +
        annotate("rect",
                 xmin = xlim1[1], xmax = xlim1[2],
                 ymin = ylim1[1], ymax = ylim1[2],
                 fill = NA, colour = roi_col1, linewidth = roi_linewidth) +
        annotate("rect",
                 xmin = xlim2[1], xmax = xlim2[2],
                 ymin = ylim2[1], ymax = ylim2[2],
                 fill = NA, colour = roi_col2, linewidth = roi_linewidth) +
        coord_equal() +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = base_size) +
        theme(
            plot.title    = element_text(hjust = 0.5, face = "bold",
                                         size = rel(1.3)),
            legend.title  = element_text(size = rel(1.05)),
            legend.text   = element_text(size = rel(1.0)),
            legend.position = "left"
        )

    # --- Build zoom insets: n_genes columns x 2 rows --------------------------
    # Shared colour limits across all genes and both ROIs
    shared_clim <- c(0, 10)
    gene_limits <- setNames(lapply(genes, function(g) shared_clim), genes)

    zoom_panels <- lapply(genes, function(g) {
        clim <- gene_limits[[g]]
        list(
            # Top row: no legend (shared colourbar goes on bottom row)
            r1 = make_zoom_inset(df_from, df_to, xlim1, ylim1,
                                 marker = g, x_col = x_col, y_col = y_col,
                                 title = g, border_col = roi_col1,
                                 focal_size = focal_size,
                                 neighbour_size = neighbour_size,
                                 neighbour_colour = nb_colour,
                                 show_density = show_density,
                                 show_legend = FALSE,
                                 colour_limits = clim,
                                 rasterise = rasterise, raster_dpi = raster_dpi,
                                 base_size = base_size),
            # Bottom row: no legend either (shared colourbar goes in middle row)
            r2 = make_zoom_inset(df_from, df_to, xlim2, ylim2,
                                 marker = g, x_col = x_col, y_col = y_col,
                                 title = NULL, border_col = roi_col2,
                                 focal_size = focal_size,
                                 neighbour_size = neighbour_size,
                                 neighbour_colour = nb_colour,
                                 show_density = show_density,
                                 show_legend = FALSE,
                                 colour_limits = clim,
                                 rasterise = rasterise, raster_dpi = raster_dpi,
                                 base_size = base_size)
        )
    })

    # Build single shared colourbar for the middle row (spans all gene columns)
    shared_clim <- c(0, 10)
    dummy <- data.frame(x = 1:2, y = 1:2, val = shared_clim)
    cbar_plot <- ggplot2::ggplot(dummy, ggplot2::aes(x = x, y = y, colour = val)) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::scale_colour_gradientn(
            colours = c("black", "darkred", "red", "orange", "yellow"),
            name = "Expression", limits = shared_clim,
            guide = ggplot2::guide_colourbar(
                direction = "horizontal",
                title.position = "left",
                barwidth = grid::unit(8, "cm"),
                barheight = grid::unit(0.35, "cm")
            )
        ) +
        ggplot2::theme_void(base_size = base_size) +
        ggplot2::theme(
            legend.position = "top",
            legend.title = ggplot2::element_text(size = base_size - 2),
            legend.text = ggplot2::element_text(size = base_size - 3)
        )

    # Compose inset grid: top row / colourbar row / bottom row
    top_row  <- Reduce(`|`, lapply(zoom_panels, `[[`, "r1"))
    bot_row  <- Reduce(`|`, lapply(zoom_panels, `[[`, "r2"))
    inset_grid <- top_row / cbar_plot / bot_row +
        patchwork::plot_layout(heights = c(10, 1.5, 10))

    # --- Compose final figure -------------------------------------------------
    combined <- p_main + inset_grid + plot_layout(widths = width_ratio)

    if (!draw_arrows) return(combined)

    # --- Connector lines from each ROI box to its inset -----------------------
    # Drawn as part of the RETURNED object (cowplot canvas), not as a grid
    # side-effect: a grid.segments() overlay is lost when the patchwork is
    # re-rendered by knitr / on assignment. Positions use estimated panel bounds
    # in npc [0, 1].
    wr <- width_ratio / sum(width_ratio)
    mp_left   <- 0.12
    mp_right  <- mp_left + wr[1] * 0.83
    mp_bottom <- 0.06
    mp_top    <- 0.95

    x_range <- range(df[[x_col]], na.rm = TRUE)
    y_range <- range(df[[y_col]], na.rm = TRUE)
    data_to_npc <- function(dx, dy) {
        c(x = mp_left   + (dx - x_range[1]) / diff(x_range) * (mp_right - mp_left),
          y = mp_bottom + (dy - y_range[1]) / diff(y_range) * (mp_top - mp_bottom))
    }

    # Inset block position estimates (heights = c(10, 1, 10): each inset row = 10/21).
    inset_left   <- mp_right + wr[2] * 0.05
    inset_h      <- mp_top - mp_bottom
    top_center_y <- mp_top    - (10 / 21 * inset_h) / 2
    bot_center_y <- mp_bottom + (10 / 21 * inset_h) / 2

    # From the right edge of each ROI box to the left edge of its inset row.
    start1 <- data_to_npc(xlim1[2], mean(ylim1)); start1["y"] <- start1["y"] + line_y_nudge
    start2 <- data_to_npc(xlim2[2], mean(ylim2)); start2["y"] <- start2["y"] + line_y_nudge
    end1   <- c(x = inset_left, y = top_center_y)
    end2   <- c(x = inset_left, y = bot_center_y)

    # Each connector defaults to its ROI/inset colour so the pairing is obvious.
    acol1 <- if (is.null(arrow_col1)) roi_col1 else arrow_col1
    acol2 <- if (is.null(arrow_col2)) roi_col2 else arrow_col2

    cowplot::ggdraw(patchwork::patchworkGrob(combined)) +
        cowplot::draw_line(x = c(start1[["x"]], end1[["x"]]),
                           y = c(start1[["y"]], end1[["y"]]),
                           colour = acol1, linewidth = arrow_linewidth, linetype = "dotted") +
        cowplot::draw_line(x = c(start2[["x"]], end2[["x"]]),
                           y = c(start2[["y"]], end2[["y"]]),
                           colour = acol2, linewidth = arrow_linewidth, linetype = "dotted")
}

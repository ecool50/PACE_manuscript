# =============================================================================
# pace_response_figures.R
# Reusable figures for a focal gene vs neighbour-density, optionally split by a
# grouping (e.g. ICI response). General over gene / focal type / neighbour type /
# group -- the mel SPP1-macrophage<-Tumour figures are one instantiation.
# Auto-sourced by scripts/zzz.R (lives in scripts/helpers/). Assumes ggplot2,
# dplyr and patchwork are attached at call time (as in the analysis notebooks).
#
# Functions:
#   pace_response_slopes()        per-group PACE shrunken slope, read from the fit
#   plot_binned_slope_by_group()  POPULATION view: binned means +/- SE + PACE slope
#   spatial_density_expr_map()    one tissue map: neighbour-density KDE + focal cells
#   density_expr_scatter()        one within-image scatter (expr vs density, lm + stats)
#   response_panel_grid()         compose N maps (top) + N scatters (bottom)
# =============================================================================

# -----------------------------------------------------------------------------
# pace_response_slopes(): the PACE shrunken slope of `gene` for focal<-neighbour,
# per response group, read straight from the fitted model (mv$fit$U) -- NOT an
# OLS / geom_smooth fit. The reference group is the focal::neighbour BLUP; the
# alternate (responder) group adds the focal::<resp_term>:neighbour interaction.
#   fit       : mv$fit (carries $U = [term x gene] shrunken slopes)
#   gene      : gene name (column of fit$U)
#   focal     : focal cell type
#   neighbour : neighbour cell type
#   resp_term : condition term, e.g. "ResponderPD" (NULL -> single slope, both equal)
# Returns a named numeric c(<ref_name>, <alt_name>).
# -----------------------------------------------------------------------------
pace_response_slopes <- function(fit, gene, focal, neighbour, resp_term = NULL,
                                 ref_name = "ref", alt_name = "alt") {
  if (!gene %in% colnames(fit$U)) stop("gene not in fit$U: ", gene)
  u <- fit$U[, gene]; names(u) <- rownames(fit$U)
  base_term <- paste0(focal, "::", neighbour)
  if (!base_term %in% names(u)) stop("slope term not found: ", base_term)
  s_ref <- unname(u[[base_term]])
  s_alt <- s_ref
  if (!is.null(resp_term)) {
    resp_int <- paste0(focal, "::", resp_term, ":", neighbour)
    if (resp_int %in% names(u)) s_alt <- s_ref + unname(u[[resp_int]])
  }
  stats::setNames(c(s_ref, s_alt), c(ref_name, alt_name))
}

# -----------------------------------------------------------------------------
# plot_binned_slope_by_group(): population view of `expr` vs neighbour `density`,
# split by `group`. Solid line + ribbon = binned means +/- SE (the observed data);
# dashed line = the PACE slope (slope_by_group), anchored at each group's data
# centroid and drawn across the panel.
#   density        : neighbour density per focal cell (numeric)
#   expr           : focal-cell expression, e.g. log1p(CP10K) (numeric)
#   group          : grouping per focal cell (factor/character)
#   slope_by_group : named numeric, PACE slope per group level (see pace_response_slopes)
#   group_colours  : named colour per group level
# -----------------------------------------------------------------------------
plot_binned_slope_by_group <- function(density, expr, group, slope_by_group,
                                       group_colours, group_labels = NULL,
                                       n_bins = 10, min_bin_n = 10,
                                       x_range = NULL, y_max = NULL,
                                       title = NULL, subtitle = NULL, legend_title = NULL,
                                       x_lab = "Neighbour density (Gaussian kernel)",
                                       y_lab = "Expression (CP10K, log1p)") {
  d <- data.frame(density = as.numeric(density), expr = as.numeric(expr),
                  group = factor(group))
  d <- d[is.finite(d$density) & is.finite(d$expr), ]
  if (is.null(x_range)) x_range <- c(0, max(d$density))

  # solid: binned means +/- SE over density quantile bins; drop sparse bins
  brks <- unique(stats::quantile(d$density, seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  d$bin <- cut(d$density, breaks = brks, include.lowest = TRUE)
  binned <- d |>
    dplyr::group_by(group, bin) |>
    dplyr::summarise(x = mean(density), mean_e = mean(expr),
                     se = stats::sd(expr) / sqrt(dplyr::n()), n = dplyr::n(),
                     .groups = "drop") |>
    dplyr::filter(n >= min_bin_n)

  # dashed: PACE slope anchored at each group's data centroid
  centro <- d |>
    dplyr::group_by(group) |>
    dplyr::summarise(dx = mean(density), ey = mean(expr), .groups = "drop")
  pace_lines <- do.call(rbind, lapply(seq_len(nrow(centro)), function(i) {
    g <- as.character(centro$group[i])
    data.frame(group = centro$group[i], x = x_range,
               y = centro$ey[i] + slope_by_group[[g]] * (x_range - centro$dx[i]))
  }))

  if (is.null(y_max)) y_max <- max(binned$mean_e + binned$se) * 1.05

  ggplot() +
    geom_ribbon(data = binned,
                aes(x, ymin = mean_e - se, ymax = mean_e + se, fill = group), alpha = 0.18) +
    geom_line(data = binned, aes(x, mean_e, colour = group), linewidth = 0.9) +
    geom_point(data = binned, aes(x, mean_e, colour = group), size = 1.8) +
    geom_line(data = pace_lines, aes(x, y, colour = group),
              linetype = "dashed", linewidth = 1.1) +
    scale_colour_manual(values = group_colours, name = legend_title,
                        labels = if (is.null(group_labels)) ggplot2::waiver() else group_labels) +
    scale_fill_manual(values = group_colours, guide = "none") +
    coord_cartesian(xlim = x_range, ylim = c(0, y_max)) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
}

# -----------------------------------------------------------------------------
# neighbour_density_kde(): smooth 2-D kernel-density estimate of `neighbour_type`
# cell positions for one image, scaled to counts. Returned as a long data.frame
# (x, y, z). Compute it once for several images to share a density colour scale.
# -----------------------------------------------------------------------------
neighbour_density_kde <- function(cell_df, neighbour_type, celltype_col = "celltype",
                                  x_col = "x", y_col = "y",
                                  ngrid = 200, kde_h = 60, padding = 20) {
  nb <- cell_df[as.character(cell_df[[celltype_col]]) == neighbour_type, , drop = FALSE]
  xr <- range(cell_df[[x_col]]) + c(-padding, padding)
  yr <- range(cell_df[[y_col]]) + c(-padding, padding)
  k  <- MASS::kde2d(nb[[x_col]], nb[[y_col]], n = ngrid, h = c(kde_h, kde_h), lims = c(xr, yr))
  kde <- expand.grid(x = k$x, y = k$y)
  kde$z <- as.vector(k$z) * nrow(nb)
  kde
}

# -----------------------------------------------------------------------------
# spatial_density_expr_map(): one tissue map for a single image -- a smooth
# neighbour-type density KDE (white -> navy) with focal-type cells overlaid and
# coloured by `expr_col` on a black-body heat ramp; 100 um scale bar; theme_void.
# Pass shared `expr_limits` / `density_limits` to keep colour scales identical
# across panels, and a precomputed `kde` (from neighbour_density_kde) to avoid
# recomputing it (e.g. when you derived shared density limits from it).
#   cell_df : data.frame of one image's cells (columns: x_col, y_col, celltype_col, expr_col)
# -----------------------------------------------------------------------------
spatial_density_expr_map <- function(cell_df, focal_type, neighbour_type, expr_col,
                                     celltype_col = "celltype", x_col = "x", y_col = "y",
                                     kde = NULL, expr_limits = NULL, density_limits = NULL,
                                     title = NULL, scalebar_len = 100,
                                     expr_name = "expression (log CP10K)",
                                     density_name = "neighbour density",
                                     ngrid = 200, kde_h = 60, padding = 20) {
  fc <- cell_df[as.character(cell_df[[celltype_col]]) == focal_type, , drop = FALSE]
  x  <- cell_df[[x_col]]; y <- cell_df[[y_col]]
  if (is.null(kde))
    kde <- neighbour_density_kde(cell_df, neighbour_type, celltype_col, x_col, y_col,
                                 ngrid, kde_h, padding)

  # black-body heat ramp (0 = black ... high = pale yellow); no purple
  heat_cols <- c("#000000", "#4D0000", "#8B0000", "#E41A1C", "#FF7F00", "#FFD92F", "#FFFFB2")
  emax <- if (is.null(expr_limits)) max(fc[[expr_col]], na.rm = TRUE) else expr_limits[2]
  heat_values <- scales::rescale(c(0, 2, 3.2, 4.2, 5.2, 6.2, emax), from = c(0, emax))

  sb_x0 <- range(x)[1] + 0.06 * diff(range(x))
  sb_y0 <- range(y)[1] + 0.06 * diff(range(y))

  ggplot() +
    geom_raster(data = kde, aes(x, y, fill = z), interpolate = TRUE) +
    geom_point(data = fc, aes(.data[[x_col]], .data[[y_col]], colour = .data[[expr_col]]),
               size = 2.4, alpha = 0.97, shape = 16) +
    scale_fill_gradientn(
      colours = c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
      limits = density_limits, name = density_name,
      guide = guide_colourbar(order = 2, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    scale_colour_gradientn(
      colours = heat_cols, values = heat_values, limits = expr_limits, name = expr_name,
      guide = guide_colourbar(order = 1, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    annotate("segment", x = sb_x0, xend = sb_x0 + scalebar_len, y = sb_y0, yend = sb_y0,
             colour = "black", linewidth = 1.4) +
    annotate("text", x = sb_x0 + scalebar_len / 2, y = sb_y0 + 0.035 * diff(range(y)),
             label = paste0(scalebar_len, " µm"), size = 2.8, colour = "black", vjust = 0) +
    coord_equal(expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.position  = "bottom", legend.box = "horizontal",
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          plot.title   = element_text(size = 11, hjust = 0.5, margin = margin(b = 4)),
          plot.margin  = margin(6, 6, 6, 6))
}

# -----------------------------------------------------------------------------
# density_expr_scatter(): one within-image scatter of `expr` vs neighbour
# `density`, with an lm fit (+95% CI) and an in-panel slope/p/R^2/n box. Pass
# shared `x_limits`/`y_limits` to align panels.
# -----------------------------------------------------------------------------
density_expr_scatter <- function(density, expr, point_colour = "#1f77b4",
                                 title = NULL, x_limits = NULL, y_limits = NULL,
                                 x_lab = "Neighbour density (Gaussian kernel)",
                                 y_lab = "Expression (log CP10K)") {
  d  <- data.frame(density = as.numeric(density), expr = as.numeric(expr))
  d  <- d[is.finite(d$density) & is.finite(d$expr), ]
  f  <- stats::lm(expr ~ density, data = d); co <- summary(f)$coefficients
  ann <- sprintf("slope = %+.3f\np = %s\nR² = %.3f\nn = %d",
                 co[2, 1], format.pval(co[2, 4], digits = 2, eps = 1e-3),
                 summary(f)$r.squared, nrow(d))
  if (is.null(x_limits)) x_limits <- range(d$density)
  if (is.null(y_limits)) y_limits <- range(d$expr)
  ggplot(d, aes(density, expr)) +
    geom_point(colour = point_colour, alpha = 0.5, size = 1.1) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey60", alpha = 0.3) +
    annotate("label", x = x_limits[1], y = y_limits[2], hjust = 0, vjust = 1, label = ann,
             size = 2.9, label.size = 0, fill = scales::alpha("white", 0.7)) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10))
}

# -----------------------------------------------------------------------------
# response_panel_grid(): compose a list of map ggplots (top row) and a list of
# scatter ggplots (bottom row) into one figure, collecting the (shared) map
# legends once at the bottom. `maps` and `scatters` must be equal-length lists.
# -----------------------------------------------------------------------------
response_panel_grid <- function(maps, scatters, map_height = 1.4) {
  top <- Reduce(`+`, maps)
  bot <- Reduce(`+`, scatters)
  (top / bot) +
    patchwork::plot_layout(heights = c(map_height, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
}

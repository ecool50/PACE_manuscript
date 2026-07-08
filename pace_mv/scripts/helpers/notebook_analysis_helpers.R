# =============================================================================
# notebook_analysis_helpers.R
# -----------------------------------------------------------------------------
# Helper functions used by the analysis notebooks (notebooks/analysis_bc.qmd
# and notebooks/analysis_mel.qmd). Auto-sourced by `.load_all("scripts")`
# (see scripts/zzz.R), so the notebooks call these without defining them inline.
#
# Contents:
#   area_fraction()          isotropic edge correction for the E^tech field (both)
#   rank_pair_drivers()      MCSD driver ranking, spatial slope          (BC)
#   rank_pair_drivers_rxs()  MCSD driver ranking, disease-modulated slope (Mel)
#   density_bin_boxplot()    gene expression across neighbour-density deciles (BC fig)
#   top_gene()               top MCSD gene for a pair key                 (BC fig)
#   fp_curve()               cumulative contamination-FP curve            (BC ROC)
#
# NOTE: density_bin_boxplot(), top_gene() and fp_curve() read notebook globals
# (df, Y, K_bio, pair_drivers, tumour_spec) at call time by design — they are
# only invoked from later notebook chunks where those objects already exist.
# =============================================================================

# ---------------------------------------------------------------------------
# area_fraction(): isotropic edge correction for the E^tech spillover field.
# For each cell, returns area(disc(cell, r) ∩ image_rectangle) / (π r²) — the
# fraction of its radius-r neighbourhood disc that lies inside the imaged
# region. Computed by angular quadrature over n_angles directions: in each
# direction the disc is clipped at min(r, distance to nearest wall), and the
# clipped-disc area is summed as (1/2) Σ r_eff(theta)² dtheta = π Σ r_eff² / N.
# Dependency-free re-implementation of the thesis-era spatstat borderEdge()
# (edge = "isotropic"); same quantity used in compute_abundance.R.
# ---------------------------------------------------------------------------
area_fraction <- function(coords, radius, x_min, x_max, y_min, y_max,
                          n_angles = 1000L) {
  theta  <- seq(0, 2 * pi, length.out = n_angles + 1L)[-1L]
  cos_th <- cos(theta)
  sin_th <- sin(theta)

  n  <- nrow(coords)
  af <- numeric(n)
  for (i in seq_len(n)) {
    x0 <- coords[i, 1]
    y0 <- coords[i, 2]

    # Ray length from the cell to each of the four walls, per direction.
    dist_right  <- ifelse(cos_th > 0, (x_max - x0) / cos_th, Inf)
    dist_left   <- ifelse(cos_th < 0, (x_min - x0) / cos_th, Inf)
    dist_top    <- ifelse(sin_th > 0, (y_max - y0) / sin_th, Inf)
    dist_bottom <- ifelse(sin_th < 0, (y_min - y0) / sin_th, Inf)

    # Clip the disc at the nearest wall in each direction.
    r_eff   <- pmin(radius, pmin(dist_right, dist_left, dist_top, dist_bottom))
    af[i]   <- (pi * sum(r_eff^2) / n_angles) / (pi * radius^2)
  }
  af
}

# ---------------------------------------------------------------------------
# rank_pair_drivers(): MCSD-rank the driver genes for one focal<-neighbour pair.
#
# For a (focal, neighbour) pair, joins the mash-shrunken spatial slopes to the
# per-gene focal specificity, scores each gene by
#     MCSD = shrunken_slope² * specificity² * focal_mean,
# keeps genes with lfsr < 0.05, and ranks them. Also returns per-gene spatial
# vs residual variance shares (R2_S).
#
# Returns list(scores = <tibble>, status = "significant" | "honestly null" |
# "dropped").  "dropped" = the pair's kernel column was removed by the
# drop-sparse-K step (no slope to rank).
# ---------------------------------------------------------------------------
rank_pair_drivers <- function(focal_type, neighbour_type,
                              shrunken_slopes, gene_decomp, fit,
                              cells_by_celltype, alpha_by_gene, gene_names) {

  # Mash-shrunken spatial slopes for this pair (BC: term == neighbour name).
  pair_slopes <- shrunken_slopes |>
    filter(focal == focal_type,
           neighbour == neighbour_type,
           term == neighbour_type) |>
    distinct(gene, .keep_all = TRUE)

  # Per-gene focal specificity + mean expression in the focal cell type.
  pair_spec <- gene_decomp |>
    filter(focal == focal_type) |>
    select(gene, spec, focal_mean) |>
    distinct(gene, .keep_all = TRUE)

  z_col <- paste0(focal_type, "::", neighbour_type)   # random-effect design column
  pair_dropped <- nrow(pair_slopes) == 0 ||
                  nrow(pair_spec) == 0 ||
                  !z_col %in% colnames(fit$re_meta$Z)
  if (pair_dropped) {
    return(list(scores = pair_slopes[0, ], status = "dropped"))
  }

  # Within-focal variance of the neighbour kernel: the multiplier that turns a
  # slope b into a spatial variance share b² * var(K).
  focal_cells   <- cells_by_celltype[[focal_type]]
  neighbour_var <- var(as.numeric(fit$re_meta$Z[focal_cells, z_col]), na.rm = TRUE)

  # Mean fitted expression per gene in the focal cells (for the residual term).
  gene_mean <- Matrix::colMeans(fit$mu[focal_cells, , drop = FALSE])
  names(gene_mean) <- gene_names

  scores <- pair_slopes |>
    inner_join(pair_spec, by = "gene") |>
    mutate(
      MCSD    = estimate_shrunk^2 * spec^2 * pmax(focal_mean, 0),
      mu_bar  = gene_mean[gene],
      alpha   = alpha_by_gene[gene],
      V_S     = estimate_shrunk^2 * neighbour_var,               # spatial variance
      V_resid = log(1 + (1 + pmax(alpha, 0)) / pmax(mu_bar, 1e-6)),  # NB1 residual
      V_total = V_S + V_resid,
      R2_S    = V_S / pmax(V_total, 1e-12)                        # spatial share
    ) |>
    filter(lfsr < 0.05) |>
    arrange(desc(MCSD)) |>
    mutate(rank = row_number()) |>
    rename(b_clean = estimate_shrunk) |>
    select(rank, gene, MCSD, b_clean, spec, focal_mean, R2_S,
           mu_bar, alpha, V_S, V_resid, V_total, lfsr, sd_shrunk)

  status <- if (nrow(scores) >= 3) "significant" else "honestly null"
  list(scores = scores, status = status)
}

# ---------------------------------------------------------------------------
# rank_pair_drivers_rxs(): MCSD-rank drivers of the DISEASE-MODULATED spatial
# response for one focal<-neighbour pair (Mel V_R×S version).
#
# The slope of interest is the `Responder:<neighbour>` interaction
# (resp_term:neighbour, e.g. "ResponderPD:Tumour"). Genes are scored by
#     MCSD = shrunken_interaction_slope² * specificity² * focal_mean,
# filtered at lfsr < 0.05. Also returns the main-spatial (V_S) and
# disease-modulated (V_R×S) per-gene variance shares.
#
# Returns list(scores = <tibble>, status = "significant" | "honestly null").
# ---------------------------------------------------------------------------
rank_pair_drivers_rxs <- function(focal_type, neighbour_type, resp_term,
                                  shrunken_slopes, gene_decomp, fit, df,
                                  cells_by_celltype, alpha_by_gene, gene_names) {

  interaction_term <- paste0(resp_term, ":", neighbour_type)  # e.g. ResponderPD:Tumour

  # Mash-shrunken disease-modulated slopes for this pair.
  pair_slopes <- shrunken_slopes |>
    filter(focal == focal_type,
           neighbour == neighbour_type,
           term == interaction_term) |>
    distinct(gene, .keep_all = TRUE)

  # Per-gene focal specificity + mean expression.
  pair_spec <- gene_decomp |>
    filter(focal == focal_type) |>
    select(gene, spec, focal_mean) |>
    distinct(gene, .keep_all = TRUE)

  if (nrow(pair_slopes) == 0 || nrow(pair_spec) == 0) {
    return(list(scores = pair_slopes[0, ], status = "dropped"))
  }

  Z <- fit$re_meta$Z
  focal_cells <- cells_by_celltype[[focal_type]]

  # Neighbour kernel and responder indicator within the focal cells.
  neighbour_kernel <- as.numeric(Z[focal_cells, paste0(focal_type, "::", neighbour_type)])
  resp_col <- match(paste0(focal_type, "::", resp_term), colnames(Z))
  responder <- if (is.na(resp_col)) df$.resp_dummy[focal_cells]
               else as.numeric(Z[focal_cells, resp_col])

  neighbour_var    <- var(neighbour_kernel, na.rm = TRUE)              # for V_S
  interaction_var  <- var(responder * neighbour_kernel, na.rm = TRUE)  # for V_R×S

  # Mean fitted expression per gene; raw main-spatial BLUP per gene.
  gene_mean <- Matrix::colMeans(fit$mu[focal_cells, , drop = FALSE])
  names(gene_mean) <- gene_names
  main_blup_row <- match(paste0(focal_type, "::", neighbour_type), rownames(fit$U))
  main_blup <- if (is.na(main_blup_row)) setNames(rep(0, length(gene_names)), gene_names)
               else setNames(as.numeric(fit$U[main_blup_row, ]), gene_names)

  scores <- pair_slopes |>
    inner_join(pair_spec, by = "gene") |>
    mutate(
      MCSD    = estimate_shrunk^2 * spec^2 * pmax(focal_mean, 0),
      mu_bar  = gene_mean[gene],
      alpha   = alpha_by_gene[gene],
      u_raw   = main_blup[gene],
      V_S     = u_raw^2 * neighbour_var,                 # main spatial variance
      V_RxS   = estimate_shrunk^2 * interaction_var,     # disease-modulated variance
      V_resid = log(1 + (1 + pmax(alpha, 0)) / pmax(mu_bar, 1e-6)),
      V_total = V_S + V_RxS + V_resid,
      R2_S    = V_S   / pmax(V_total, 1e-12),
      R2_RxS  = V_RxS / pmax(V_total, 1e-12)
    ) |>
    filter(lfsr < 0.05) |>
    arrange(desc(MCSD)) |>
    mutate(rank = row_number()) |>
    rename(b_clean = estimate_shrunk) |>
    select(rank, gene, MCSD, b_clean, u_raw, spec, focal_mean,
           R2_S, R2_RxS, mu_bar, alpha, V_S, V_RxS, V_resid, V_total, lfsr, sd_shrunk)

  status <- if (nrow(scores) >= 3) "significant" else "honestly null"
  list(scores = scores, status = status)
}

# ---------------------------------------------------------------------------
# density_bin_boxplot(): boxplot of a focal gene's RAW counts across bins of the
# NUMBER of neighbour-type cells within `radius` um of each focal cell. Green
# boxes / white fill, a filled dot = per-bin mean, outliers faint. Reads notebook
# globals df (x, y, celltype) and Y at call time. (BC notebook figure.)
#   radius : neighbourhood radius in um for counting neighbours (default h_bio).
#   breaks : neighbour-count cut points; the open top bin is relabelled with the
#            observed maximum count.
# ---------------------------------------------------------------------------
density_bin_boxplot <- function(focal_type, neighbour_type, gene_name,
                                radius = 30, breaks = c(0, 2, 4, 6, 8, Inf),
                                box_colour = "#4F8B5E") {
  focal_idx <- which(df$celltype == focal_type)
  nbr_idx   <- which(df$celltype == neighbour_type)

  # number of neighbour-type cells within `radius` um of each focal cell
  fr <- dbscan::frNN(x = as.matrix(df[nbr_idx, c("x", "y")]), eps = radius,
                     query = as.matrix(df[focal_idx, c("x", "y")]))
  n_neighbours <- lengths(fr$id)

  # Long format: one row per (focal cell, gene). Supports one or many genes;
  # multiple genes are shown as facets.
  expr_mat  <- Y[focal_idx, gene_name, drop = FALSE]         # cells x genes (RAW counts)
  plot_data <- data.frame(
    n_neighbours = rep(n_neighbours, times = length(gene_name)),
    gene         = factor(rep(gene_name, each = length(focal_idx)), levels = gene_name),
    expr         = as.numeric(expr_mat)
  )
  plot_data$bin <- cut(plot_data$n_neighbours, breaks = breaks, include.lowest = TRUE)
  plot_data <- plot_data[!is.na(plot_data$bin), ]
  # relabel the open top bin "(8,Inf]" with the observed maximum count
  lev <- levels(plot_data$bin)
  lev[length(lev)] <- sprintf("(%g,%d]", breaks[length(breaks) - 1L], max(plot_data$n_neighbours))
  levels(plot_data$bin) <- lev

  p <- ggplot(plot_data, aes(bin, expr)) +
    geom_boxplot(colour = box_colour, fill = "white", linewidth = 0.6,
                 outlier.colour = box_colour, outlier.alpha = 0.25, outlier.size = 1.2) +
    stat_summary(fun = mean, geom = "point", colour = box_colour, size = 3.2) +
    labs(x = paste("Number of", neighbour_type, "neighbours"),
         y = paste(focal_type, "RAW counts")) +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

  # Single gene: keep the original title/axis. Multiple: facet by gene.
  if (length(gene_name) == 1L) {
    p + labs(title = paste(gene_name, "expression"),
             y = paste(focal_type, gene_name, "counts"))
  } else {
    p + facet_wrap(~ gene, scales = "free_y") +
      labs(title = sprintf("%s expression across %s-neighbour density",
                           focal_type, neighbour_type))
  }
}

# top MCSD gene for a pair key (reads notebook global pair_drivers). (BC figure.)
top_gene <- function(pair_key) pair_drivers[[pair_key]]$scores$gene[1]

# ---------------------------------------------------------------------------
# fp_curve(): cumulative contamination-false-positive proportion among the
# top-k MCSD drivers of one focal type. A gene is a likely bleed FP if it is
# not focal-specific (spec < 0.5) AND Tumour-owned (tumour_spec > 0.5). Reads
# notebook globals pair_drivers and tumour_spec at call time. (BC ROC.)
# ---------------------------------------------------------------------------
fp_curve <- function(focal_type) {
  scores <- pair_drivers[[paste0(focal_type, "_Tumour")]]$scores
  if (nrow(scores) == 0) return(NULL)

  ranked <- scores |>
    left_join(tumour_spec, by = "gene") |>
    # FP = not focal-specific AND Tumour-owned (likely bleed)
    mutate(is_fp = (spec < 0.5) & (tumour_spec > 0.5)) |>
    arrange(desc(MCSD))

  tibble(
    focal   = focal_type,
    k       = seq_len(nrow(ranked)),
    fp_prop = cumsum(ranked$is_fp) / seq_len(nrow(ranked))
  )
}

# ---------------------------------------------------------------------------
# plot_tau_diagnostic(fit): two-panel view of the EB prior variances tau.
#   A: heatmap of log10(tau) per (focal cell type x term). The (Intercept)
#      COLUMN is cell-type identity (large); neighbour-slope columns are mostly
#      ~0 with a few flagships. (Cohort-agnostic: BC celltype block = 10 terms x
#      9 focals; Mel = 16 terms x 7 focals incl. Responder interactions.)
#   B: all tau ranked on a log scale, coloured by intercept vs spatial slope.
# tau_blocks[[ct]] is stored as [term x focal]; we transpose to [focal x term]
# so the intercept lands in its own column (this is the indexing the earlier
# ad-hoc plot got wrong, which made the large values fall on a diagonal).
# ---------------------------------------------------------------------------
plot_tau_diagnostic <- function(fit, block = "celltype", tau_floor = 1e-7) {
  blocks <- fit$re_meta$blocks
  bi <- which(vapply(blocks, function(b) b$group_col, character(1)) == block)
  if (!length(bi)) stop("plot_tau_diagnostic: no RE block with group_col = '", block, "'")
  tb <- fit$tau_blocks[[bi]]                       # [term x focal]

  long <- as.data.frame(as.table(t(tb)), stringsAsFactors = FALSE)  # t(tb) = [focal x term]
  names(long) <- c("focal", "term", "tau")
  long$term        <- factor(long$term,  levels = rownames(tb))     # (Intercept) first
  long$focal       <- factor(long$focal, levels = rev(colnames(tb)))# first focal at top
  long$is_intercept <- long$term == "(Intercept)"
  long$log10_tau    <- log10(pmax(long$tau, tau_floor))
  long$rank         <- rank(long$tau, ties.method = "first")
  long$kind         <- ifelse(long$is_intercept, "intercept", "spatial slope")

  n_int <- sum(long$is_intercept)
  n_lt  <- sum(long$tau < 5e-4)
  med   <- stats::median(long$tau)
  pal   <- c(intercept = "#B2182B", `spatial slope` = "#2166AC")

  pA <- ggplot2::ggplot(long, ggplot2::aes(term, focal, fill = log10_tau)) +
    ggplot2::geom_tile(colour = "grey92") +
    ggplot2::geom_text(data = long[long$is_intercept, ],          # label the (Intercept) COLUMN
                       ggplot2::aes(label = sprintf("%.2f", tau)), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                  midpoint = -3, name = expression(log[10] ~ tau)) +
    ggplot2::labs(title = bquote("A. Variance component " * tau * " per (focal x term)"),
                  subtitle = "(Intercept) column = cell-type identity (large); slopes = spatial proximity (mostly ~0, a few flagships)",
                  x = "term  (Intercept + neighbour kernel)", y = "focal cell type") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank())

  pB <- ggplot2::ggplot(long, ggplot2::aes(rank, tau, colour = kind)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_hline(yintercept = med, linetype = 2, colour = "grey40") +
    ggplot2::annotate("text", x = 1, y = med, vjust = -0.6, hjust = 0,
                      label = sprintf("median = %.5f", med), size = 3, colour = "grey30") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::labs(title = bquote("B. All " * .(nrow(long)) * " " * tau * " ranked (log scale)"),
                  subtitle = sprintf("%d / %d below 5e-4; large tail = %d intercepts + flagship pairs",
                                     n_lt, nrow(long), n_int),
                  x = "rank", y = expression(tau ~ "(log scale)")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = c(0.22, 0.85))

  if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pA, pB, widths = c(1.1, 1))
  else list(heatmap = pA, ranked = pB)
}

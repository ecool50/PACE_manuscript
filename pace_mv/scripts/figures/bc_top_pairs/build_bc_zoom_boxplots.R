## BC top-two-pairs revision: ZOOMED tumour-ring maps + per-density-bin BOXPLOTS.
## Replaces the previous scatters with boxplots; adds one zoomed tumour "ring" panel.
## Col 1 = ADH1B / Stromal<-Tumour;  Col 2 = MRC1 / Macrophage<-Tumour.
## Style reused from mel build_K_composite.R / build_L_exemplar.R:
##   tumour-density background = smooth 2D KDE of Tumour positions (MASS::kde2d, white->navy);
##   focal cells overlaid as points coloured by gene CP10K on the black-body heat ramp;
##   theme_void maps; 100um scale bar; cairo_pdf; British spelling; log CP10K.
suppressMessages({
  library(ggplot2); library(patchwork); library(MASS); library(scales)
})
set.seed(1)

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/bc_top_pairs"; cache_dir <- "data/figure_cache"
if (!dir.exists(fig_dir)) fig_dir <- "."   # allow running from inside the dir
dat <- readRDS(file.path(cache_dir, "bc_adh1b_mrc1_data.rds"))

stromal <- dat[dat$cellType == "Stromal", ]
mac     <- dat[dat$cellType == "Macrophage", ]
tum     <- dat[dat$cellType == "Tumour", ]
n_str <- nrow(stromal); n_mac <- nrow(mac); n_tum <- nrow(tum)

## black-body heat ramp (shared style)
heat_cols <- c("#000000","#4D0000","#8B0000","#E41A1C","#FF7F00","#FFD92F","#FFFFB2")
heat_vals <- function(rng) scales::rescale(
  c(0, 0.15, 0.3, 0.45, 0.6, 0.8, 1) * rng[2], from = c(0, rng[2]))
dens_fill <- c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B")

adh_rng <- range(stromal$ADH1B_cp10k)
mrc_rng <- range(mac$MRC1_cp10k)

## ============================================================================
## (2) ZOOMED TUMOUR RING — chosen ROI window
## ============================================================================
## Selected by scanning ~800um windows for a compact, central tumour-density peak
## (duct/DCIS-like mass) with focal Stromal & Macrophage cells spanning low->high
## density. Window is fully interior to the section.
ROI_CX <- 3480; ROI_CY <- 1210; ROI_W <- 800
roi_x  <- c(ROI_CX - ROI_W/2, ROI_CX + ROI_W/2)
roi_y  <- c(ROI_CY - ROI_W/2, ROI_CY + ROI_W/2)
in_roi <- function(z) z$x >= roi_x[1] & z$x <= roi_x[2] & z$y >= roi_y[1] & z$y <= roi_y[2]

tum_roi <- tum[in_roi(tum), ]
str_roi <- stromal[in_roi(stromal), ]
mac_roi <- mac[in_roi(mac), ]

## recompute tumour KDE on cells in the window
make_kde <- function(tcells, lims, ngrid = 180, h = c(55, 55)) {
  k  <- MASS::kde2d(tcells$x, tcells$y, n = ngrid, h = h, lims = lims)
  df <- expand.grid(x = k$x, y = k$y)
  df$z <- as.vector(k$z) * nrow(tcells) * (h[1] * h[2])
  df
}
kde_roi <- make_kde(tum_roi, lims = c(roi_x, roi_y), h = c(75, 75))
roi_dens_rng <- range(kde_roi$z)

scalebar_roi <- function(len = 100) {
  x0 <- roi_x[1] + 0.06 * ROI_W; y0 <- roi_y[1] + 0.06 * ROI_W
  list(
    annotate("segment", x = x0, xend = x0 + len, y = y0, yend = y0,
             colour = "black", linewidth = 1.4),
    annotate("text", x = x0 + len/2, y = y0 + 0.03 * ROI_W,
             label = "100 µm", size = 2.8, colour = "black", vjust = 0))
}

make_zoom <- function(focal, gene_col, rng, gene_title, title) {
  ## draw focal cells low->high gene last so bright ones sit on top
  focal <- focal[order(focal[[gene_col]]), ]
  ggplot() +
    geom_raster(data = kde_roi, aes(x, y, fill = z), interpolate = TRUE) +
    geom_point(data = focal, aes(x, y, colour = .data[[gene_col]]),
               size = 3.2, alpha = 0.98, shape = 16) +
    scale_fill_gradientn(
      colours = dens_fill, limits = roi_dens_rng, name = "Tumour density",
      guide = guide_colourbar(order = 3, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    scale_colour_gradientn(
      colours = heat_cols, values = heat_vals(rng), limits = rng, name = gene_title,
      guide = guide_colourbar(direction = "horizontal",
                              title.position = "top", barwidth = 7, barheight = 0.5)) +
    scalebar_roi() +
    coord_equal(expand = FALSE, xlim = roi_x, ylim = roi_y) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom", legend.box = "horizontal",
      legend.title     = element_text(size = 8), legend.text = element_text(size = 7),
      plot.title       = element_text(size = 11, hjust = 0.5, margin = margin(b = 4)),
      plot.margin      = margin(6, 6, 6, 6))
}
z_adh <- make_zoom(str_roi, "ADH1B_cp10k", adh_rng,
                   "Stromal ADH1B (log CP10K)", "ADH1B — Stromal←Tumour") +
  guides(colour = guide_colourbar(order = 1, direction = "horizontal",
                                  title.position = "top", barwidth = 7, barheight = 0.5))
z_mrc <- make_zoom(mac_roi, "MRC1_cp10k", mrc_rng,
                   "Macrophage MRC1 (log CP10K)", "MRC1 — Macrophage←Tumour") +
  guides(colour = guide_colourbar(order = 2, direction = "horizontal",
                                  title.position = "top", barwidth = 7, barheight = 0.5))

## ============================================================================
## CONTEXT: full-section tumour-density map with DOTTED ROI rectangle
## ============================================================================
padding <- 20; ngrid_full <- 220
fx <- range(dat$x) + c(-padding, padding)
fy <- range(dat$y) + c(-padding, padding)
kde_full <- make_kde(tum, lims = c(fx, fy), ngrid = ngrid_full, h = c(60, 60))
full_dens_rng <- range(kde_full$z)

scalebar_full <- function(len = 100) {
  x0 <- fx[1] + 0.05 * diff(fx); y0 <- fy[1] + 0.05 * diff(fy)
  list(
    annotate("segment", x = x0, xend = x0 + len, y = y0, yend = y0,
             colour = "black", linewidth = 1.4),
    annotate("text", x = x0 + len/2, y = y0 + 0.02 * diff(fy),
             label = "100 µm", size = 2.8, colour = "black", vjust = 0))
}

locator <- ggplot() +
  geom_raster(data = kde_full, aes(x, y, fill = z), interpolate = TRUE) +
  annotate("rect", xmin = roi_x[1], xmax = roi_x[2], ymin = roi_y[1], ymax = roi_y[2],
           fill = NA, colour = "black", linetype = "dotted", linewidth = 0.9) +
  scale_fill_gradientn(
    colours = dens_fill, limits = full_dens_rng, name = "Tumour density",
    guide = guide_colourbar(direction = "horizontal", title.position = "top",
                            barwidth = 8, barheight = 0.5)) +
  scalebar_full() +
  coord_equal(expand = FALSE) +
  labs(title = sprintf("Full BC section — ROI window (%.0f, %.0f), %d µm (dotted)",
                       ROI_CX, ROI_CY, ROI_W)) +
  theme_void(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position  = "bottom",
    legend.title     = element_text(size = 8), legend.text = element_text(size = 7),
    plot.title       = element_text(size = 10, hjust = 0.5, margin = margin(b = 4)),
    plot.margin      = margin(6, 6, 6, 6))

ggsave(file.path(fig_dir, "BC_section_ROI.pdf"), locator,
       width = 8, height = 6.2, device = cairo_pdf)

## ============================================================================
## (1) BOXPLOTS — gene CP10K per tumour-density bin (replaces the scatters)
## ============================================================================
## ~40% of focal cells have tumour_density exactly 0, so pure quantile bins
## collapse. Use FIXED sensible breaks: a dedicated "0" bin (no tumour neighbour)
## plus 6 graded bins across the positive range (max focal density ~28), giving
## 7 populated, readable bins ordered low->high density.
DENS_BREAKS <- c(-Inf, 0, 2, 5, 9, 14, 20, Inf)
DENS_LABELS <- c("0", "0–2", "2–5", "5–9", "9–14", "14–20", ">20")
make_bins <- function(d) {
  d$dens_bin <- cut(d$tumour_density, breaks = DENS_BREAKS,
                    labels = DENS_LABELS, right = TRUE)
  list(d = d, breaks = DENS_BREAKS)
}
fitstat <- function(d, gene_col) {
  f <- lm(d[[gene_col]] ~ d$tumour_density); co <- summary(f)$coefficients
  list(slope = co[2,1], p = co[2,4], n = nrow(d))
}

make_boxplot <- function(d, gene_col, col, ytitle, title) {
  bb <- make_bins(d); db <- bb$d; brk <- bb$breaks
  fit <- fitstat(d, gene_col)
  ## % expressing (>0) and n per bin
  pct <- tapply(db[[gene_col]], db$dens_bin, function(v) 100 * mean(v > 0))
  cnt <- tapply(db[[gene_col]], db$dens_bin, length)
  lev <- levels(db$dens_bin)
  ypos <- min(db[[gene_col]]) - 0.06 * diff(range(db[[gene_col]]))
  lab_df <- data.frame(dens_bin = factor(lev, levels = lev),
                       lab = sprintf("%.0f%%", pct[lev]),
                       y = ypos)
  ## per-bin MEAN trend: the depletion is a frequency/mean shift; with this
  ## zero-inflated gene the box median sits at 0 in most bins, so we overlay the
  ## per-bin mean (diamond + connecting line) to make the downward trend explicit.
  mean_df <- data.frame(dens_bin = factor(lev, levels = lev),
                        m = tapply(db[[gene_col]], db$dens_bin, mean)[lev])
  ann <- sprintf("slope = %+.3f\np = %s\nn = %d",
                 fit$slope, format.pval(fit$p, digits = 2, eps = 1e-300), fit$n)
  ggplot(db, aes(dens_bin, .data[[gene_col]])) +
    geom_boxplot(fill = col, colour = "grey25", alpha = 0.55,
                 outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.4) +
    geom_line(data = mean_df, aes(dens_bin, m, group = 1),
              inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    geom_point(data = mean_df, aes(dens_bin, m),
               inherit.aes = FALSE, colour = "black", fill = "white",
               shape = 23, size = 2.4, stroke = 0.8) +
    geom_text(data = lab_df, aes(dens_bin, y, label = lab),
              inherit.aes = FALSE, size = 2.5, colour = "grey30", vjust = 1) +
    annotate("label", x = 0.6, y = max(db[[gene_col]]), hjust = 0, vjust = 1,
             label = ann, size = 2.9, fill = alpha("white", 0.7)) +
    labs(x = "Tumour-neighbour density bin (low → high)  [% = focal cells expressing]",
         y = ytitle, title = paste0(title, "  (◇ = bin mean)")) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10),
          axis.text.x = element_text(size = 6.5),
          plot.margin = margin(6, 6, 14, 6))
}

bx_adh <- make_boxplot(stromal, "ADH1B_cp10k", "#1f77b4",
                       "Stromal ADH1B (log CP10K)", "ADH1B in Stromal vs tumour density")
bx_mrc <- make_boxplot(mac, "MRC1_cp10k", "#d62728",
                       "Macrophage MRC1 (log CP10K)", "MRC1 in Macrophage vs tumour density")

## ============================================================================
## COMPOSITE: row1 = two zoom maps; row2 = two boxplots; + locator strip on top.
## ============================================================================
zoom_row <- (z_adh + z_mrc) + plot_layout(guides = "collect")
box_row  <- (bx_adh + bx_mrc)

comp <- (locator / zoom_row / box_row) +
  plot_layout(heights = c(0.75, 1.4, 1.0)) &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig_dir, "BC_ADH1B_MRC1_zoom_boxplots.pdf"), comp,
       width = 11, height = 15, device = cairo_pdf)

## ============================================================================
## VALIDATION / REPORT
## ============================================================================
fa <- fitstat(stromal, "ADH1B_cp10k"); fm <- fitstat(mac, "MRC1_cp10k")
nb_str <- table(make_bins(stromal)$d$dens_bin)
nb_mac <- table(make_bins(mac)$d$dens_bin)
cat("\n=== ROI WINDOW ===\n")
cat(sprintf("centre (%.0f, %.0f), size %d µm | x[%.0f,%.0f] y[%.0f,%.0f]\n",
            ROI_CX, ROI_CY, ROI_W, roi_x[1], roi_x[2], roi_y[1], roi_y[2]))
cat(sprintf("focal in window: Stromal=%d  Macrophage=%d  Tumour=%d\n",
            nrow(str_roi), nrow(mac_roi), nrow(tum_roi)))
cat(sprintf("ROI KDE z range [%.2f, %.2f] (%s)\n", roi_dens_rng[1], roi_dens_rng[2],
            ifelse(diff(roi_dens_rng) > 1e-6, "varies", "FLAT")))
cat(sprintf("in-window cor(gene,td): Stromal ADH1B=%.2f  Macrophage MRC1=%.2f\n",
            cor(str_roi$ADH1B_cp10k, str_roi$tumour_density),
            cor(mac_roi$MRC1_cp10k, mac_roi$tumour_density)))
cat("\n=== BIN BREAKS (fixed, tumour density) ===\n")
cat("breaks:", paste(DENS_BREAKS, collapse = ", "), "\n")
cat("labels:", paste(DENS_LABELS, collapse = " | "), "\n")
cat("Stromal n/bin:   ", paste(nb_str, collapse = " "), "\n")
cat("Macrophage n/bin:", paste(nb_mac, collapse = " "), "\n")
cat("\n=== OVERALL FITS (sanity vs expected) ===\n")
cat(sprintf("ADH1B/Stromal:   slope=%+.4f (exp -0.122)  p=%.3g (exp 5.7e-130)  n=%d\n",
            fa$slope, fa$p, fa$n))
cat(sprintf("MRC1/Macrophage: slope=%+.4f (exp -0.198)  p=%.3g (exp 2.2e-187)  n=%d\n",
            fm$slope, fm$p, fm$n))
cat("\n=== SAVED ===\nBC_ADH1B_MRC1_zoom_boxplots.pdf\nBC_section_ROI.pdf\n")

## K (RAW COUNTS): Composite PD-vs-RESP SPP1 demonstration.
## Identical to build_K_composite.R but uses RAW SPP1_count for BOTH the map point
## colour AND the scatter y-axis, and recomputes within-image lm on raw counts.
## Shared SPP1 colour limit = [0, p99] where p99 = 99th pctile of mac SPP1_count
## across the THREE images (32156_17, 32157_18, 32152_14); oob=squish caps outliers.
suppressMessages({
  library(ggplot2); library(patchwork); library(MASS); library(viridisLite)
  library(scales)
})

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/mel_spp1"; cache_dir <- "data/figure_cache"
allcells <- readRDS("/tmp/allcells_cands.rds")
PD <- "32156_17"; RESP <- "32157_18"

## shared p99 colour limit must be computed across THREE images, so source 32152_14 counts.
spp1_p99 <- readRDS(file.path(cache_dir, "spp1_count_p99.rds"))  # written by build_L_exemplar_rawcounts.R
SPP1_LIM <- c(0, spp1_p99)

## ===================== TOP: tissue maps =====================
make_kde <- function(im, ngrid = 200, padding = 20) {
  s   <- allcells[[im]]
  tum <- s[s$celltype_c == "Tumour", ]
  xr  <- range(s$x) + c(-padding, padding)
  yr  <- range(s$y) + c(-padding, padding)
  h <- c(60, 60)
  k <- MASS::kde2d(tum$x, tum$y, n = ngrid, h = h, lims = c(xr, yr))
  df <- expand.grid(x = k$x, y = k$y)
  df$z <- as.vector(k$z) * nrow(tum)
  df
}
kde_pd  <- make_kde(PD)
kde_rsp <- make_kde(RESP)

get_mac <- function(im) {
  s <- allcells[[im]]
  mac <- s[s$celltype_c == "macrophage", ]
  mac$imageID <- im
  mac
}
mac_pd  <- get_mac(PD)
mac_rsp <- get_mac(RESP)

dens_rng <- range(c(kde_pd$z, kde_rsp$z))
## pure black-body heat ramp (0=black ... high=yellow); NO purple. Spread evenly over [0, p99].
heat_cols   <- c("#000000","#4D0000","#8B0000","#E41A1C","#FF7F00","#FFD92F","#FFFFB2")
heat_values <- scales::rescale(seq(0, SPP1_LIM[2], length.out = length(heat_cols)),
                               from = SPP1_LIM)

scalebar_layer <- function(im, len = 100) {
  s  <- allcells[[im]]
  xr <- range(s$x); yr <- range(s$y)
  x0 <- xr[1] + 0.06 * diff(xr)
  y0 <- yr[1] + 0.06 * diff(yr)
  list(
    annotate("segment", x = x0, xend = x0 + len, y = y0, yend = y0,
             colour = "black", linewidth = 1.4),
    annotate("text", x = x0 + len / 2, y = y0 + 0.035 * diff(yr),
             label = "100 µm", size = 2.8, colour = "black", vjust = 0)
  )
}

make_map <- function(im, kde, mac, title) {
  ggplot() +
    geom_raster(data = kde, aes(x, y, fill = z), interpolate = TRUE) +
    geom_point(data = mac, aes(x, y, colour = SPP1_count),
               size = 2.4, alpha = 0.97, shape = 16) +
    scale_fill_gradientn(
      colours = c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
      limits = dens_rng, name = "Tumour density",
      guide = guide_colourbar(order = 2, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    scale_colour_gradientn(
      colours = heat_cols, values = heat_values, limits = SPP1_LIM,
      oob = scales::squish,
      name = "Macrophage SPP1 (raw counts)",
      guide = guide_colourbar(order = 1, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    scalebar_layer(im) +
    coord_equal(expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      plot.title       = element_text(size = 11, hjust = 0.5, margin = margin(b = 4)),
      plot.margin      = margin(6, 6, 6, 6))
}
m_pd  <- make_map(PD,  kde_pd,  mac_pd,  "PD  –  32156_17")
m_rsp <- make_map(RESP, kde_rsp, mac_rsp, "non-PD  –  32157_18")

## ===================== BOTTOM: scatters (raw counts) =====================
fitstat <- function(mac) {
  f <- lm(SPP1_count ~ tumour_density, data = mac); co <- summary(f)$coefficients
  list(slope = co[2,1], p = co[2,4], r2 = summary(f)$r.squared, n = nrow(mac))
}
fp <- fitstat(mac_pd); fr <- fitstat(mac_rsp)

xr <- range(c(mac_pd$tumour_density, mac_rsp$tumour_density))
yr <- range(c(mac_pd$SPP1_count, mac_rsp$SPP1_count))

make_scatter <- function(mac, fit, im, resp_lab, col) {
  ann <- sprintf("slope = %+.3f\np = %s\nR² = %.3f\nn = %d",
                 fit$slope, format.pval(fit$p, digits = 2, eps = 1e-3), fit$r2, fit$n)
  ggplot(mac, aes(tumour_density, SPP1_count)) +
    geom_point(colour = col, alpha = 0.5, size = 1.1) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                fill = "grey60", alpha = 0.3) +
    annotate("label", x = xr[1], y = yr[2], hjust = 0, vjust = 1, label = ann,
             size = 2.9, label.size = 0, fill = alpha("white", 0.7)) +
    coord_cartesian(xlim = xr, ylim = yr) +
    labs(x = "Tumour-neighbour density (Gaussian, h = 30 µm)",
         y = "Macrophage SPP1 (raw counts)",
         title = sprintf("%s  (%s)", im, resp_lab)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10))
}
s_pd  <- make_scatter(mac_pd,  fp, PD,  "PD", "#d62728")
s_rsp <- make_scatter(mac_rsp, fr, RESP, "non-PD",         "#1f77b4")

top <- (m_pd + m_rsp)
bot <- (s_pd + s_rsp)
comp <- (top / bot) +
  plot_layout(heights = c(1.4, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig_dir, "K_PD_vs_RESP_maps_and_scatters_rawcounts.pdf"), comp,
       width = 11, height = 10, device = cairo_pdf)

cat("\n=== VALIDATION (raw counts) ===\n")
cat(sprintf("Shared SPP1 colour limits (raw): [%.3f, %.3f]  (p99 across 3 images)\n", SPP1_LIM[1], SPP1_LIM[2]))
cat(sprintf("PD   max mac SPP1_count = %g  (n over p99: %d)\n", max(mac_pd$SPP1_count), sum(mac_pd$SPP1_count > SPP1_LIM[2])))
cat(sprintf("RESP max mac SPP1_count = %g  (n over p99: %d)\n", max(mac_rsp$SPP1_count), sum(mac_rsp$SPP1_count > SPP1_LIM[2])))
cat(sprintf("\nPD   %s: slope=%+.4f p=%.4g R2=%.3f n=%d\n", PD,  fp$slope, fp$p, fp$r2, fp$n))
cat(sprintf("RESP %s: slope=%+.4f p=%.4g R2=%.3f n=%d\n", RESP, fr$slope, fr$p, fr$r2, fr$n))
cat("\n=== SAVED ===\nK_PD_vs_RESP_maps_and_scatters_rawcounts.pdf\n")

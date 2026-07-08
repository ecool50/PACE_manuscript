## J (RAW COUNTS): Publication-style single-image spatial maps (PD vs RESP).
## Identical to build_J_pub_maps.R but macrophages coloured by RAW SPP1_count.
## Shared SPP1 colour limit = [0, p99] across the THREE images (read from
## spp1_count_p99.rds, written by build_L_exemplar_rawcounts.R). oob=squish caps outliers.
suppressMessages({
  library(ggplot2); library(patchwork); library(MASS); library(scales)
})

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/mel_spp1"; cache_dir <- "data/figure_cache"
allcells <- readRDS("/tmp/allcells_cands.rds")
PD <- "32156_17"; RESP <- "32157_18"

spp1_p99 <- readRDS(file.path(cache_dir, "spp1_count_p99.rds"))
SPP1_LIM <- c(0, spp1_p99)

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
  s[s$celltype_c == "macrophage", ]
}
mac_pd  <- get_mac(PD)
mac_rsp <- get_mac(RESP)

dens_rng <- range(c(kde_pd$z, kde_rsp$z))

heat_cols   <- c("#000000","#4D0000","#8B0000","#E41A1C","#FF7F00","#FFD92F","#FFFFB2")
heat_values <- scales::rescale(seq(0, SPP1_LIM[2], length.out = length(heat_cols)), from = SPP1_LIM)

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

m_pd  <- make_map(PD,  kde_pd,  mac_pd,  "Non-responder (PD)  –  32156_17")
m_rsp <- make_map(RESP, kde_rsp, mac_rsp, "Responder  –  32157_18")

ggsave(file.path(fig_dir, "J_map_PD_32156_17_rawcounts.pdf"),  m_pd,
       width = 5.2, height = 5.6, device = cairo_pdf)
ggsave(file.path(fig_dir, "J_map_RESP_32157_18_rawcounts.pdf"), m_rsp,
       width = 5.2, height = 5.6, device = cairo_pdf)

comb <- (m_pd + m_rsp) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "J_map_PD_vs_RESP_rawcounts.pdf"), comb,
       width = 9.5, height = 5.8, device = cairo_pdf)

cat("\n=== VALIDATION (raw counts) ===\n")
cat(sprintf("Shared SPP1 colour limit (raw): [%.3f, %.3f]\n", SPP1_LIM[1], SPP1_LIM[2]))
cat(sprintf("PD   mac: n=%d  max SPP1_count=%g  over-p99=%d\n", nrow(mac_pd), max(mac_pd$SPP1_count), sum(mac_pd$SPP1_count > SPP1_LIM[2])))
cat(sprintf("RESP mac: n=%d  max SPP1_count=%g  over-p99=%d\n", nrow(mac_rsp), max(mac_rsp$SPP1_count), sum(mac_rsp$SPP1_count > SPP1_LIM[2])))
cat("\n=== SAVED ===\n")
cat("J_map_PD_32156_17_rawcounts.pdf\nJ_map_RESP_32157_18_rawcounts.pdf\nJ_map_PD_vs_RESP_rawcounts.pdf\n")

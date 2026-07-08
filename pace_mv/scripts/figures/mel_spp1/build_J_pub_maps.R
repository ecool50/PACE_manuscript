## J: Publication-style single-image spatial maps (PD vs RESP).
## Background = smooth 2D KDE of TUMOUR cells (white->navy).
## Macrophages = points coloured by SPP1 CP10K (pure black-body heat ramp, no purple).
## theme_void canvas, scale bar (100 um), shared colour/fill limits across panels.
## Data: /tmp/allcells_cands.rds (validated; mac SPP1 matches pd_resp_maps_data.rds).
suppressMessages({
  library(ggplot2); library(patchwork); library(MASS); library(scales)
})

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/mel_spp1"; cache_dir <- "data/figure_cache"
allcells <- readRDS("/tmp/allcells_cands.rds")
PD <- "32156_17"; RESP <- "32157_18"

## ---------- compute smooth tumour-density KDE per image ----------
## Use a shared spatial bandwidth and a fine grid; n=200x200.
make_kde <- function(im, ngrid = 200, padding = 20) {
  s   <- allcells[[im]]
  tum <- s[s$celltype_c == "Tumour", ]
  xr  <- range(s$x) + c(-padding, padding)
  yr  <- range(s$y) + c(-padding, padding)
  ## fixed bandwidth (microns) so the smoothing is comparable across panels
  h <- c(60, 60)
  k <- MASS::kde2d(tum$x, tum$y, n = ngrid, h = h,
                   lims = c(xr, yr))
  df <- expand.grid(x = k$x, y = k$y)
  df$z <- as.vector(k$z)
  ## normalise to a per-image density that reflects local tumour fraction:
  ## scale by number of tumour cells so high z == many tumour cells nearby
  df$z <- df$z * nrow(tum)
  df
}
kde_pd  <- make_kde(PD)
kde_rsp <- make_kde(RESP)

## ---------- macrophage tables ----------
get_mac <- function(im) {
  s <- allcells[[im]]
  s[s$celltype_c == "macrophage", ]
}
mac_pd  <- get_mac(PD)
mac_rsp <- get_mac(RESP)

## ---------- shared scales across BOTH panels ----------
spp1_rng <- range(c(mac_pd$SPP1_cp10k, mac_rsp$SPP1_cp10k))
dens_rng <- range(c(kde_pd$z, kde_rsp$z))

## pure black-body heat ramp (0=black ... high=yellow); NO purple (replaces inferno).
## stops tuned so ~2=dark red, ~4=red/orange, ~6=orange/yellow over [0, spp1_rng[2]].
heat_cols   <- c("#000000","#4D0000","#8B0000","#E41A1C","#FF7F00","#FFD92F","#FFFFB2")
heat_values <- scales::rescale(c(0, 2, 3.2, 4.2, 5.2, 6.2, spp1_rng[2]),
                               from = c(0, spp1_rng[2]))

## ---------- scale-bar helper (100 um) ----------
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

## ---------- map builder ----------
make_map <- function(im, kde, mac, title) {
  ggplot() +
    geom_raster(data = kde, aes(x, y, fill = z), interpolate = TRUE) +
    geom_point(data = mac, aes(x, y, colour = SPP1_cp10k),
               size = 2.4, alpha = 0.97, shape = 16) +
    scale_fill_gradientn(
      colours = c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B"),
      limits = dens_rng, name = "Tumour density",
      guide = guide_colourbar(order = 2, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    scale_colour_gradientn(
      colours = heat_cols, values = heat_values, limits = spp1_rng,
      name = "Macrophage SPP1 (log CP10K)",
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
      plot.title       = element_text(size = 11, hjust = 0.5,
                                      margin = margin(b = 4)),
      plot.margin      = margin(6, 6, 6, 6))
}

m_pd  <- make_map(PD,  kde_pd,  mac_pd,
                  "Non-responder (PD)  –  32156_17")
m_rsp <- make_map(RESP, kde_rsp, mac_rsp,
                  "Responder  –  32157_18")

## ---------- standalone panels ----------
ggsave(file.path(fig_dir, "J_map_PD_32156_17.pdf"),  m_pd,
       width = 5.2, height = 5.6, device = cairo_pdf)
ggsave(file.path(fig_dir, "J_map_RESP_32157_18.pdf"), m_rsp,
       width = 5.2, height = 5.6, device = cairo_pdf)

## ---------- combined side-by-side, shared legends ----------
comb <- (m_pd + m_rsp) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "J_map_PD_vs_RESP.pdf"), comb,
       width = 9.5, height = 5.8, device = cairo_pdf)

## ---------- validation report ----------
rep_one <- function(im, mac) {
  s   <- allcells[[im]]
  nt  <- sum(s$celltype_c == "Tumour")
  nm  <- nrow(mac)
  pos <- mean(mac$SPP1_count > 0) * 100
  kde <- if (im == PD) kde_pd else kde_rsp
  zr  <- range(kde$z)
  cat(sprintf("%s | n_mac=%d  n_tum=%d  %%mac SPP1>0=%.1f%%  KDE z range=[%.2f, %.2f] (smooth: %s)\n",
              im, nm, nt, pos, zr[1], zr[2],
              ifelse(diff(zr) > 1e-6, "yes-varies", "FLAT")))
}
cat("\n=== VALIDATION ===\n")
rep_one(PD,  mac_pd)
rep_one(RESP, mac_rsp)
cat(sprintf("\nShared SPP1 colour limits: [%.3f, %.3f]\n", spp1_rng[1], spp1_rng[2]))
cat(sprintf("Shared tumour-density fill limits: [%.3f, %.3f]\n", dens_rng[1], dens_rng[2]))
cat("\n=== SAVED ===\n")
cat("J_map_PD_32156_17.pdf\nJ_map_RESP_32157_18.pdf\nJ_map_PD_vs_RESP.pdf\n")

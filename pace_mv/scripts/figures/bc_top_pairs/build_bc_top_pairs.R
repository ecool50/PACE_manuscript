## BC top-two-pairs demonstration, in the Mel SPP1 (build_K_composite.R) style.
## Col 1 = ADH1B / Stromal<-Tumour;  Col 2 = MRC1 / Macrophage<-Tumour.
## TOP row  = tumour-density KDE (white->navy) + focal cells coloured by gene CP10K
##            on a black-body heat ramp; theme_void; 100um scale bar.
## BOTTOM   = within-section scatter (gene CP10K vs tumour-neighbour density);
##            lm + 95% CI; in-panel slope/p/R2/n.
## One physical BC section; log CP10K.
suppressMessages({
  library(ggplot2); library(patchwork); library(MASS); library(scales)
})

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
fig_dir   <- "plots/bc_top_pairs"
cache_dir <- "data/figure_cache"
set.seed(1)

ymv <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds")
fit <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
Y <- ymv$Y; df <- ymv$df; Kb <- fit$K_bio

stopifnot(nrow(Y) == nrow(df), nrow(Kb) == nrow(df))

ct    <- as.character(df$cellType)
tdens <- Kb[, "Tumour"]

## CP10K (log) from raw counts
df$ADH1B_cp10k <- log1p(1e4 * Y[, "ADH1B"] / df$nCount)
df$MRC1_cp10k  <- log1p(1e4 * Y[, "MRC1"]  / df$nCount)
df$tumour_density <- tdens
df$cellType_c <- ct

## ---- per-cell data table used for everything ----
dat <- data.frame(
  x = df$x, y = df$y, cellType = ct, nCount = df$nCount,
  tumour_density = tdens,
  ADH1B_cp10k = df$ADH1B_cp10k, MRC1_cp10k = df$MRC1_cp10k
)
saveRDS(dat, file.path(cache_dir, "bc_adh1b_mrc1_data.rds"))

stromal <- dat[dat$cellType == "Stromal", ]
mac     <- dat[dat$cellType == "Macrophage", ]
tum     <- dat[dat$cellType == "Tumour", ]

n_total <- nrow(dat); n_str <- nrow(stromal); n_mac <- nrow(mac); n_tum <- nrow(tum)

## ===================== tumour-density KDE background =====================
padding <- 20; ngrid <- 220
xr <- range(dat$x) + c(-padding, padding)
yr <- range(dat$y) + c(-padding, padding)
k  <- MASS::kde2d(tum$x, tum$y, n = ngrid, h = c(60, 60), lims = c(xr, yr))
kde <- expand.grid(x = k$x, y = k$y)
## scale to expected tumour cells per (h x h) window so the fill is interpretable
kde$z <- as.vector(k$z) * n_tum * (60 * 60)
dens_rng <- range(kde$z)

## ===================== map colour scales (each gene own limits) =========
heat_cols <- c("#000000","#4D0000","#8B0000","#E41A1C","#FF7F00","#FFD92F","#FFFFB2")
heat_vals <- function(rng) scales::rescale(
  c(0, 0.15, 0.3, 0.45, 0.6, 0.8, 1) * rng[2], from = c(0, rng[2]))

adh_rng <- range(stromal$ADH1B_cp10k)
mrc_rng <- range(mac$MRC1_cp10k)

dens_fill <- c("white", "#DEEBF7", "#9ECAE1", "#4292C6", "#08519C", "#08306B")

## subsample for the MAP only (>25k) -- both focal types are <25k here
map_subsample <- function(d, n = 25000) if (nrow(d) > n) d[sample(nrow(d), n), ] else d
stromal_map <- map_subsample(stromal)
mac_map     <- map_subsample(mac)
cat(sprintf("Map cells: Stromal %d (of %d), Macrophage %d (of %d)\n",
            nrow(stromal_map), n_str, nrow(mac_map), n_mac))

scalebar_layer <- function(len = 100) {
  x0 <- xr[1] + 0.06 * diff(xr); y0 <- yr[1] + 0.06 * diff(yr)
  list(
    annotate("segment", x = x0, xend = x0 + len, y = y0, yend = y0,
             colour = "black", linewidth = 1.4),
    annotate("text", x = x0 + len/2, y = y0 + 0.025 * diff(yr),
             label = "100 µm", size = 2.8, colour = "black", vjust = 0))
}

make_map <- function(focal, gene_col, rng, gene_title, title) {
  ggplot() +
    geom_raster(data = kde, aes(x, y, fill = z), interpolate = TRUE) +
    geom_point(data = focal, aes(x, y, colour = .data[[gene_col]]),
               size = 1.5, alpha = 0.95, shape = 16) +
    scale_fill_gradientn(
      colours = dens_fill, limits = dens_rng, name = "Tumour density",
      guide = guide_colourbar(order = 3, direction = "horizontal",
                              title.position = "top", barwidth = 8, barheight = 0.5)) +
    scale_colour_gradientn(
      colours = heat_cols, values = heat_vals(rng), limits = rng,
      name = gene_title,
      guide = guide_colourbar(direction = "horizontal",
                              title.position = "top", barwidth = 7, barheight = 0.5)) +
    scalebar_layer() +
    coord_equal(expand = FALSE) +
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
## distinct colour guide order so two gene bars don't merge
m_adh <- make_map(stromal_map, "ADH1B_cp10k", adh_rng,
                  "Stromal ADH1B (log CP10K)", "ADH1B — Stromal←Tumour") +
  guides(colour = guide_colourbar(order = 1, direction = "horizontal",
                                  title.position = "top", barwidth = 7, barheight = 0.5))
m_mrc <- make_map(mac_map, "MRC1_cp10k", mrc_rng,
                  "Macrophage MRC1 (log CP10K)", "MRC1 — Macrophage←Tumour") +
  guides(colour = guide_colourbar(order = 2, direction = "horizontal",
                                  title.position = "top", barwidth = 7, barheight = 0.5))

## ===================== scatters (all cells, lm + 95% CI) =================
fitstat <- function(d, gene_col) {
  f <- lm(d[[gene_col]] ~ d$tumour_density); co <- summary(f)$coefficients
  list(slope = co[2,1], p = co[2,4], r2 = summary(f)$r.squared, n = nrow(d))
}
fa <- fitstat(stromal, "ADH1B_cp10k")
fm <- fitstat(mac,     "MRC1_cp10k")

make_scatter <- function(d, gene_col, fit, col, ytitle, title) {
  ann <- sprintf("slope = %+.3f\np = %s\nR² = %.3f\nn = %d",
                 fit$slope, format.pval(fit$p, digits = 2, eps = 1e-300),
                 fit$r2, fit$n)
  xr2 <- range(d$tumour_density); yr2 <- range(d[[gene_col]])
  ggplot(d, aes(.data$tumour_density, .data[[gene_col]])) +
    geom_point(colour = col, alpha = 0.4, size = 0.9) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black",
                fill = "grey60", alpha = 0.3) +
    annotate("label", x = xr2[1], y = yr2[2], hjust = 0, vjust = 1, label = ann,
             size = 2.9, label.size = 0, fill = alpha("white", 0.7)) +
    labs(x = "Tumour-neighbour density (Gaussian frNN)", y = ytitle, title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10))
}
s_adh <- make_scatter(stromal, "ADH1B_cp10k", fa, "#1f77b4",
                      "Stromal ADH1B (log CP10K)", "ADH1B in Stromal vs tumour density")
s_mrc <- make_scatter(mac, "MRC1_cp10k", fm, "#d62728",
                      "Macrophage MRC1 (log CP10K)", "MRC1 in Macrophage vs tumour density")

## ===================== compose =====================
top  <- (m_adh + m_mrc)
bot  <- (s_adh + s_mrc)
comp <- (top / bot) +
  plot_layout(heights = c(1.4, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig_dir, "BC_ADH1B_MRC1_maps_scatters.pdf"), comp,
       width = 11, height = 10, device = cairo_pdf)

## ===================== validation report =====================
cat("\n=== VALIDATION ===\n")
cat(sprintf("nrow match: K_bio=%d df=%d Y=%d\n", nrow(Kb), nrow(df), nrow(Y)))
cat(sprintf("n total=%d | Stromal=%d | Macrophage=%d | Tumour=%d\n",
            n_total, n_str, n_mac, n_tum))
cat(sprintf("ADH1B/Stromal:   slope=%+.4f p=%.3g R2=%.4f n=%d  (NEG? %s)\n",
            fa$slope, fa$p, fa$r2, fa$n, fa$slope < 0))
cat(sprintf("MRC1/Macrophage: slope=%+.4f p=%.3g R2=%.4f n=%d  (NEG? %s)\n",
            fm$slope, fm$p, fm$r2, fm$n, fm$slope < 0))
cat(sprintf("ADH1B colour limits: [%.3f, %.3f]\n", adh_rng[1], adh_rng[2]))
cat(sprintf("MRC1  colour limits: [%.3f, %.3f]\n", mrc_rng[1], mrc_rng[2]))
cat(sprintf("Tumour-density fill limits: [%.2f, %.2f] (%s)\n",
            dens_rng[1], dens_rng[2], ifelse(diff(dens_rng) > 1e-6, "varies", "FLAT")))
cat("\n=== SAVED ===\nBC_ADH1B_MRC1_maps_scatters.pdf + bc_adh1b_mrc1_data.rds\n")

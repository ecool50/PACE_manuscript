## L: Strong single-image responder exemplar (32152_14).
## Clean case: within-image macrophage SPP1 ~ tumour-neighbour-density slope = +0.16,
##             p = 1.5e-11, R2 = 0.18, n = 227.
## TOP    = styled tissue map (smooth tumour-density KDE white->navy + macrophages
##          coloured by SPP1 CP10K on the pure black-body heat ramp; theme_void; 100um bar).
##          Same colour/fill limits as the K maps for comparability.
## BOTTOM = within-image scatter (SPP1 CP10K vs tumour-neighbour density; lm + 95% CI;
##          in-panel slope/p/R2/n).
## 32152_14 is NOT in /tmp/allcells_cands.rds, so its data is recomputed from the h5ad
## exactly as build_withSD.R does (px=0.12028; tumour-neighbour density frNN eps=90,
## Gaussian h_bio=30; SPP1 CP10K = log1p(1e4*count/nCount)). KDE: MASS::kde2d (J style).
suppressMessages({
  library(zellkonverter); library(SingleCellExperiment)
  library(dbscan); library(Matrix)
  library(ggplot2); library(patchwork); library(MASS); library(scales)
})
set.seed(1)
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/mel_spp1"; cache_dir <- "data/figure_cache"
px_size <- 0.12028
h_bio   <- 30
radius  <- 90
IM      <- "32152_14"

## ---------- 1. LOAD h5ad (same as build_withSD.R) ----------
sce <- readH5AD("data/simvi_melanoma/Melanoma_5612.h5ad", reader = "R")
assayNames(sce)[assayNames(sce) == "X"] <- "counts"
rownames(sce) <- make.names(rownames(sce))
stopifnot("SPP1" %in% rownames(sce))

cd <- as.data.frame(colData(sce))
cd$imageID <- paste(cd$SPID, cd$fov, sep = "_")

## collapse tumour subtypes
ct <- as.character(cd$celltype)
ct[grepl("^Tumor_", ct)] <- "Tumour"
cd$celltype_c <- ct

## subset to the exemplar image (do this BEFORE density so neighbours are within-image)
keep <- cd$imageID == IM
stopifnot(any(keep))
sce <- sce[, keep]; cd <- cd[keep, ]

## coords (microns)
sp <- reducedDim(sce, "spatial")
cd$x <- sp[, 1] * px_size
cd$y <- sp[, 2] * px_size

## ---------- 2. SPP1 CP10K (over all genes present) ----------
cnt <- assay(sce, "counts")
nCount <- Matrix::colSums(cnt)
spp1_count <- as.numeric(cnt["SPP1", ])
cd$nCount     <- nCount
cd$SPP1_count <- spp1_count
cd$SPP1_cp10k <- log1p(1e4 * spp1_count / nCount)

## ---------- 3. tumour-neighbour density (frNN eps=90, Gaussian h_bio=30) ----------
coords <- as.matrix(cd[, c("x", "y")])
is_tum <- cd$celltype_c == "Tumour"
nn <- frNN(coords, eps = radius)
dens <- numeric(nrow(cd))
for (j in seq_len(nrow(cd))) {
  nb <- nn$id[[j]]
  if (length(nb) == 0) { dens[j] <- 0; next }
  sel <- is_tum[nb]
  if (!any(sel)) { dens[j] <- 0; next }
  d <- nn$dist[[j]][sel]
  dens[j] <- sum(exp(-(d^2) / (h_bio^2)))
}
cd$tumour_density <- dens

## ---------- 4. macrophage table + fit ----------
mac <- cd[cd$celltype_c == "macrophage", ]
f  <- lm(SPP1_cp10k ~ tumour_density, data = mac)
co <- summary(f)$coefficients
fit <- list(slope = co[2, 1], p = co[2, 4],
            r2 = summary(f)$r.squared, n = nrow(mac))

## ---------- VALIDATION GATE (must reproduce +0.16, p~1.5e-11, n=227) ----------
cat("\n================ VALIDATION (32152_14) ================\n")
cat(sprintf("n_mac = %d   (expect 227)\n", fit$n))
cat(sprintf("slope = %+.4f (expect ~+0.16)\n", fit$slope))
cat(sprintf("p     = %.3e  (expect ~1.5e-11)\n", fit$p))
cat(sprintf("R2    = %.3f  (expect ~0.18)\n", fit$r2))
ok <- (abs(fit$slope - 0.16) < 0.03) && (fit$n == 227) && (fit$p < 1e-9)
cat(sprintf("VALIDATION %s\n", ifelse(ok, "PASS", "FAIL")))
if (!ok) stop("32152_14 numbers did not reproduce; aborting before plotting.")

## ---------- 5. shared scales = K-map limits (for comparability) ----------
## Recompute exactly as build_K_composite.R: SPP1 from PD+RESP macrophages,
## tumour density from MASS::kde2d * n_tum on the same two images.
allcells <- readRDS("/tmp/allcells_cands.rds")
PD <- "32156_17"; RESP <- "32157_18"
make_kde_K <- function(im, ngrid = 200, padding = 20) {
  s   <- allcells[[im]]
  tum <- s[s$celltype_c == "Tumour", ]
  xr  <- range(s$x) + c(-padding, padding)
  yr  <- range(s$y) + c(-padding, padding)
  k   <- MASS::kde2d(tum$x, tum$y, n = ngrid, h = c(60, 60), lims = c(xr, yr))
  as.vector(k$z) * nrow(tum)
}
spp1_rng <- range(c(allcells[[PD]][allcells[[PD]]$celltype_c == "macrophage", ]$SPP1_cp10k,
                    allcells[[RESP]][allcells[[RESP]]$celltype_c == "macrophage", ]$SPP1_cp10k))
dens_rng <- range(c(make_kde_K(PD), make_kde_K(RESP)))
cat(sprintf("\nK-shared SPP1 limits: [%.3f, %.3f]   dens limits: [%.4f, %.4f]\n",
            spp1_rng[1], spp1_rng[2], dens_rng[1], dens_rng[2]))

## ---------- 6. tumour-density KDE for the exemplar (J style) ----------
make_kde <- function(s, ngrid = 200, padding = 20) {
  tum <- s[s$celltype_c == "Tumour", ]
  xr  <- range(s$x) + c(-padding, padding)
  yr  <- range(s$y) + c(-padding, padding)
  k   <- MASS::kde2d(tum$x, tum$y, n = ngrid, h = c(60, 60), lims = c(xr, yr))
  df  <- expand.grid(x = k$x, y = k$y)
  df$z <- as.vector(k$z) * nrow(tum)
  df
}
kde <- make_kde(cd)

## ---------- 7. pure black-body heat ramp (matches J/K) ----------
heat_cols   <- c("#000000","#4D0000","#8B0000","#E41A1C","#FF7F00","#FFD92F","#FFFFB2")
heat_values <- scales::rescale(c(0, 2, 3.2, 4.2, 5.2, 6.2, spp1_rng[2]),
                               from = c(0, spp1_rng[2]))

scalebar_layer <- function(s, len = 100) {
  xr <- range(s$x); yr <- range(s$y)
  x0 <- xr[1] + 0.06 * diff(xr)
  y0 <- yr[1] + 0.06 * diff(yr)
  list(
    annotate("segment", x = x0, xend = x0 + len, y = y0, yend = y0,
             colour = "black", linewidth = 1.4),
    annotate("text", x = x0 + len / 2, y = y0 + 0.035 * diff(yr),
             label = "100 µm", size = 2.8, colour = "black", vjust = 0))
}

## ---------- TOP: tissue map ----------
m <- ggplot() +
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
  scalebar_layer(cd) +
  coord_equal(expand = FALSE) +
  labs(title = sprintf("Responder — %s (within-image slope %+.2f, p = %s)",
                       IM, fit$slope, format.pval(fit$p, digits = 2, eps = 1e-3))) +
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

## ---------- BOTTOM: scatter ----------
ann <- sprintf("slope = %+.3f\np = %s\nR² = %.3f\nn = %d",
               fit$slope, format.pval(fit$p, digits = 2, eps = 1e-3), fit$r2, fit$n)
s_plot <- ggplot(mac, aes(tumour_density, SPP1_cp10k)) +
  geom_point(colour = "#1f77b4", alpha = 0.5, size = 1.1) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black",
              fill = "grey60", alpha = 0.3) +
  annotate("label", x = min(mac$tumour_density), y = max(mac$SPP1_cp10k),
           hjust = 0, vjust = 1, label = ann,
           size = 2.9, label.size = 0, fill = alpha("white", 0.7)) +
  labs(x = "Tumour-neighbour density (Gaussian, h = 30 µm)",
       y = "Macrophage SPP1 (log CP10K)",
       title = sprintf("%s  (responder)", IM)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10))

## ---------- compose: map on top, scatter bottom, map legends at bottom ----------
comp <- (m / s_plot) +
  plot_layout(heights = c(1.4, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig_dir, "L_strong_exemplar_32152_14.pdf"), comp,
       width = 6.0, height = 9.2, device = cairo_pdf)

cat("\n=== SAVED ===\nL_strong_exemplar_32152_14.pdf\n")
cat(sprintf("final: slope=%+.4f p=%.3e R2=%.3f n=%d\n",
            fit$slope, fit$p, fit$r2, fit$n))

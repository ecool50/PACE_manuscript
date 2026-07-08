## Paired PD-vs-RESP single-image demo: macrophage SPP1 vs tumour-neighbour density.
## RESP = 32157_18 (focal dense-vs-sparse contrast + positive sig slope)
## PD   = 32156_17 (expected negative sig slope)
## Validated: recomputed slopes match mac_spp1_withSD.rds exactly (cor=1.0).
suppressMessages({
  library(ggplot2); library(patchwork); library(scales); library(akima)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"); fig_dir <- "plots/mel_spp1"; cache_dir <- "data/figure_cache"
allcells <- readRDS("/tmp/allcells_cands.rds")
PD <- "32156_17"; RESP <- "32157_18"

theme_set(theme_bw(base_size=11))

## ----- per-image macrophage data + fit stats -----
get_mac <- function(im) {
  s <- allcells[[im]]
  mac <- s[s$celltype_c=="macrophage", ]
  mac$imageID <- im
  mac
}
mac_pd  <- get_mac(PD)
mac_rsp <- get_mac(RESP)

fitstat <- function(mac) {
  f <- lm(SPP1_cp10k ~ tumour_density, data=mac); co <- summary(f)$coefficients
  list(slope=co[2,1], p=co[2,4], r2=summary(f)$r.squared, n=nrow(mac))
}
fp <- fitstat(mac_pd); fr <- fitstat(mac_rsp)

## save per-image macrophage data used
saveRDS(list(PD=mac_pd, RESP=mac_rsp,
             fits=list(PD=fp, RESP=fr),
             selection=c(PD=PD, RESP=RESP)),
        file.path(cache_dir, "pd_resp_maps_data.rds"))

## ----- shared colour / shading scales (across BOTH panels) -----
spp1_rng <- range(c(mac_pd$SPP1_cp10k, mac_rsp$SPP1_cp10k))
dens_all <- c(allcells[[PD]]$tumour_density, allcells[[RESP]]$tumour_density)
dens_rng <- range(dens_all)

## ----- interpolated continuous density field per image -----
make_field <- function(im) {
  s <- allcells[[im]]
  # interpolate per-cell density onto a grid over the cell bounding box
  ii <- interp(s$x, s$y, s$tumour_density,
               nx=160, ny=160, duplicate="mean")
  df <- expand.grid(x=ii$x, y=ii$y); df$z <- as.vector(ii$z)
  df <- df[!is.na(df$z), ]
  df
}
fld_pd  <- make_field(PD)
fld_rsp <- make_field(RESP)

## ----- tissue map builder -----
make_map <- function(im, mac, fld, fit, resp_lab) {
  s <- allcells[[im]]
  nonmac <- s[s$celltype_c != "macrophage", ]
  ttl <- sprintf("%s  (%s)\nslope = %+.3f, p = %s", im, resp_lab,
                 fit$slope, format.pval(fit$p, digits=2, eps=1e-3))
  ggplot() +
    geom_raster(data=fld, aes(x, y, fill=z), interpolate=TRUE) +
    geom_point(data=nonmac, aes(x, y), colour="grey75", size=0.25, alpha=0.35) +
    geom_point(data=mac, aes(x, y, colour=SPP1_cp10k), size=1.7, alpha=0.95) +
    scale_fill_gradient(low="grey92", high="grey15",
                        limits=dens_rng, name="Tumour-neighbour\ndensity",
                        guide=guide_colourbar(order=1)) +
    scale_colour_viridis_c(option="inferno", limits=spp1_rng,
                           name="Macrophage\nSPP1 (log CP10K)",
                           guide=guide_colourbar(order=2)) +
    coord_equal() +
    labs(x="x (µm)", y="y (µm)", title=ttl) +
    theme(plot.title=element_text(size=9.5, lineheight=1.05),
          legend.key.height=unit(0.7,"cm"), legend.text=element_text(size=7),
          legend.title=element_text(size=8))
}
m_pd  <- make_map(PD,  mac_pd,  fld_pd,  fp, "non-responder, PD")
m_rsp <- make_map(RESP, mac_rsp, fld_rsp, fr, "responder")

## ----- scatters (shared axes) -----
xr <- range(c(mac_pd$tumour_density, mac_rsp$tumour_density))
yr <- range(c(mac_pd$SPP1_cp10k, mac_rsp$SPP1_cp10k))
make_scatter <- function(mac, fit, im, resp_lab, col) {
  ann <- sprintf("slope = %+.3f\np = %s\nR² = %.3f\nn = %d",
                 fit$slope, format.pval(fit$p, digits=2, eps=1e-3), fit$r2, fit$n)
  ggplot(mac, aes(tumour_density, SPP1_cp10k)) +
    geom_point(colour=col, alpha=0.5, size=1.1) +
    geom_smooth(method="lm", formula=y~x, colour="black", fill="grey60", alpha=0.3) +
    annotate("label", x=xr[1], y=yr[2], hjust=0, vjust=1, label=ann,
             size=2.9, label.size=0, fill=alpha("white",0.7)) +
    coord_cartesian(xlim=xr, ylim=yr) +
    labs(x="Tumour-neighbour density (Gaussian, h=30 µm)",
         y="Macrophage SPP1 (log CP10K)",
         title=sprintf("%s  (%s)", im, resp_lab))
}
s_pd  <- make_scatter(mac_pd,  fp, PD,  "non-responder, PD", "#d62728")
s_rsp <- make_scatter(mac_rsp, fr, RESP, "responder",         "#1f77b4")

## ----- compose 2x2 -----
top <- (m_pd + m_rsp) + plot_layout(guides="collect")
bot <- s_pd + s_rsp
comp <- (top / bot) +
  plot_annotation(
    title="Macrophage SPP1 vs tumour proximity: responder vs non-responder (single images)") &
  theme(plot.title=element_text(size=12.5, face="bold"))

ggsave(file.path(fig_dir, "I_PD_vs_RESP_panels.pdf"), comp,
       width=11, height=9.5, device=cairo_pdf)

cat("\n=== SAVED ===\n")
cat("I_PD_vs_RESP_panels.pdf\n")
cat("pd_resp_maps_data.rds\n")
cat(sprintf("\nPD  %s: slope=%+.4f p=%.4f R2=%.3f n=%d\n", PD, fp$slope, fp$p, fp$r2, fp$n))
cat(sprintf("RESP %s: slope=%+.4f p=%.4f R2=%.3f n=%d\n", RESP, fr$slope, fr$p, fr$r2, fr$n))

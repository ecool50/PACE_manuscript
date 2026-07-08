## decomp_gated_eta2_bc.R
## ---------------------------------------------------------------------------
## Principled "neighbour attribution of focal spatial variance" decomposition
## that resolves the sum-to-1-vs-phantom tension.
##
## THE TENSION (proved real): if a column is normalised to sum to 1 OVER
## NEIGHBOURS ONLY, then a focal with near-zero spatial signal (Dendritic, Mast)
## is forced to allocate 100% of nothing to some neighbour -> a "phantom"
## dominant pair (Dendritic<-T_Cell = 98% of a spatial block that is 0.01% of
## total variance). Mass conservation over neighbours alone GUARANTEES this.
## You cannot have both "sum-to-1 over neighbours" and "no phantom".
##
## THE RESOLUTION: keep an honest sum-to-1 decomposition by adding the term
## that the neighbour-only normalisation hides -- the NON-SPATIAL residual.
## We report TWO transparent variance terms, both ANOVA sum-of-squares ratios:
##
##   (1) eta2_spatial,F = SS_spatial,F / SS_total,F          (per-focal GATE)
##         = fraction of focal F's total expression variance that is
##           neighbour-induced spatial variance. Columns of the FULL table
##           (neighbours + "Non-spatial" residual) SUM TO 1 exactly.
##
##   (2) prop_k|F = V_pair,kF / SS_spatial,F                 (within-spatial share)
##         = of the spatial change, the proportion from neighbour k.
##           This is what the collaborator means by "sum to 1".
##
##   HEADLINE metric = gated proportion:  g_F * prop_k|F
##     g_F = eta2_spatial,F / (eta2_spatial,F + tau)   (James-Stein gate)
##     tau = median nonzero eta2_spatial  (focal at the median is half-gated)
##   This is a transparent product of two SS-ratios: a per-focal reliability
##   gate times the within-neighbour proportion. Phantom focals have g_F ~ 0
##   so the WHOLE column dims and no spurious dominant neighbour survives.
##
## Variance terms used (all SS, no resampling/Shapley/RWA/Pillai):
##   V_pair,kF = sum_g  spec_g^2 * b_kg^2 * Var(N_k | exposed)
##   SS_spatial,F = sum_k V_pair,kF
##   SS_total,F   = sum_g spec_g^2 * Total_g
##   spec^2 weight = the discriminator that suppresses precise contamination
##     (KRT7/PTGDS/CCDC80). se^2 dropped: it is ~8% of energy and the
##     inflating genes are PRECISE, not noisy (mashr shrinks only ~16%).
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix); library(dbscan); library(ggplot2); library(tidyr); library(dplyr)
})
say <- function(...) cat(sprintf(...), "\n")

## ---- load canonical fit ----
mv  <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
fit <- mv$fit
U   <- fit$U
TYPES <- sort(unique(as.character(mv$cell_meta$celltype)))
ct    <- as.character(mv$cell_meta$celltype)
gf    <- mv$decomposition$gene_focal_4block
gs    <- mv$gene_set

## ---- raw K_bio neighbour-proximity kernel (global frNN; BC is one section) ----
ymv <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds")
df  <- ymv$df
xy  <- as.matrix(df[, c("x", "y")])
n   <- nrow(df)
h_bio <- 30
radius <- 90
fr <- dbscan::frNN(xy, eps = radius)
row_idx  <- rep.int(seq_len(n), lengths(fr$id))
col_idx  <- unlist(fr$id,   use.names = FALSE)
dist_val <- unlist(fr$dist, use.names = FALSE)
nbr_type <- match(ct[col_idx], TYPES)
keep <- !is.na(nbr_type)
prox <- as.matrix(sparseMatrix(
  i = row_idx[keep], j = nbr_type[keep],
  x = exp(-(dist_val[keep])^2 / h_bio^2),
  dims = c(n, length(TYPES))
))
colnames(prox) <- TYPES
rm(fr, row_idx, col_idx, dist_val, nbr_type, keep); invisible(gc())

## ---- variance components per focal ----
V_pair      <- matrix(0, length(TYPES), length(TYPES),
                      dimnames = list(neighbour = TYPES, focal = TYPES))
SS_total_F  <- setNames(numeric(length(TYPES)), TYPES)

for (focal in TYPES) {
  focal_cells <- which(ct == focal)
  gf_focal    <- gf[gf$focal == focal, ]
  spec        <- gf_focal$spec[match(gs, gf_focal$gene)];  spec[is.na(spec)]   <- 0
  total_var   <- gf_focal$Total[match(gs, gf_focal$gene)]; total_var[is.na(total_var)] <- 0
  SS_total_F[focal] <- sum(spec^2 * total_var)
  for (neighbour in setdiff(TYPES, focal)) {
    slope_row <- paste0(focal, "::", neighbour)
    if (!(slope_row %in% rownames(U))) next
    exposure  <- prox[focal_cells, neighbour]
    exposed   <- exposure > 1e-8
    var_cond  <- if (sum(exposed) > 2) var(exposure[exposed]) else 0
    V_pair[neighbour, focal] <- sum(spec^2 * U[slope_row, ]^2 * var_cond)
  }
}

SS_spatial_F <- colSums(V_pair)

## ---- the two transparent SS-ratios ----
eta2_neighbour <- sweep(V_pair, 2, SS_total_F, "/")          # of TOTAL variance (sums-to-1 w/ residual)
residual_F     <- 1 - colSums(eta2_neighbour)                # Non-spatial term
eta2_spatial_F <- SS_spatial_F / SS_total_F                  # per-focal spatial signal
prop_within    <- sweep(V_pair, 2, ifelse(SS_spatial_F > 0, SS_spatial_F, 1), "/")

## ---- James-Stein reliability gate (signal-to-noise driven shrinkage) ----
tau  <- median(eta2_spatial_F[eta2_spatial_F > 0])
gate <- eta2_spatial_F / (eta2_spatial_F + tau)

## ---- HEADLINE metric: gated within-spatial proportion ----
gated_prop <- sweep(prop_within, 2, gate, "*")

## ---- exact eta2 table that literally sums to 1 (neighbours + residual) ----
eta2_full <- rbind(eta2_neighbour, `Non-spatial` = residual_F)

## ---- report ----
say("=== per-focal gate (tau = %.4f) ===", tau)
for (focal in TYPES)
  say("  %-15s eta2_spatial=%6.3f%%  gate=%.3f  residual=%.2f%%",
      focal, 100 * eta2_spatial_F[focal], gate[focal], 100 * residual_F[focal])

say("\n=== eta2_full table (neighbour shares of TOTAL var + Non-spatial; columns sum to 100%%) ===")
print(round(100 * eta2_full, 3))
say("column sums (check): %s", paste(sprintf("%.0f", colSums(100 * eta2_full)), collapse = " "))

say("\n=== headline gated proportion: real pairs stay prominent ===")
for (p in list(c("Tumour","Stromal"), c("Tumour","Macrophage"),
               c("Tumour","Endothelial"), c("Myoepithelial","Tumour")))
  say("  %-13s<-%-13s  gated = %.0f%%  (prop_within=%.0f%% x gate=%.2f)",
      p[1], p[2], 100 * gated_prop[p[1], p[2]], 100 * prop_within[p[1], p[2]], gate[p[2]])

say("\n=== headline gated proportion: phantoms collapse ===")
for (focal in c("Dendritic_Cell", "Mast", "T_Cell")) {
  top <- names(which.max(gated_prop[, focal]))
  say("  %-15s top neighbour %-12s = %.1f%%  (was %.0f%% under neighbour-only norm)",
      focal, top, 100 * gated_prop[top, focal], 100 * prop_within[top, focal])
}

saveRDS(list(V_pair = V_pair, eta2_neighbour = eta2_neighbour, residual_F = residual_F,
             eta2_spatial_F = eta2_spatial_F, prop_within = prop_within,
             gate = gate, tau = tau, gated_prop = gated_prop, eta2_full = eta2_full),
        "streaming/neighbour_gated_eta2_bc.rds")

## ---- figure: headline gated-proportion heatmap with gate annotated on columns ----
gate_lab <- sprintf("%s\n(gate %.2f)", TYPES, gate[TYPES])
names(gate_lab) <- TYPES

plot_df <- as.data.frame(as.table(gated_prop)) |>
  rename(value = Freq) |>
  mutate(
    focal     = factor(focal, levels = TYPES),
    neighbour = factor(neighbour, levels = rev(TYPES)),
    value     = ifelse(as.character(focal) == as.character(neighbour), NA, value)
  )

p <- ggplot(plot_df, aes(focal, neighbour, fill = value)) +
  geom_tile(colour = "grey92", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(value) | value < 0.03, "", sprintf("%.0f", 100 * value))),
            size = 3) +
  scale_x_discrete(labels = gate_lab) +
  scale_fill_viridis_c(option = "inferno", na.value = "white", limits = c(0, 1),
                       labels = scales::percent, name = "Gated\nspatial share") +
  coord_equal() +
  labs(
    x = "Focal cell type (gate = per-focal spatial signal-to-noise)",
    y = "Neighbour",
    title = "BC 313: gated neighbour attribution of focal spatial variance",
    subtitle = "cell = (within-spatial proportion, sums-to-1) x (eta2-spatial reliability gate); phantom focals dim out"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle = element_text(size = 8))

ggsave("streaming/figures/decomp_gated_eta2_bc.pdf", p, width = 8.4, height = 7, device = cairo_pdf)
say("\nwrote streaming/figures/decomp_gated_eta2_bc.pdf")
say("wrote streaming/neighbour_gated_eta2_bc.rds")

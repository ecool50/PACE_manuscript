## decomp_canonical_bc.R — BC canonical decomposition figures.
##   - Focal stacked bar (no Responder, V_spatial only)
##   - Pairwise V_spatial heatmap
## Uses locked plot functions: plot_stacked_bar_focal_no_resp, make_pair_tile.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(stringr); library(forcats)
  library(ggplot2); library(patchwork); library(scales); library(Matrix)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
PLOT_DIR <- "plots/decomposition"
source("scripts/core/06a-plots-decomp.R")
source("scripts/helpers/pace_pair_variance_pratt.R")

mv <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
if (is.null(colnames(mv$fit$mu))) colnames(mv$fit$mu) <- mv$gene_set

## --- Focal 4-block decomp (BC has no Responder) ---
## BC stores its 4-block percentages in gene_focal_4block (no V_RxS column).
gf <- mv$decomposition$gene_focal_4block

## --- CANONICAL = unweighted POOLED aggregate (Goldstein-Browne-Rasbash 2002):
## pct_m = sum_j V_{m,j} / sum_j Total_j. Each gene contributes its own total
## variance -> clean "fraction of total variance" estimand. The former
## mean-of-ratios (Jensen-biased toward low-variance genes) and the inverse-
## residual-variance-weighted variant (endogenous weight -> inflates non-residual
## blocks) are emitted to the sensitivity CSV only, NOT used for the figure. ---
EPS <- 1e-9
.agg_bc <- function(d, scheme) {
  if (scheme == "mean")
    return(tibble::tibble(
      `Cell type`          = mean(d[["Cell type %"]],     na.rm=TRUE)/100,
      `Spatial cell state` = mean(d[["Spatial state %"]], na.rm=TRUE)/100,
      `Spillover`          = mean(d[["Spillover %"]],     na.rm=TRUE)/100,
      `Residuals`          = mean(d[["Residual %"]],      na.rm=TRUE)/100))
  w   <- if (scheme == "weighted") 1/pmax(d$V_disp, EPS) else rep(1, nrow(d))
  tot <- sum(w * d$Total)
  tibble::tibble(
    `Cell type`          = sum(w*d$celltype_offset_sq)/tot,
    `Spatial cell state` = sum(w*d$V_state)/tot,
    `Spillover`          = sum(w*d$V_spill)/tot,
    `Residuals`          = sum(w*d$V_disp)/tot)
}
blks_all <- gf %>% group_by(focal) %>% group_modify(~ .agg_bc(.x, "pooled")) %>%
  ungroup() %>% pivot_longer(-focal, names_to = "block", values_to = "pct_total")
## sensitivity: weighted-pooled vs unweighted-pooled vs mean-of-ratios
.sens <- bind_rows(
  gf %>% group_by(focal) %>% group_modify(~ .agg_bc(.x,"weighted")) %>% mutate(scheme="weighted_pooled"),
  gf %>% group_by(focal) %>% group_modify(~ .agg_bc(.x,"pooled"))   %>% mutate(scheme="unweighted_pooled"),
  gf %>% group_by(focal) %>% group_modify(~ .agg_bc(.x,"mean"))     %>% mutate(scheme="mean_of_ratios")) %>% ungroup()
dir.create(PLOT_DIR, recursive=TRUE, showWarnings=FALSE)
utils::write.csv(.sens, file.path(PLOT_DIR, "decomp_bc_aggregation_sensitivity_percell.csv"), row.names=FALSE)
cat(sprintf("[saved] %s\n", file.path(PLOT_DIR, "decomp_bc_aggregation_sensitivity_percell.csv")))

p_focal <- plot_stacked_bar_focal_no_resp(
  blks_all = blks_all,
  block_levels = c("Cell type","Spatial cell state","Spillover","Residuals"),
  title = "BC (Janesick Xenium, V_spatial)")
ggsave(file.path(PLOT_DIR, "decomp_canonical_focal_bc_percell.pdf"),
        p_focal, width = 14, height = 6.5, device = cairo_pdf)
cat(sprintf("[saved] %s\n",
            file.path(PLOT_DIR, "decomp_canonical_focal_bc_percell.pdf")))

## --- Pairwise V_spatial heatmap via Pratt ---
pratt <- pace_pair_variance_pratt(
  mv, cond_prefix = NULL, cohort_label = "BC", block_label = "Spatial")
focal_spatial <- blks_all %>%
  filter(block == "Spatial cell state") %>%
  select(focal, focal_S = pct_total) %>%
  mutate(focal_S = focal_S * 100)

pair_pp <- pratt$pair_long %>%
  left_join(focal_spatial, by = "focal") %>%
  mutate(pair_pct = within_focal_share_pct * focal_S / 100,
         pair_pct = pmax(pair_pct, 0))

make_pair_tile <- function(d, ttl) {
  d %>%
    mutate(focal     = str_replace_all(focal, "[\\._]", " "),
           neighbour = str_replace_all(neighbour, "[\\._]", " "),
           pair_pct  = ifelse(focal == neighbour, NA_real_, pair_pct)) %>%
    ggplot(aes(x = focal, y = neighbour, fill = pair_pct)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", na.value = "white",
                          labels = function(x) sprintf("%.2f%%", x)) +
    labs(x = "Focal", y = "Neighbour", fill = "Variance\nshare",
          title = ttl) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    coord_equal()
}
p_pair <- make_pair_tile(pair_pp,
                          "BC: Spatial Cell State Contribution (V_spatial)")

## --- Top 5 genes per top 2 pairs (composite as canonical pairwise figure) ---
## Top-2 pairs MUST be ranked by pair_pct (the heatmap metric) so the left
## (heatmap) and right (per-pair top genes) panels of the composite figure are
## ordered consistently. Ranking by sum(top-N MCSD) gave a different order
## because MCSD weights effect size + spec + abundance, not statistical
## variance contribution.
N <- 5
pair_lookup <- pair_pp %>%
  filter(focal != neighbour) %>%
  mutate(pair = paste(focal, neighbour, sep = "_"),
         pair_pretty = sprintf("%s -> %s", focal, neighbour)) %>%
  select(pair, pair_pretty, pair_pct)
pair_order <- pair_lookup %>% arrange(desc(pair_pct)) %>%
  slice_head(n = 2) %>% pull(pair_pretty)

rows <- list()
for (p in names(mv$mcsd_canonical)) {
  s <- mv$mcsd_canonical[[p]]$scores
  if (is.null(s) || !nrow(s)) next
  pl <- pair_lookup$pair_pretty[match(p, pair_lookup$pair)]
  if (is.na(pl) || !(pl %in% pair_order)) next   # keep only the top-2 pairs
  rows[[p]] <- head(s, N) %>%
    mutate(pair = p, pair_pretty = pl,
            total_MCSD = sum(MCSD, na.rm = TRUE))
}
gene_df <- bind_rows(rows)
gene_df <- gene_df %>%
  filter(pair_pretty %in% pair_order) %>%
  mutate(pair_pretty = factor(pair_pretty, levels = pair_order),
         direction = ifelse(b_clean > 0, "Tum up", "Tum down"),
         stars = case_when(lfsr < 1e-10 ~ "***", lfsr < 1e-5 ~ "**",
                            lfsr < 0.05 ~ "*", TRUE ~ ""))
gene_df$label_y <- factor(
  paste(gene_df$pair_pretty,
        sprintf("%s  (%.0f%%)", gene_df$gene,
                100 * gene_df$MCSD / gene_df$total_MCSD), sep = "::"),
  levels = unlist(lapply(levels(gene_df$pair_pretty), function(pl) {
    sub <- gene_df %>% filter(pair_pretty == pl) %>% arrange(MCSD)
    paste(sub$pair_pretty,
          sprintf("%s  (%.0f%%)", sub$gene, 100 * sub$MCSD / sub$total_MCSD),
          sep = "::")
  })))

p_forest <- ggplot(gene_df,
                    aes(x = b_clean, y = label_y, color = direction, fill = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = b_clean - 1.96 * sd_shrunk,
                     xmax = b_clean + 1.96 * sd_shrunk),
                 height = 0.25, linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.3) +
  geom_text(aes(label = stars), color = "grey15", hjust = -0.5, size = 3.4,
            fontface = "bold") +
  scale_color_manual(values = c("Tum up" = "#B2182B", "Tum down" = "#2166AC"),
                       name = NULL) +
  scale_fill_manual(values = c("Tum up" = "#B2182B", "Tum down" = "#2166AC"),
                      name = NULL) +
  scale_y_discrete(labels = function(x) sub("^.*::", "", x)) +
  facet_wrap(~ pair_pretty, ncol = 1, scales = "free_y") +
  labs(title = "Top 5 genes per top pair",
       x = "PACE shrunken slope", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        strip.text = element_text(face = "bold", size = 11, hjust = 0),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.y = element_text(size = 9.5),
        legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

fig_pair <- p_pair + p_forest + plot_layout(widths = c(1, 1.1))
ggsave(file.path(PLOT_DIR, "decomp_canonical_pairwise_bc_percell.pdf"),
        fig_pair, width = 15, height = 9, device = cairo_pdf)
cat(sprintf("[saved] %s\n",
            file.path(PLOT_DIR, "decomp_canonical_pairwise_bc_percell.pdf")))

cat("\n=== Focal 4-block decomposition table (BC) ===\n")
print(blks_all %>% pivot_wider(names_from = block, values_from = pct_total) %>%
       mutate_if(is.numeric, ~round(. * 100, 2)))

cat("\n=== Top 10 pair_pct (V_spatial) ===\n")
print(pair_pp %>% arrange(desc(pair_pct)) %>%
       mutate(pair_pct = round(pair_pct, 3)) %>%
       select(focal, neighbour, V_pair_pratt, within_focal_share_pct, pair_pct) %>%
       head(10))

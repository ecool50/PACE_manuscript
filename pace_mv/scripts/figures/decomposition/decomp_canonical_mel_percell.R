## decomp_canonical_mel.R - Mel canonical decomposition figures using the
## locked manuscript functions:
##   - plot_stacked_bar_focal()     from 06a-plots-decomp.R
##   - make_pair_tile()              from figures/mvpql_pair_heatmap.R
##   - compute_pair_pct_long()       from figures/_compute_pair_pct.R
##
## Operates on Mel E^tech canonical fit (mvpql_percell_hc.rds).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(stringr); library(forcats)
  library(ggplot2); library(patchwork); library(scales); library(Matrix)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")

source("scripts/core/06a-plots-decomp.R")
source("scripts/helpers/pace_pair_variance_pratt.R")    ## LOCKED Pratt 1987 pair-variance (2026-05-28)

PLOT_DIR <- "plots/decomposition"
mv <- readRDS("data/simvi_melanoma/sweeps/mvpql_mel_streaming_nb1.rds")  # NB1 canonical (re-pointed 2026-06-24)
if (is.null(colnames(mv$fit$mu))) colnames(mv$fit$mu) <- mv$gene_set

## ============================================================================
## Build blks_all (focal × block × pct_total) for plot_stacked_bar_focal()
## ============================================================================
gf <- mv$decomposition$gene_focal_5block

## gene_focal_5block already stores per-gene pct_total per block (as %, not
## proportion). Aggregate to focal-mean and convert to proportion to match
## what plot_stacked_bar_focal expects (pct_total in 0-1).
## --- CANONICAL = unweighted POOLED aggregate (Goldstein-Browne-Rasbash 2002):
## pct_m = sum_j V_{m,j} / sum_j Total_j. The mean-of-ratios and inverse-
## residual-weighted variants go to the sensitivity CSV only, NOT the figure. ---
EPS <- 1e-9
.agg_mel <- function(d, scheme) {
  if (scheme == "mean")
    return(tibble::tibble(
      `Cell type`               = mean(d[["Cell type %"]],               na.rm=TRUE)/100,
      `Spatial cell state`      = mean(d[["Spatial state %"]],           na.rm=TRUE)/100,
      `Responder spatial state` = mean(d[["Responder spatial state %"]], na.rm=TRUE)/100,
      `Spillover`               = mean(d[["Spillover %"]],               na.rm=TRUE)/100,
      `Residuals`               = mean(d[["Residual %"]],                na.rm=TRUE)/100))
  w   <- if (scheme == "weighted") 1/pmax(d$V_disp, EPS) else rep(1, nrow(d))
  tot <- sum(w * d$Total)
  tibble::tibble(
    `Cell type`               = sum(w*d$celltype_offset_sq)/tot,
    `Spatial cell state`      = sum(w*d$V_state_baseline)/tot,
    `Responder spatial state` = sum(w*d$V_state_responder)/tot,
    `Spillover`               = sum(w*d$V_spill)/tot,
    `Residuals`               = sum(w*d$V_disp)/tot)
}
## CANONICAL block bar = SINGLE-FRAME (observed scale), pooled by denom -> BC-consistent.
## (The link-scale .agg_mel above is retained only for the sensitivity CSV below.)
sf <- mv$decomposition$gene_focal_single_frame
.agg_mel_sf <- function(d) {
  tot <- sum(d$denom)
  tibble::tibble(
    `Cell type`               = sum(d$`Cell type %`/100         * d$denom)/tot,
    `Spatial cell state`      = sum(d$`Spatial %`/100           * d$denom)/tot,
    `Responder spatial state` = sum(d$`Responder spatial %`/100 * d$denom)/tot,
    `Spillover`               = sum(d$`Spillover %`/100         * d$denom)/tot,
    `Residuals`               = sum(d$`Residual %`/100          * d$denom)/tot)
}
blks_all <- sf %>% group_by(focal) %>% group_modify(~ .agg_mel_sf(.x)) %>%
  ungroup() %>% pivot_longer(-focal, names_to = "block", values_to = "pct_total")
.sens <- bind_rows(
  gf %>% group_by(focal) %>% group_modify(~ .agg_mel(.x,"weighted")) %>% mutate(scheme="weighted_pooled"),
  gf %>% group_by(focal) %>% group_modify(~ .agg_mel(.x,"pooled"))   %>% mutate(scheme="unweighted_pooled"),
  gf %>% group_by(focal) %>% group_modify(~ .agg_mel(.x,"mean"))     %>% mutate(scheme="mean_of_ratios")) %>% ungroup()
dir.create(PLOT_DIR, recursive=TRUE, showWarnings=FALSE)
utils::write.csv(.sens, file.path(PLOT_DIR, "decomp_mel_aggregation_sensitivity_percell.csv"), row.names=FALSE)
cat(sprintf("[saved] %s\n", file.path(PLOT_DIR, "decomp_mel_aggregation_sensitivity_percell.csv")))

p_focal <- plot_stacked_bar_focal(
  blks_all = blks_all,
  block_levels = c("Cell type","Spatial cell state","Responder spatial state",
                    "Spillover","Residuals"),
  title = "Mel (CosMx, Responder = PD vs RESP)",
  zoom_blocks = c("Spatial cell state","Responder spatial state","Spillover")
)
ggsave(file.path(PLOT_DIR, "decomp_canonical_focal_mel_percell.pdf"),
        p_focal, width = 14, height = 6.5, device = cairo_pdf)
cat(sprintf("[saved] %s\n",
            file.path(PLOT_DIR, "decomp_canonical_focal_mel_percell.pdf")))

## ============================================================================
## Pairwise V_RxS heatmap via compute_pair_pct_long + make_pair_tile
## (from mvpql_pair_heatmap.R)
## ============================================================================
make_pair_tile <- function(d, ttl,
                            fill_lbl = "Variance share",
                            subtitle = NULL) {
  d |>
    dplyr::mutate(focal     = stringr::str_replace_all(focal, "[\\._]", " "),
                  neighbour = stringr::str_replace_all(neighbour, "[\\._]", " "),
                  pair_pct  = ifelse(focal == neighbour, NA_real_, pair_pct)) |>
    ggplot(aes(x = focal, y = neighbour, fill = pair_pct)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", na.value = "white",
                          labels = function(x) sprintf("%.2f%%", x)) +
    labs(x = "Focal", y = "Neighbour", fill = fill_lbl,
          title = ttl, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          plot.title       = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle    = element_text(hjust = 0.5, colour = "grey30",
                                            size = 9, margin = margin(b = 6))) +
    coord_equal()
}

## --- Pratt 1987 V_pair (locked canonical pair-level decomp) ---
## V_pair(c,t,g) = u_{c,t,g} · (Σ_K(c) u_{c,g})_t -- sums exactly to V_spatial(c,g)
## Express as % of pooled within-celltype link-scale variance via the focal
## V_RxS share (from gene_focal_5block) to match the manuscript metric.
pratt_res <- pace_pair_variance_pratt(
  mv, cond_prefix = "ResponderPD", cohort_label = "Mel", block_label = "RxS")
focal_rxs_pct <- blks_all |>
  dplyr::filter(block == "Responder spatial state") |>
  dplyr::select(focal, focal_rxs = pct_total) |>
  dplyr::mutate(focal_rxs = focal_rxs * 100)  ## back to percent

mel_pp <- pratt_res$pair_long |>
  dplyr::left_join(focal_rxs_pct, by = "focal") |>
  dplyr::mutate(
    block = "responder_spatial",
    ## within_focal_share_pct is signed % within focal V_RxS;
    ## multiply by focal V_RxS as % of total variance to get % of pooled total.
    pair_pct = within_focal_share_pct * focal_rxs / 100
  ) |>
  dplyr::select(focal, neighbour, block, pair_pct)

cat("\n[Mel] Top 5 (focal -> neighbour) Resp × Spatial pairs (PRATT canonical):\n")
print(mel_pp |> dplyr::filter(focal != neighbour) |>
        dplyr::arrange(dplyr::desc(pair_pct)) |> head(5))

p_mel_resp <- make_pair_tile(
  mel_pp,
  "Melanoma: Responder × Spatial Contribution (PD vs RESP)",
  subtitle = "% of pooled within-celltype link-scale variance attributable to (PD × neighbour) interaction")

## --- Top 5 genes per top 2 pairs (composite as canonical pairwise figure) ---
## Top-2 pairs MUST be ranked by pair_pct (the heatmap metric) so the left
## (heatmap) and right (per-pair top genes) panels of the composite figure are
## ordered consistently. Ranking by sum(top-N MCSD) gave a different order
## because MCSD weights effect size + spec + abundance, not statistical
## variance contribution.
N <- 5
pair_lookup <- mel_pp %>%
  dplyr::filter(focal != neighbour) %>%
  dplyr::mutate(pair = paste(focal, neighbour, sep = "_"),
                 pair_pretty = sprintf("%s -> %s", focal, neighbour)) %>%
  dplyr::select(pair, pair_pretty, pair_pct)
pair_order <- pair_lookup %>% dplyr::arrange(dplyr::desc(pair_pct)) %>%
  dplyr::slice_head(n = 2) %>% dplyr::pull(pair_pretty)

## FIX (2026-06-11): source top genes for each top RxS pair from shrunken_long (the
## ResponderPD:neighbour slope), NOT mcsd_canonical (which stores X->Tumour pairs and
## doesn't cover the Tumour->X / Fibroblast->Endothelial pairs in pair_order). Rank by
## confect = max(0, |slope| - 1.645*sd), lfsr<0.05.
Zc <- 1.645
rows_g <- list()
for (pl in pair_order) {
  parts <- strsplit(pl, " -> ", fixed = TRUE)[[1]]; f <- parts[1]; h <- parts[2]
  s <- mv$shrunken_long %>%
    dplyr::filter(focal == f, neighbour == h, term == paste0("ResponderPD:", h)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(b_clean = estimate_shrunk,
                  MCSD    = pmax(0, abs(estimate_shrunk) - Zc * sd_shrunk)) %>%
    dplyr::filter(lfsr < 0.05, MCSD > 0) %>%
    dplyr::arrange(dplyr::desc(MCSD)) %>% head(N)
  if (!nrow(s)) next
  rows_g[[pl]] <- s %>% dplyr::mutate(pair = pl, pair_pretty = pl, total_MCSD = sum(MCSD, na.rm = TRUE))
}
gene_df <- dplyr::bind_rows(rows_g)
if (nrow(gene_df) == 0) {
  ggsave(file.path(PLOT_DIR, "decomp_canonical_pairwise_mel_percell.pdf"), p_mel_resp,
         width = 9, height = 7, device = cairo_pdf)
  cat("[saved heatmap-only] no significant RxS genes for the top pairs\n"); quit(save = "no")
}
gene_df <- gene_df %>%
  dplyr::filter(pair_pretty %in% pair_order) %>%
  dplyr::mutate(pair_pretty = factor(pair_pretty, levels = pair_order),
                 direction = ifelse(b_clean > 0, "PD up", "PD down"),
                 stars = dplyr::case_when(lfsr < 1e-10 ~ "***",
                                            lfsr < 1e-5 ~ "**",
                                            lfsr < 0.05 ~ "*", TRUE ~ ""))
gene_df$label_y <- factor(
  paste(gene_df$pair_pretty,
        sprintf("%s  (%.0f%%)", gene_df$gene,
                100 * gene_df$MCSD / gene_df$total_MCSD), sep = "::"),
  levels = unlist(lapply(levels(gene_df$pair_pretty), function(pl) {
    sub <- gene_df %>% dplyr::filter(pair_pretty == pl) %>% dplyr::arrange(MCSD)
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
  scale_color_manual(values = c("PD up" = "#B2182B", "PD down" = "#2166AC"),
                       name = NULL) +
  scale_fill_manual(values = c("PD up" = "#B2182B", "PD down" = "#2166AC"),
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

fig_pair <- p_mel_resp + p_forest + plot_layout(widths = c(1, 1.1))
ggsave(file.path(PLOT_DIR, "decomp_canonical_pairwise_mel_percell.pdf"),
        fig_pair, width = 15, height = 9, device = cairo_pdf)
cat(sprintf("[saved] %s\n",
            file.path(PLOT_DIR, "decomp_canonical_pairwise_mel_percell.pdf")))

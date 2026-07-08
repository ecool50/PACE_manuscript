## decomp_eta2_total_bc.R — full ANOVA eta^2 partition, normalized by TOTAL variance.
## Each focal column = [neighbour spatial % of total] + [Non-spatial residual % of total],
## summing to 1. Magnitude is INTRINSIC: a focal with no real spatial signal has ~0 neighbour
## mass (phantom impossible). Spatial part = cleaned (spec>=0.5) RWA shares x cleaned spatial%.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
mv <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
gf <- mv$decomposition$gene_focal_4block
share <- readRDS("streaming/neighbour_share_rwa_cleaned_bc.rds")   # within-spatial shares (sum 1)
TYPES <- rownames(share)
sclean <- gf |> group_by(focal) |> summarise(s = 100*sum(V_state[spec>=0.5])/sum(Total), .groups="drop")
ford <- sclean |> arrange(desc(s)) |> pull(focal)

## p_k = neighbour's % of TOTAL variance = within-spatial share x cleaned spatial %
long <- as.data.frame(as.table(share)) |> rename(w=Freq) |>
  left_join(sclean, by="focal") |>
  mutate(p = w * s,                                   # % of total
         p = ifelse(as.character(focal)==as.character(neighbour), NA, p))
nonsp <- sclean |> transmute(focal, neighbour="Non-spatial", p = 100 - s)  # residual+celltype+spillover
d <- bind_rows(long |> select(focal,neighbour,p), nonsp) |>
  mutate(focal=factor(focal, levels=ford),
         neighbour=factor(neighbour, levels=c(rev(TYPES),"Non-spatial")),
         kind=ifelse(neighbour=="Non-spatial","nonsp","spatial"))

cap <- quantile(d$p[d$kind=="spatial"], 0.99, na.rm=TRUE)
sp <- d |> filter(kind=="spatial"); ns <- d |> filter(kind=="nonsp")
p <- ggplot() +
  geom_tile(data=sp, aes(focal, neighbour, fill=pmin(p,cap)), colour="grey92", linewidth=0.3) +
  geom_text(data=sp, aes(focal, neighbour, label=ifelse(is.na(p)|p<0.02,"",sprintf("%.2f",p))), size=2.7) +
  geom_tile(data=ns, aes(focal, neighbour), fill="grey85", colour="grey70", linewidth=0.3) +
  geom_text(data=ns, aes(focal, neighbour, label=sprintf("%.1f",p)), size=2.6, colour="grey30") +
  scale_fill_viridis_c(option="magma", na.value="white", name="% of total\nvariance", limits=c(0,cap)) +
  labs(x="Focal (column sums to 100%)", y="Neighbour",
       title="BC 313: ANOVA eta^2 spatial partition (% of TOTAL variance, columns sum to 1)",
       subtitle="neighbour cells = spatial variance that neighbour explains; grey = non-spatial residual. Magnitude is intrinsic: no-signal focals are dark -> no phantom.") +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.subtitle=element_text(size=8))
ggsave("streaming/figures/decomp_eta2_total_bc.pdf", p, width=9, height=7.5, device=cairo_pdf)
cat("wrote streaming/figures/decomp_eta2_total_bc.pdf; spatial color cap =", round(cap,3),"% of total\n")

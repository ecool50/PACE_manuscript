## decomp_nested_cleaned_bc.R — nested cleaned spatial decomposition figure.
## TOP bar  = Level-1: each focal's cleaned spatial-block % (magnitude; spec>=0.5 genes).
## HEATMAP  = Level-2: RWA neighbour shares WITHIN the spatial block (columns sum to 1).
## Read together: a column sums to 1, but the bar above tells you how big the block it divides is.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr)})
mv <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
gf <- mv$decomposition$gene_focal_4block
share <- readRDS("streaming/neighbour_share_rwa_cleaned_bc.rds")

## Level-1 cleaned spatial-block % per focal; order focals by it (magnitude descending)
fs <- gf |> group_by(focal) |>
  summarise(spatial = 100*sum(V_state[spec>=0.5])/sum(Total), .groups="drop") |>
  arrange(desc(spatial))
ford <- fs$focal

## ---- TOP: Level-1 magnitude bar ----
fs <- fs |> mutate(focal=factor(focal, levels=ford),
                   low = spatial < 0.1)
pTop <- ggplot(fs, aes(focal, spatial, fill=low)) +
  geom_col(width=0.8) +
  geom_text(aes(label=sprintf("%.2f", spatial)), vjust=-0.3, size=2.8) +
  scale_fill_manual(values=c("FALSE"="#1F3B57","TRUE"="grey70"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.25))) +
  labs(y="Spatial\nblock %", x=NULL,
       title="BC 313: nested contamination-cleaned spatial decomposition",
       subtitle="TOP = Level-1 spatial-block magnitude per focal (grey = <0.1%, interpret column with caution)  |  HEATMAP = Level-2 neighbour shares within the block (columns sum to 1)") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.subtitle=element_text(size=8))

## ---- BOTTOM: Level-2 column-sum-to-1 neighbour shares ----
d <- as.data.frame(as.table(share)) |> rename(s=Freq) |>
  mutate(focal=factor(focal, levels=ford), neighbour=factor(neighbour, levels=rev(ford)),
         s=ifelse(as.character(focal)==as.character(neighbour), NA, s))
pHM <- ggplot(d, aes(focal, neighbour, fill=s)) +
  geom_tile(colour="grey92", linewidth=0.3) +
  geom_text(aes(label=ifelse(is.na(s)|s<0.02,"",sprintf("%.0f",100*s))), size=3) +
  scale_fill_viridis_c(option="magma", na.value="white", limits=c(0,1),
                       labels=scales::percent, name="Neighbour\nshare") +
  labs(x="Focal (column sums to 1)", y="Neighbour") +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1))

comp <- pTop / pHM + plot_layout(heights=c(1,4))
ggsave("streaming/figures/decomp_nested_cleaned_bc.pdf", comp, width=9, height=8, device=cairo_pdf)
cat("wrote streaming/figures/decomp_nested_cleaned_bc.pdf\n")

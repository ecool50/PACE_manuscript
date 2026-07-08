## decomp_shrunk_shares_bc.R — column-sum-to-1 neighbour shares, shrunk toward uniform by
## spatial-block magnitude. Weak focals -> ~uniform (no false dominant neighbour); strong
## focals -> data. Fix is IN the proportions; columns sum to 1 over neighbours (no extra row).
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
mv <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
gf <- mv$decomposition$gene_focal_4block
share <- readRDS("streaming/neighbour_share_rwa_cleaned_bc.rds")   # RWA cleaned within-spatial shares
TYPES <- rownames(share)
sF <- gf |> group_by(focal) |> summarise(s = 100*sum(V_state[spec>=0.5])/sum(Total), .groups="drop")
alpha <- median(sF$s)                                              # data-derived shrinkage scale
cat(sprintf("shrinkage scale alpha = median focal spatial%% = %.3f\n", alpha))

ford <- sF |> arrange(desc(s)) |> pull(focal)
shr <- matrix(0, length(TYPES), length(TYPES), dimnames=list(neighbour=TYPES, focal=TYPES))
for (F in TYPES) {
  w <- share[,F]; kept <- which(w > 0); K <- length(kept); if(K<1) next
  s <- sF$s[sF$focal==F]; lam <- s/(s+alpha)
  shr[kept,F] <- lam*w[kept] + (1-lam)/K                          # shrink toward uniform-over-kept
}
saveRDS(shr, "streaming/neighbour_share_shrunk_bc.rds")
cat("\nlambda (data weight) per focal:\n")
for (F in ford) { s<-sF$s[sF$focal==F]; cat(sprintf("  %-15s s=%.2f%%  lambda=%.2f\n", F, s, s/(s+alpha))) }

d <- as.data.frame(as.table(shr)) |> rename(p=Freq) |>
  mutate(focal=factor(focal,levels=ford), neighbour=factor(neighbour,levels=rev(ford)),
         p=ifelse(as.character(focal)==as.character(neighbour),NA,p))
pg <- ggplot(d, aes(focal, neighbour, fill=p)) +
  geom_tile(colour="grey92", linewidth=0.3) +
  geom_text(aes(label=ifelse(is.na(p)|p<0.02,"",sprintf("%.0f",100*p))), size=3) +
  scale_fill_viridis_c(option="magma", na.value="white", limits=c(0,1), labels=scales::percent, name="Share") +
  coord_equal() +
  labs(x="Focal (columns sum to 1)", y="Neighbour",
       title="BC 313: neighbour shares shrunk toward uniform by spatial magnitude",
       subtitle=sprintf("weak-signal focals -> ~uniform (no false dominant neighbour); strong focals -> data. alpha=median spatial%%=%.2f", alpha)) +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.subtitle=element_text(size=8))
ggsave("streaming/figures/decomp_shrunk_shares_bc.pdf", pg, width=8.2, height=7, device=cairo_pdf)
cat("\nwrote streaming/figures/decomp_shrunk_shares_bc.pdf\n")

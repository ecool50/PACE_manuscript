## panelB_spp1_superplot.R — GATE robustness check for SPP1 Macrophage->Fibroblast (RxS).
## Per-patient slope SuperPlot (Lord 2020) + leave-one-patient-out contrast stability.
## Viz aid: per-patient raw slopes (external regression OK as viz aid; PACE slope is the headline).
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
  library(SingleCellExperiment); library(SummarizedExperiment); library(zellkonverter)
})
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
mv <- readRDS("data/simvi_melanoma/sweeps/mvpql_percell_hc.rds")
cm <- mv$cell_meta
fibprox <- mv$fit$re_meta$Z[, "Macrophage::Fibroblast"]      # within-image-centered fibroblast kernel

## raw SPP1 from h5ad
sce <- zellkonverter::readH5AD("data/simvi_melanoma/Melanoma_5612.h5ad", reader="R", verbose=FALSE)
Yall <- as.matrix(assay(sce,"X")); rm(sce)
cid <- paste0("Cell", cm$cell_id); mi <- match(cid, colnames(Yall))
spp1 <- rep(NA_real_, nrow(cm)); ok <- !is.na(mi)
spp1[ok] <- Yall["SPP1", mi[ok]]
norm <- log1p(spp1 / pmax(cm$nCount,1) * 1e3)               # log normalized SPP1

d <- tibble(img=as.character(cm$imageID), cond=as.character(cm$Condition),
            ct=as.character(cm$celltype), fib=as.numeric(fibprox), spp1=norm) |>
     filter(ct=="Macrophage", is.finite(spp1), is.finite(fib))

MINC <- 20
per_pat <- d |> group_by(img,cond) |>
  filter(n() >= MINC, sd(fib)>0) |>
  summarise(slope = coef(lm(spp1 ~ fib))[2], n=n(), .groups="drop")
cat(sprintf("per-patient slopes: %d patients (PD=%d, RESP=%d); min cells=%d\n",
            nrow(per_pat), sum(per_pat$cond=="PD"), sum(per_pat$cond=="RESP"), MINC))
print(per_pat |> arrange(cond, slope), n=30)
wt <- suppressWarnings(wilcox.test(slope ~ cond, per_pat))
cat(sprintf("\nWilcoxon PD vs RESP per-patient slopes: p=%.3f ; median PD=%.3f RESP=%.3f\n",
            wt$p.value, median(per_pat$slope[per_pat$cond=="PD"]), median(per_pat$slope[per_pat$cond=="RESP"])))

## LOPO: drop each patient, recompute median(RESP)-median(PD) contrast
full_c <- median(per_pat$slope[per_pat$cond=="RESP"]) - median(per_pat$slope[per_pat$cond=="PD"])
lopo <- sapply(seq_len(nrow(per_pat)), function(i){ s<-per_pat[-i,]
  median(s$slope[s$cond=="RESP"]) - median(s$slope[s$cond=="PD"]) })
cat(sprintf("LOPO contrast (RESP-PD): full=%.3f ; range over drops=[%.3f, %.3f] ; sign flips=%d\n",
            full_c, min(lopo), max(lopo), sum(sign(lopo)!=sign(full_c))))

## ---- plots ----
TB<-"#4E79A7"; RD<-"#E15759"; pal<-c(PD=RD, RESP=TB)
pA <- ggplot(per_pat, aes(cond, slope, color=cond)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") +
  geom_jitter(aes(size=n), width=0.12, alpha=0.85) +
  stat_summary(fun=median, geom="crossbar", width=0.4, fatten=2, color="black") +
  scale_color_manual(values=pal, guide="none") + scale_size_area(max_size=7, name="n macroph.") +
  labs(title=sprintf("A  Per-patient SPP1~Fib-proximity slope (Wilcoxon p=%.3f)",wt$p.value),
       x=NULL, y="per-patient slope") + theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())
pB <- ggplot(d, aes(fib, spp1, color=cond)) +
  geom_point(alpha=0.12, size=0.5) + geom_smooth(method="lm", se=TRUE, linewidth=0.9) +
  scale_color_manual(values=pal, name=NULL) +
  labs(title="B  Cell-level: SPP1 vs Fibroblast proximity", x="Fibroblast proximity (within-image)", y="log-norm SPP1") +
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(), legend.position=c(.98,.98), legend.justification=c(1,1))
pC <- ggplot(tibble(x=seq_along(lopo), c=lopo), aes(x,c)) +
  geom_hline(yintercept=full_c, linetype="dashed", color=TB) + geom_hline(yintercept=0, color="grey60") +
  geom_point(size=2, color="black") +
  labs(title="C  Leave-one-patient-out contrast (RESP - PD)", x="patient dropped", y="contrast") +
  theme_bw(base_size=11)+theme(panel.grid.minor=element_blank())
ggsave("plots/diagnostics/panelB_spp1_macfib_gate.pdf", (pA|pB)/pC, width=11, height=8, device=cairo_pdf)
cat("\n[saved] plots/diagnostics/panelB_spp1_macfib_gate.pdf\n")

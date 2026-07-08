## decomp_anova_ss_bc.R — ANOVA-style SS-ratio neighbour decomposition (Ellis spec).
## Columns sum to 1. NO resampling/Shapley/RWA — just variance terms.
## V_pair_k = [sum_g (b_kg^2 + se_kg^2)] * w_k,  w_k = proximity-variance weight.
##   Lever 1 (denominator = spatial block): share_k = V_pair_k / sum_j V_pair_j.
##   Lever 2 (de-emphasize abundance): w_k = Var(N_k | N_k>0)  (exposed cells only)
##                                     vs realized w_k = Var(N_k) (all cells).
suppressPackageStartupMessages({library(Matrix); library(dbscan); library(dplyr); library(ggplot2); library(tidyr)})
say <- function(...) cat(sprintf(...),"\n")
mv <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
fit <- mv$fit; U <- fit$U; seU <- fit$se_U
TYPES <- sort(unique(as.character(mv$cell_meta$celltype)))
ct <- as.character(mv$cell_meta$celltype)
ymv <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds"); df <- ymv$df

## raw (uncentered) K_bio kernel per cell x type (h_bio=30, RAD=90; BC = one section -> global)
xy <- as.matrix(df[,c("x","y")]); n <- nrow(df); H <- 30; RAD <- 90
fr <- dbscan::frNN(xy, eps=RAD)
iv <- rep.int(seq_len(n), lengths(fr$id)); jv <- unlist(fr$id, use.names=FALSE); dv <- unlist(fr$dist, use.names=FALSE)
jc <- match(ct[jv], TYPES); ok <- !is.na(jc)
A <- as.matrix(sparseMatrix(i=iv[ok], j=jc[ok], x=exp(-(dv[ok])^2/H^2), dims=c(n,length(TYPES)))); colnames(A) <- TYPES
rm(fr,iv,jv,dv,jc,ok); invisible(gc())

realized <- exposed <- matrix(0,length(TYPES),length(TYPES),dimnames=list(neighbour=TYPES,focal=TYPES))
for (F in TYPES) {
  cF <- which(ct==F)
  for (k in setdiff(TYPES,F)) {
    row <- paste0(F,"::",k); if(!(row %in% rownames(U))) next
    Bk <- sum(U[row,]^2 + seU[row,]^2)               # gene-summed effect energy
    a <- A[cF, k]; ex <- a > 1e-8
    realized[k,F] <- Bk * var(a)
    exposed[k,F]  <- Bk * (if(sum(ex)>2) var(a[ex]) else 0)
  }
}
norm <- function(M){ cs<-colSums(M); sweep(M,2,ifelse(cs>0,cs,1),"/") }
Rn <- norm(realized); En <- norm(exposed)
saveRDS(list(realized=Rn, exposed=En), "streaming/neighbour_share_anova_bc.rds")

## ANOVA table (exposure-conditional = both levers), columns sum to 1
say("=== ANOVA SS-ratio table: share of focal spatial variance by neighbour (Lever 1+2, exposure-conditional) ===")
tab <- as.data.frame(round(100*En,1)); say("(columns sum to 100%%)")
print(tab)
say("\ncolumn sums (check): %s", paste(sprintf("%.0f",colSums(100*En)),collapse=" "))

## heatmap (exposure-conditional)
d <- as.data.frame(as.table(En)) |> rename(s=Freq) |>
  mutate(focal=factor(focal,levels=TYPES), neighbour=factor(neighbour,levels=rev(TYPES)),
         s=ifelse(as.character(focal)==as.character(neighbour),NA,s))
p <- ggplot(d, aes(focal, neighbour, fill=s)) + geom_tile(colour="grey92", linewidth=0.3) +
  geom_text(aes(label=ifelse(is.na(s)|s<0.02,"",sprintf("%.0f",100*s))), size=3) +
  scale_fill_viridis_c(option="inferno", na.value="white", limits=c(0,1), labels=scales::percent, name="Spatial\nSS share") +
  coord_equal() +
  labs(x="Focal (column sums to 1)", y="Neighbour",
       title="BC 313: ANOVA SS-ratio spatial decomposition (Ellis spec)",
       subtitle="share = b^2 * Var(N|exposed) per neighbour / spatial block; columns sum to 1; exposure-conditional de-emphasizes abundance") +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.subtitle=element_text(size=8))
ggsave("streaming/figures/decomp_anova_ss_bc.pdf", p, width=8.2, height=7, device=cairo_pdf)
say("\nwrote streaming/figures/decomp_anova_ss_bc.pdf")

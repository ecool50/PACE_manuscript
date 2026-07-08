## decomp_anova_ss_weighted_bc.R — Ellis ANOVA SS-ratio + spec^2 reliability weight.
## V_pair_k = sum_g spec_g^2 * (b_kg^2 + se_kg^2) * Var(N_k | exposed); share_k col-normalized.
## spec^2 = the discriminating statistic that un-drowns precise contamination. Shows the
## per-focal gene ranking BEFORE (raw V_state) vs AFTER (spec^2-weighted) to prove signal recovery.
suppressPackageStartupMessages({library(Matrix); library(dbscan); library(dplyr); library(ggplot2)})
say <- function(...) cat(sprintf(...),"\n")
mv <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
fit <- mv$fit; U <- fit$U; seU <- fit$se_U
TYPES <- sort(unique(as.character(mv$cell_meta$celltype))); ct <- as.character(mv$cell_meta$celltype)
gf <- mv$decomposition$gene_focal_4block; gs <- mv$gene_set
ymv <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds"); df <- ymv$df

## raw K_bio kernel per cell x type (h_bio=30, RAD=90, global)
xy <- as.matrix(df[,c("x","y")]); n <- nrow(df); H <- 30
fr <- dbscan::frNN(xy, eps=90); iv <- rep.int(seq_len(n), lengths(fr$id)); jv <- unlist(fr$id,use.names=FALSE); dv <- unlist(fr$dist,use.names=FALSE)
jc <- match(ct[jv],TYPES); ok <- !is.na(jc)
A <- as.matrix(sparseMatrix(i=iv[ok],j=jc[ok],x=exp(-(dv[ok])^2/H^2),dims=c(n,length(TYPES)))); colnames(A)<-TYPES
rm(fr,iv,jv,dv,jc,ok); invisible(gc())

share <- matrix(0,length(TYPES),length(TYPES),dimnames=list(neighbour=TYPES,focal=TYPES))
say("=== per-focal gene ranking: BEFORE (raw V_state) vs AFTER (spec^2-weighted) ===")
for (F in TYPES) {
  cF <- which(ct==F); spec <- gf$spec[gf$focal==F][match(gs, gf$gene[gf$focal==F])]; spec[is.na(spec)] <- 0
  nbr <- setdiff(TYPES,F)
  vcond <- setNames(sapply(nbr, function(k){a<-A[cF,k];ex<-a>1e-8; if(sum(ex)>2) var(a[ex]) else 0}), nbr)
  # per gene: spatial V contribution summed over neighbours
  geneV <- numeric(length(gs))
  Vpair <- setNames(numeric(length(nbr)), nbr)
  for (k in nbr) { row<-paste0(F,"::",k); if(!(row%in%rownames(U))) next
    contr <- (U[row,]^2 + seU[row,]^2) * vcond[k]      # per gene, this neighbour
    geneV <- geneV + contr
    Vpair[k] <- sum(spec^2 * contr) }                   # spec^2-weighted pair variance
  geneV_wt <- spec^2 * geneV
  if (sum(Vpair)>0) share[nbr,F] <- Vpair/sum(Vpair)
  if (F %in% c("Stromal","Macrophage","Endothelial","Tumour")) {
    raw5 <- gs[order(geneV,decreasing=TRUE)][1:5]; wt5 <- gs[order(geneV_wt,decreasing=TRUE)][1:5]
    sp <- function(g) spec[match(g,gs)]
    say("  %-12s BEFORE: %s", F, paste(sprintf("%s[%.2f]",raw5,sp(raw5)),collapse=", "))
    say("  %-12s AFTER : %s", F, paste(sprintf("%s[%.2f]",wt5,sp(wt5)),collapse=", "))
  }
}
saveRDS(share, "streaming/neighbour_share_anova_weighted_bc.rds")

d <- as.data.frame(as.table(share)) |> rename(s=Freq) |>
  mutate(focal=factor(focal,levels=TYPES), neighbour=factor(neighbour,levels=rev(TYPES)),
         s=ifelse(as.character(focal)==as.character(neighbour),NA,s))
p <- ggplot(d, aes(focal, neighbour, fill=s)) + geom_tile(colour="grey92", linewidth=0.3) +
  geom_text(aes(label=ifelse(is.na(s)|s<0.02,"",sprintf("%.0f",100*s))), size=3) +
  scale_fill_viridis_c(option="inferno", na.value="white", limits=c(0,1), labels=scales::percent, name="Spatial\nSS share") +
  coord_equal() +
  labs(x="Focal (column sums to 1)", y="Neighbour",
       title="BC 313: spec^2-weighted ANOVA SS decomposition (signal un-drowned)",
       subtitle="V_pair = sum_g spec^2 * b^2 * Var(N|exposed); columns sum to 1; contamination genes down-weighted -> real programs surface") +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.subtitle=element_text(size=8))
ggsave("streaming/figures/decomp_anova_ss_weighted_bc.pdf", p, width=8.2, height=7, device=cairo_pdf)
say("\nwrote streaming/figures/decomp_anova_ss_weighted_bc.pdf")

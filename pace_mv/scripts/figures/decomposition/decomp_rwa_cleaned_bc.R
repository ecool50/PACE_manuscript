## decomp_rwa_cleaned_bc.R — contamination-cleaned decomposition.
## Keep only spec>=0.5 genes (focal is the gene's top expresser -> not neighbour bleed),
## recompute Spatial-state % and RWA neighbour shares on that set. Raw (all genes) = upper-bound.
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
say <- function(...) cat(sprintf(...),"\n")
mv  <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
fit <- mv$fit; U <- fit$U; seU <- fit$se_U; Z <- fit$re_meta$Z
TYPES <- sort(unique(as.character(mv$cell_meta$celltype)))
gf <- mv$decomposition$gene_focal_4block
sym_sqrt <- function(M){e<-eigen(M,symmetric=TRUE);e$vectors%*%diag(sqrt(pmax(e$values,0)),ncol(M))%*%t(e$vectors)}

rwa_share <- function(F, keep_genes) {
  colI <- paste0(F,"::(Intercept)"); cF <- which(Z[,colI]!=0)
  nb_cols <- paste0(F,"::",setdiff(TYPES,F)); nb_cols <- nb_cols[nb_cols %in% colnames(Z)]
  N <- as.matrix(Z[cF, nb_cols, drop=FALSE]); keepc <- which(apply(N,2,var)>1e-12)
  if(length(keepc)<2) return(NULL)
  N <- N[,keepc,drop=FALSE]; cols <- nb_cols[keepc]; nbr <- sub(paste0("^",F,"::"),"",cols)
  Sig <- stats::cov(N); D <- diag(Sig); Dh <- sqrt(D); R <- Sig/outer(Dh,Dh); Rh <- sym_sqrt(R)
  gi <- match(keep_genes, mv$gene_set); gi <- gi[!is.na(gi)]
  Uc <- U[cols, gi, drop=FALSE]; Sc <- seU[cols, gi, drop=FALSE]
  B  <- Uc %*% t(Uc) + diag(rowSums(Sc^2), nrow=length(cols))
  Gst <- (Dh*B)*rep(Dh, each=length(Dh))
  w <- as.numeric((Rh*Rh) %*% diag(Rh %*% Gst %*% Rh))
  setNames(pmax(w,0)/sum(pmax(w,0)), nbr)
}

share_clean <- matrix(0,length(TYPES),length(TYPES),dimnames=list(neighbour=TYPES,focal=TYPES))
say("%-15s %7s %10s %10s | %-26s %-26s", "focal","n_keep","Spatial%raw","Spatial%cln","top-2 RAW","top-2 CLEANED")
raw <- readRDS("streaming/neighbour_share_rwa_bc.rds")
for (F in TYPES) {
  d <- gf |> filter(focal==F)
  keep <- d$gene[d$spec>=0.5]
  sp_raw <- 100*sum(d$V_state)/sum(d$Total)
  sp_cln <- 100*sum(d$V_state[d$spec>=0.5])/sum(d$Total)
  sh <- rwa_share(F, keep)
  if(!is.null(sh)) share_clean[names(sh),F] <- sh
  t2 <- function(x){x<-x[x>0]; if(!length(x)) return("-"); v<-sort(x,decreasing=TRUE)[1:min(2,length(x))]; paste(sprintf("%s %.0f%%",names(v),100*v),collapse=", ")}
  say("%-15s %7d %10.2f %10.2f | %-26s %-26s", F, length(keep), sp_raw, sp_cln, t2(raw[,F]), t2(share_clean[,F]))
}
saveRDS(share_clean, "streaming/neighbour_share_rwa_cleaned_bc.rds")

## cleaned heatmap
d <- as.data.frame(as.table(share_clean)) |> rename(s=Freq) |>
  mutate(focal=factor(focal,levels=TYPES), neighbour=factor(neighbour,levels=rev(TYPES)),
         s=ifelse(as.character(focal)==as.character(neighbour),NA,s))
p <- ggplot(d, aes(focal, neighbour, fill=s)) + geom_tile(colour="grey92", linewidth=0.3) +
  geom_text(aes(label=ifelse(is.na(s)|s<0.02,"",sprintf("%.0f",100*s))), size=3) +
  scale_fill_viridis_c(option="magma", na.value="white", limits=c(0,1), labels=scales::percent, name="Share") +
  coord_equal() + labs(x="Focal (column sums to 1)", y="Neighbour",
    title="BC 313: contamination-cleaned neighbour attribution (RWA, spec>=0.5 genes only)",
    subtitle="decomposes only genes the focal type owns -> removes bleed-driven spatial variance") +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("streaming/figures/neighbour_share_rwa_cleaned_bc.pdf", p, width=8.2, height=7, device=cairo_pdf)
say("\nwrote streaming/figures/neighbour_share_rwa_cleaned_bc.pdf")

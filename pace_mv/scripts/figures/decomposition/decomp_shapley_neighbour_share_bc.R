## decomp_shapley_neighbour_share_bc.R
## Column-normalized neighbour attribution of focal spatial variance, via the SHAPLEY / LMG
## relative-importance decomposition (Lindeman-Merenda-Gold 1980; Gromping 2006, relaimpo).
## Per focal F, decompose V_state(F) = sum_g b_g' Sigma b_g among the HETEROTYPIC neighbour
## predictors; Shapley value of neighbour k = average marginal contribution over all orderings.
## Non-negative, sums to total -> columns sum to 1. Fairly splits variance shared by collinear
## (e.g. abundant) neighbours, so abundance does not double-count.
## Closed form used: sum over genes of [variance of eta explained by predictor-subset S]
##   = trace( Sigma_SS^{-1} G_SS ),  G = Sigma B Sigma,  B = sum_g (b_g b_g') + diag(se^2).
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(tidyr); library(dplyr)})
say <- function(...) cat(sprintf(...),"\n")
mv  <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
fit <- mv$fit; U <- fit$U; seU <- fit$se_U; Z <- fit$re_meta$Z
TYPES <- sort(unique(as.character(mv$cell_meta$celltype)))

## value function: variance explained by predictor subset S (gene-summed)
val_fun <- function(S, Sig, G) if (length(S)==0) 0 else
  sum(diag(solve(Sig[S,S,drop=FALSE], G[S,S,drop=FALSE])))

## exact Shapley over K predictors (K<=8 here -> cheap), with memoised value function
shapley <- function(Sig, G) {
  K <- nrow(Sig); ids <- seq_len(K); phi <- numeric(K)
  cache <- new.env(hash=TRUE)
  v <- function(S){ if (length(S)==0) return(0); key <- paste(S,collapse=","); got <- cache[[key]]
    if (!is.null(got)) return(got); out <- val_fun(S,Sig,G); cache[[key]] <- out; out }
  for (k in ids) {
    others <- setdiff(ids,k)
    for (m in 0:length(others)) {
      wm <- factorial(m)*factorial(K-m-1)/factorial(K)
      cb <- if (m==0) matrix(integer(0),nrow=0,ncol=1) else combn(others,m)
      for (ci in seq_len(ncol(cb))) { S <- cb[,ci]; phi[k] <- phi[k] + wm*(v(c(S,k))-v(S)) }
    }
  }
  setNames(phi, rownames(Sig))
}

share <- matrix(0,length(TYPES),length(TYPES),dimnames=list(neighbour=TYPES,focal=TYPES))
for (F in TYPES) {
  colI <- paste0(F,"::(Intercept)"); if(!(colI %in% colnames(Z))) next
  cF <- which(Z[,colI]!=0)
  nb_cols <- paste0(F,"::",setdiff(TYPES,F)); nb_cols <- nb_cols[nb_cols %in% colnames(Z)]
  N <- as.matrix(Z[cF, nb_cols, drop=FALSE])
  keep <- which(apply(N,2,var) > 1e-12)          # drop zero-variance (dropped) pairs
  if (length(keep) < 2) next
  N <- N[,keep,drop=FALSE]; cols <- nb_cols[keep]
  nbr <- sub(paste0("^",F,"::"),"",cols)
  Sig <- stats::cov(N)
  Ucn <- U[cols,,drop=FALSE]; Scn <- seU[cols,,drop=FALSE]
  B   <- Ucn %*% t(Ucn) + diag(rowSums(Scn^2), nrow=length(cols))   # K x K slope cross-product
  G   <- Sig %*% B %*% Sig
  phi <- pmax(shapley(Sig, G), 0)                 # LMG values are >=0 up to numerical noise
  share[nbr, F] <- phi / sum(phi)                 # column-normalize -> sums to 1
}
saveRDS(share, "streaming/neighbour_share_shapley_bc.rds")

say("=== Shapley neighbour share per focal (column sums to 1, heterotypic) — top 3 ===")
for (F in TYPES) {
  v <- sort(share[,F], decreasing=TRUE); v <- v[v>0][1:3]
  say("  %-15s : %s", F, paste(sprintf("%s %.0f%%", names(v), 100*v), collapse=", "))
}

d <- as.data.frame(as.table(share)) |> rename(s=Freq) |>
  mutate(focal=factor(focal,levels=TYPES), neighbour=factor(neighbour,levels=rev(TYPES)),
         s=ifelse(as.character(focal)==as.character(neighbour),NA,s))
p <- ggplot(d, aes(focal, neighbour, fill=s)) +
  geom_tile(colour="grey92", linewidth=0.3) +
  geom_text(aes(label=ifelse(is.na(s)|s<0.02,"",sprintf("%.0f",100*s))), size=3) +
  scale_fill_viridis_c(option="magma", na.value="white", limits=c(0,1), labels=scales::percent,
                       name="Shapley\nshare") +
  coord_equal() +
  labs(x="Focal (column sums to 1)", y="Neighbour",
       title="BC 313-panel: Shapley/LMG attribution of focal spatial variance",
       subtitle="share of each focal's heterotypic spatial variance explained by each neighbour (fair allocation across collinear neighbours)") +
  theme_classic(base_size=11) + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("streaming/figures/neighbour_share_shapley_bc.pdf", p, width=8.2, height=7, device=cairo_pdf)
cat("wrote streaming/figures/neighbour_share_shapley_bc.pdf\n")

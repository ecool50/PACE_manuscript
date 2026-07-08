## decomp_rwa_neighbour_share_bc.R
## Relative Weights Analysis (Genizi 1993; Johnson 2000) — the O(K^3), non-negative,
## sums-to-total approximation to LMG/Shapley. Decomposes focal spatial variance
## var(eta)=u'Sigma u among neighbour predictors via the symmetric sqrt of the predictor
## CORRELATION (orthonormal counterpart). Closed form, gene-vectorized; scales to many types.
##   D = diag(Sigma); R = D^-.5 Sigma D^-.5; Rh = R^.5;  B = sum_g u_g u_g' + diag(se^2)
##   G* = D^.5 B D^.5;  weight_k = sum_j (Rh_kj)^2 * diag(Rh G* Rh)_j  ; share = weight/sum
## Decomposes the SAME total trace(B Sigma) as the Shapley version -> directly comparable.
suppressPackageStartupMessages({library(Matrix); library(dplyr)})
say <- function(...) cat(sprintf(...),"\n")
mv  <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
fit <- mv$fit; U <- fit$U; seU <- fit$se_U; Z <- fit$re_meta$Z
TYPES <- sort(unique(as.character(mv$cell_meta$celltype)))
shap <- readRDS("streaming/neighbour_share_shapley_bc.rds")   # neighbour x focal (Shapley)

sym_sqrt <- function(M){ e <- eigen(M, symmetric=TRUE); V <- e$vectors
  V %*% diag(sqrt(pmax(e$values,0)), nrow=ncol(V)) %*% t(V) }

share <- matrix(0,length(TYPES),length(TYPES),dimnames=list(neighbour=TYPES,focal=TYPES))
t0 <- Sys.time()
for (F in TYPES) {
  colI <- paste0(F,"::(Intercept)"); if(!(colI %in% colnames(Z))) next
  cF <- which(Z[,colI]!=0)
  nb_cols <- paste0(F,"::",setdiff(TYPES,F)); nb_cols <- nb_cols[nb_cols %in% colnames(Z)]
  N <- as.matrix(Z[cF, nb_cols, drop=FALSE]); keep <- which(apply(N,2,var)>1e-12)
  if(length(keep)<2) next
  N <- N[,keep,drop=FALSE]; cols <- nb_cols[keep]; nbr <- sub(paste0("^",F,"::"),"",cols)
  Sig <- stats::cov(N); D <- diag(Sig); Dh <- sqrt(D)
  R  <- Sig / outer(Dh,Dh)                              # correlation
  Rh <- sym_sqrt(R)
  Ucn <- U[cols,,drop=FALSE]; Scn <- seU[cols,,drop=FALSE]
  B   <- Ucn %*% t(Ucn) + diag(rowSums(Scn^2), nrow=length(cols))   # K x K slope cross-product
  Gst <- (Dh * B) * rep(Dh, each=length(Dh))            # D^.5 B D^.5
  proxy_var <- diag(Rh %*% Gst %*% Rh)                  # variance on each orthonormal proxy
  w <- as.numeric((Rh*Rh) %*% proxy_var)                # relative weight per neighbour
  share[nbr, F] <- pmax(w,0) / sum(pmax(w,0))
}
say("RWA over %d focals computed in %.2fs (O(K^3), scales to many types)", length(TYPES), as.numeric(Sys.time()-t0,units="secs"))
saveRDS(share, "streaming/neighbour_share_rwa_bc.rds")

## ---- agreement with Shapley ----
say("\n=== RWA vs Shapley agreement (per focal, heterotypic shares) ===")
say("%-15s %10s %10s | %-28s %-28s", "focal","cosine","spearman","RWA top-2","Shapley top-2")
cos_all <- c()
for (F in TYPES) {
  a <- share[,F]; b <- shap[,F]; ok <- (a+b)>0
  if(sum(ok)<2) next
  cs <- sum(a[ok]*b[ok])/sqrt(sum(a[ok]^2)*sum(b[ok]^2))
  sp <- suppressWarnings(cor(a[ok],b[ok],method="spearman")); cos_all <- c(cos_all,cs)
  t2 <- function(x){v<-sort(x,decreasing=TRUE);v<-v[v>0][1:2];paste(sprintf("%s %.0f%%",names(v),100*v),collapse=", ")}
  say("%-15s %10.3f %10.2f | %-28s %-28s", F, cs, sp, t2(a), t2(b))
}
say("\nmean cosine(RWA, Shapley) across focals = %.3f", mean(cos_all))

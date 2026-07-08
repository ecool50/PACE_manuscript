## pace_pratt_macro.R — Pratt decomposition levels 2 (feature/source-specific spillover)
## and 3 (global macro-block) for PACE-MV. Closed-form signed Pratt: the variance
## contribution of a part = Cov(part, whole), so contributions sum exactly to Var(whole).
##
## LEVEL 3 (global macro-block): decompose the modelled STRUCTURED link-scale predictor
## eta_hat_i = u0_i + spatial_i + spill_i across ALL cells. NOTE: eta_hat excludes the global
## intercept (beta0_g), the library/offset term, AND the observation-level residual entirely, so
## these are shares OF THE MODELLED STRUCTURE, not of a gene's total expression variance (audited
## 2026-06-20). The ~93% "cell type" is therefore "of explained structure", not "of all variance".
##   V_celltype = Cov_global(u0,   eta_hat)   [u0_i  = cell-type baseline BLUP for cell i]
##   V_spatial  = Cov_global(spat, eta_hat)   [spat_i = sum_t u^t_c N_it]
##   V_spill    = Cov_global(spill,eta_hat)   [spill_i = log(mu_i) - log(mu_bio_i), link-scale]
## They sum to Var(eta_hat). This explicitly absorbs cell-type<->spatial colocalization
## (e.g. T cells near tumour) via the covariances, unlike the law-of-total-variance four-block.
##
## LEVEL 2 (source-specific spillover): split the per-cell spillover S_i = sum_g mu_spill_ig
## by the SOURCE cell type of the ambient. Because mu_spill = (rho_i/af_i) * E^tech and E^tech is
## additive over source types, the per-cell scalar rho_i/af_i cancels in the SHARES, so
##   share^t_i = ET^t_i / sum_t ET^t_i,   ET^t_i = sum_{j in nbr(i), type(j)=t != type(i)} exp(-d/h)*totcount_j
##   S^t_i = S_i * share^t_i  (sums to S_i exactly)
## Pratt within the spillover block, per focal type: V_spill^t = Cov_focal(S^t, S).

## column covariance contributions (Pratt) for one whole vector
.pratt_cov <- function(part, whole) {
  stats::cov(part, whole)            # = Cov(part, whole); sum over parts = Var(whole)
}

## ---- LEVEL 3: global macro-block Pratt, per gene + pooled ----
pratt_macro_global <- function(fit, mu, mu_bio, genes, headline = NULL) {
  Z <- fit$re_meta$Z; U <- fit$U
  if (is.null(colnames(U))) colnames(U) <- genes
  int_cols <- grep("::\\(Intercept\\)$", colnames(Z))
  nb_cols  <- setdiff(seq_len(ncol(Z)), int_cols)
  Z_int <- Z[, int_cols, drop = FALSE]; U_int <- U[int_cols, , drop = FALSE]
  Z_nb  <- Z[, nb_cols,  drop = FALSE]; U_nb  <- U[nb_cols,  , drop = FALSE]

  C_ct <- C_sp <- C_bl <- V_tot <- numeric(length(genes))
  per_gene <- vector("list", length(genes))
  for (gi in seq_along(genes)) {
    g <- genes[gi]
    u0   <- as.numeric(Z_int %*% U_int[, gi])
    spat <- as.numeric(Z_nb  %*% U_nb[, gi])
    spill<- log(pmax(mu[, gi], 1e-12)) - log(pmax(mu_bio[, gi], 1e-12))
    eta  <- u0 + spat + spill
    vt <- stats::var(eta)
    cc <- .pratt_cov(u0, eta); cs <- .pratt_cov(spat, eta); cb <- .pratt_cov(spill, eta)
    C_ct[gi] <- cc; C_sp[gi] <- cs; C_bl[gi] <- cb; V_tot[gi] <- vt
    per_gene[[gi]] <- data.frame(gene = g,
      cell_type_pct = 100*cc/vt, spatial_pct = 100*cs/vt, spillover_pct = 100*cb/vt)
  }
  pooled <- data.frame(
    cell_type_pct = 100*sum(C_ct)/sum(V_tot),
    spatial_pct   = 100*sum(C_sp)/sum(V_tot),
    spillover_pct = 100*sum(C_bl)/sum(V_tot))
  list(pooled = pooled, per_gene = do.call(rbind, per_gene),
       headline = if (!is.null(headline)) do.call(rbind, per_gene)[match(headline, genes), ] else NULL)
}

## ---- ambient mass per source cell type: ET^t_i (n x nTypes), per-slide frNN ----
build_source_ambient <- function(df, Y, h_tech = 5, rad_mult = 3) {
  TYPES <- sort(unique(as.character(df$celltype)))
  ct <- as.character(df$celltype); n <- nrow(df)
  totcount <- as.numeric(Matrix::rowSums(Y))
  rad <- rad_mult * h_tech
  ET <- matrix(0, n, length(TYPES), dimnames = list(NULL, TYPES))
  by_slide <- split(seq_len(n), df$imageID)
  for (si in seq_along(by_slide)) {
    rows <- by_slide[[si]]; if (length(rows) < 2) next
    nn <- dbscan::frNN(as.matrix(df[rows, c("x","y")]), eps = rad)
    ct_im <- ct[rows]; tot_im <- totcount[rows]
    for (i in seq_along(rows)) {
      nb <- nn$id[[i]]; if (!length(nb)) next
      keep <- ct_im[nb] != ct_im[i]; if (!any(keep)) next
      w <- exp(-nn$dist[[i]][keep] / h_tech)
      contrib <- w * tot_im[nb[keep]]
      tt <- ct_im[nb[keep]]
      ET[rows[i], ] <- ET[rows[i], ] + tapply(contrib, factor(tt, levels = TYPES), sum)[TYPES] |>
        (\(x){ x[is.na(x)] <- 0; x })()
    }
  }
  ET
}

## ---- LEVEL 2: source-specific spillover Pratt, per focal type ----
pratt_spillover_by_source <- function(fit, mu_spill, df, ET, focals = NULL) {
  TYPES <- colnames(ET); ct <- as.character(df$celltype)
  if (is.null(focals)) focals <- TYPES
  S_all <- as.numeric(Matrix::rowSums(mu_spill))     # total per-cell spillover mass (from fit)
  ET_tot <- rowSums(ET)
  out <- list()
  for (fc in focals) {
    idx <- which(ct == fc); idx <- idx[ET_tot[idx] > 0 & S_all[idx] > 0]
    if (length(idx) < 50) next
    S <- S_all[idx]
    share <- ET[idx, , drop = FALSE] / ET_tot[idx]   # source-type shares (edge & rho cancel)
    Vtot <- stats::var(S)
    rows <- lapply(TYPES, function(t) {
      St <- S * share[, t]
      data.frame(focal = fc, source = t, V_pratt = stats::cov(St, S),
                 share_pct = 100*stats::cov(St, S)/Vtot)
    })
    out[[fc]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

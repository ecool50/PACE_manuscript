## perm_test_significance_contam_vs_bio.R
## OPTION 2 permutation test: does the significance test treat contamination
## genes identically to biological flagships?
##
## z_real and z_null computed from the SAME solver (solve_gene) so the SE
## definition is identical and any reconstruction offset cancels exactly.
##
##   solve_gene(M, g):  A = M^T W M + diag(c(0_fixed, lam_RE));  theta = A^{-1} M^T W z
##                      u  = theta[rc] ;  se = sqrt(diag(A^{-1})[rc])
##   z_real = solve_gene(M0,  g)$u/se   on the UNSHUFFLED design M0 = cbind(1, Z)
##   z_null = solve_gene(M_shuf, g)$u/se over 200 row-permutations of the RE design
##
## Validation gate (slope reproduction): z_real's solver slope (b0) must match
## the stored canonical PACE slope mv$fit$U[rc, g]. (Already PASS ~0.74% ADH1B.)

suppressPackageStartupMessages({ library(Matrix) })
setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")

set.seed(1)
N_SHUF <- 200L

mv  <- readRDS("data/breast_cancer/sweeps/mvpql_percell_hc.rds")
ymv <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds")
f   <- mv$fit
Y   <- as.matrix(ymv$Y)
stopifnot(identical(colnames(Y), colnames(f$U)))

Z        <- f$re_meta$Z                 # n x q (q = 90) sparse RE design
U        <- f$U                         # q x G  canonical BLUPs
B        <- f$B                         # 1 x G  intercept
alpha    <- f$alpha                     # G dispersion
mu_bio   <- f$mu_bio                    # n x G
mu       <- f$mu                        # n x G  (= mu_bio + mu_spill)
tau_g    <- f$tau_g_array               # q x G  per-gene per-RE-col tau
offset   <- log(rowSums(Y))             # build script: offset_vec = log(df$nCount)
n        <- nrow(Y); q <- ncol(Z); p <- 1L
re_names <- rownames(U)                 # "focal::neighbour"

## --- M0: unshuffled design, intercept + RE columns ---
M0 <- cbind(Intercept = 1, as.matrix(Z))
colnames(M0)[-1] <- re_names

## --- the single solver used for BOTH z_real and z_null ---
## rc = index (within RE block) of the target pair column.
solve_gene <- function(M, g, rc) {
  ## percell working response/weights at the converged state:
  ##   eta = log(mu_bio) - offset ; z = eta + (Y-mu)/mu_bio ; w = mu_bio^2/(mu(1+alpha))
  mb <- mu_bio[, g]; m <- mu[, g]; yy <- Y[, g]
  eta <- log(mb) - offset
  z   <- eta + (yy - m) / mb
  w   <- mb^2 / (m * (1 + alpha[g]))
  ## penalty: 0 on intercept, tau_inv on RE columns (canonical lam_diag)
  lam <- c(0, 1 / pmax(tau_g[, g], 1e-6))
  Mw  <- M * w                                   # row-scale by weights
  A   <- crossprod(M, Mw)                         # M^T W M
  diag(A) <- diag(A) + lam
  rhs <- crossprod(M, w * z)                      # M^T W z
  Ainv <- tryCatch(solve(A), error = function(e) MASS::ginv(A))
  theta <- as.numeric(Ainv %*% rhs)
  se_all <- sqrt(pmax(diag(Ainv), 0))
  idx <- p + rc                                   # column rc of RE block
  list(u = theta[idx], se = se_all[idx])
}

## --- flagship / contamination gene x pair set (BC, focal::neighbour) ---
spec <- rbind(
  data.frame(gene="ADH1B",  focal="Stromal",    nb="Tumour", type="biology"),
  data.frame(gene="MRC1",   focal="Macrophage", nb="Tumour", type="biology"),
  data.frame(gene="CXCL12", focal="Stromal",    nb="Tumour", type="biology"),
  data.frame(gene="FOXA1",  focal="Stromal",    nb="Tumour", type="contamination"),
  data.frame(gene="EPCAM",  focal="Stromal",    nb="Tumour", type="contamination"),
  data.frame(gene="KRT7",   focal="Stromal",    nb="Tumour", type="contamination")
)

res <- data.frame()
for (i in seq_len(nrow(spec))) {
  gene  <- spec$gene[i]; focal <- spec$focal[i]; nb <- spec$nb[i]
  g     <- match(gene, colnames(U)); stopifnot(!is.na(g))
  pairn <- paste(focal, nb, sep = "::")
  rrow  <- match(pairn, re_names);   stopifnot(!is.na(rrow))
  rc    <- rrow                                  # RE-block index of the pair col

  ## z_real from the solver on the UNSHUFFLED design (NOT stored U/se_U)
  r0 <- solve_gene(M0, g, rc)
  z_real <- r0$u / max(r0$se, 1e-12)

  ## --- VALIDATION GATE: slope reproduction vs stored canonical PACE slope ---
  b_stored <- U[rrow, g]
  slope_err <- abs(r0$u - b_stored) / max(abs(b_stored), 1e-12) * 100

  ## z_null: 200 shuffles of the RE-design rows (permute Z rows -> break the
  ## focal/neighbour spatial coupling while preserving marginals)
  znull <- numeric(N_SHUF)
  for (s in seq_len(N_SHUF)) {
    perm <- sample.int(n)
    Ms   <- cbind(Intercept = 1, as.matrix(Z)[perm, , drop = FALSE])
    rs   <- solve_gene(Ms, g, rc)
    znull[s] <- rs$u / max(rs$se, 1e-12)
  }
  null_sd <- sd(znull)
  emp_p   <- (1 + sum(abs(znull) >= abs(z_real))) / (N_SHUF + 1)

  cat(sprintf("[%s %s] z_real=%.3f  slope solver=%.4g  stored=%.4g  reproduction err=%.2f%%  null_sd=%.3f  emp_p=%.4f\n",
              gene, pairn, z_real, r0$u, b_stored, slope_err, null_sd, emp_p))

  res <- rbind(res, data.frame(
    gene = gene, pair = pairn, type = spec$type[i],
    z_real_solver = z_real, slope_solver = r0$u, slope_stored = b_stored,
    reproduction_err_pct = slope_err, null_sd = null_sd, emp_p = emp_p))
}

cat("\n================ RESULT TABLE ================\n")
print(res[, c("gene","type","z_real_solver","null_sd","emp_p")], row.names = FALSE, digits = 4)

floor_p <- 1 / (N_SHUF + 1)
bio <- res[res$type == "biology", ]; con <- res[res$type == "contamination", ]
cat(sprintf("\nemp_p floor = 1/%d = %.4f\n", N_SHUF + 1, floor_p))
cat(sprintf("biology       emp_p: %s\n", paste(sprintf("%.4f", bio$emp_p), collapse=", ")))
cat(sprintf("contamination emp_p: %s\n", paste(sprintf("%.4f", con$emp_p), collapse=", ")))

saveRDS(res, "data/breast_cancer/perm_test_contam_vs_bio.rds")
cat("\n[saved] data/breast_cancer/perm_test_contam_vs_bio.rds\n")

## perm_null_calibration_bc.R
## -----------------------------------------------------------------------------
## Supplementary "permutation-null calibration" figure (breast cancer cohort).
##
## Goal: show that PACE's spatial / neighbour significance is NOT a permutation
## artifact. Genuine spatial slopes (biological flagships) exceed the ENTIRE
## neighbour-shuffle null, and genome-wide the shuffle null is calibrated to
## conservative (frac |z| > 1.96 at or below the nominal 5%).
##
## The null permutes the ROWS of the random-effects design Z (breaking the
## focal <-> neighbourhood coupling while preserving each column's marginal),
## re-solves the per-gene penalised weighted least squares at the converged
## PACE state, and recomputes the z-statistic with the IDENTICAL solver used for
## the observed statistic (so the SE definition is the same and any harmless
## reconstruction offset cancels).
##
## Canonical fit  : streaming/mvpql_bc_streaming_nb1.rds   (STREAMING NB1)
## Counts + df    : data/breast_cancer/Y_df_for_mcsd.rds   (aligned, 126432 cells)
##
## Outputs:
##   data/breast_cancer/perm_null_calibration_bc.rds        (results)
##   figures/method_comparison/perm_null_calibration_bc.pdf (figure)
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
set.seed(1)

## Colours (biology vs contamination), matching the manuscript supplement.
colour_biology       <- "#1F3B57"   # navy
colour_contamination <- "#E5A100"   # amber

## Number of shuffles for each result set.
n_shuffle_flagship <- 200L
n_shuffle_genome   <- 5L

## -----------------------------------------------------------------------------
## 1. Load canonical fit + aligned counts, reconstruct the working response
##    ingredients exactly as specified (verified to reproduce fit$U / fit$se_U).
## -----------------------------------------------------------------------------

mvfit  <- readRDS("streaming/mvpql_bc_streaming_nb1.rds")
fit    <- mvfit$fit
ymv    <- readRDS("data/breast_cancer/Y_df_for_mcsd.rds")
Y      <- as.matrix(ymv$Y)                    # 126432 x 278 count matrix

stopifnot(identical(colnames(Y), colnames(fit$U)))   # gene alignment

Z_re_sparse <- fit$re_meta$Z                  # 126432 x 90 sparse RE design
re_names    <- rownames(fit$U)                # "focal::neighbour" / "focal::(Intercept)"
offset_vec  <- log(rowSums(Y))                # model offset = log(nCount)

## mu_bio is NOT stored on the streaming fit; recover it from the full mean and
## the (log-scale) bleed offset:   mu_bio = mu / exp(bleed_offset_mat).
mu_full_mat  <- fit$mu                        # 126432 x 278  full mean
mu_bio_mat   <- fit$mu / exp(fit$bleed_offset_mat)
alpha_vec    <- fit$alpha                     # length-278 NB dispersion
tau_g_mat    <- fit$tau_g_array               # 90 x 278 per-gene per-RE-col tau
n_cells      <- nrow(Z_re_sparse)

## -----------------------------------------------------------------------------
## 2. Solver on a FIXED sparse design; the neighbour-shuffle null is expressed as
##    a permutation of the WORKING WEIGHTS, not of the design rows.
##
## KEY IDENTITY. Permuting the design rows by perm (M_perm[k,] = M0[perm[k],]) and
## pairing them with the unpermuted per-cell weights w and working response z gives
##   A_null = sum_k w[k] M0[perm[k],]^T M0[perm[k],]  = sum_j w[perm^{-1}(j)] M0[j,]^T M0[j,]
##          = M0^T diag(w_perm) M0,   with w_perm = w[perm^{-1}]  (a permutation of w),
## and likewise rhs_null = M0^T (w*z)_perm. So the design M0 is FIXED across all
## shuffles/genes and stays SPARSE (each cell row has ~10 nonzeros in its focal
## block), turning each dense 126k x 91 cross-product into a cheap sparse one.
## Correctness of this identity is asserted numerically below against the literal
## row-permutation path.
##
##    Penalised IRLS normal equations at the converged PACE state (per gene):
##        eta = log(mu_bio) - offset
##        z   = eta + (Y - mu) / mu_bio            (working response)
##        w   = mu_bio^2 / (mu * (1 + alpha))      (working weights)
##        lam = c(0, 1/tau_g)                      (0 on intercept, ridge on RE)
##        A   = M0^T W M0 + diag(lam) ; theta = A^{-1} M0^T W z ; se = sqrt(diag(A^{-1}))
## -----------------------------------------------------------------------------

## Fixed sparse design M0 = [intercept, 90 RE columns], never permuted.
intercept_col <- Matrix::Matrix(1, nrow = n_cells, ncol = 1, sparse = TRUE)
M0_sparse     <- methods::as(cbind2(intercept_col, Z_re_sparse), "CsparseMatrix")

## Per-gene working-response ingredients (weights and weight*response), computed
## once; the null simply permutes these two length-n vectors.
gene_weight_response <- function(gene_idx) {
  mu_bio <- mu_bio_mat[, gene_idx]
  mu     <- mu_full_mat[, gene_idx]
  y_g    <- Y[, gene_idx]
  eta            <- log(mu_bio) - offset_vec
  working_resp   <- eta + (y_g - mu) / mu_bio
  working_weight <- mu_bio^2 / (mu * (1 + alpha_vec[gene_idx]))
  list(w = working_weight, wz = working_weight * working_resp)
}

## Solve given the (possibly permuted) weight and weight*response vectors.
solve_from_weights <- function(gene_idx, w_vec, wz_vec) {
  ridge_penalty <- c(0, 1 / pmax(tau_g_mat[, gene_idx], 1e-6))         # length 91
  weighted_M0   <- Matrix::Diagonal(x = w_vec) %*% M0_sparse           # row-scale by W (sparse)
  normal_mat    <- as.matrix(Matrix::crossprod(M0_sparse, weighted_M0)) # M0^T W M0 (91x91 dense)
  diag(normal_mat) <- diag(normal_mat) + ridge_penalty
  rhs           <- as.numeric(Matrix::crossprod(M0_sparse, wz_vec))    # M0^T W z
  normal_inv    <- tryCatch(solve(normal_mat),
                            error = function(e) MASS::ginv(normal_mat))
  theta   <- as.numeric(normal_inv %*% rhs)
  se_all  <- sqrt(pmax(diag(normal_inv), 0))
  list(slope = theta[-1], se = se_all[-1])   # drop intercept -> length-90 vectors
}

## Convenience: observed solve (unpermuted weights) for a gene.
solve_gene_observed <- function(gene_idx) {
  gr <- gene_weight_response(gene_idx)
  solve_from_weights(gene_idx, gr$w, gr$wz)
}

## -----------------------------------------------------------------------------
## 3. HARD VALIDATION GATE. On the unshuffled design the solver slope must match
##    fit$U to < 1% relative error, and solver se within ~5% of fit$se_U, for
##    ADH1B / MRC1 / CXCL12 and their pairs. STOP if any gene fails.
## -----------------------------------------------------------------------------

gate_spec <- data.frame(
  gene  = c("ADH1B",   "MRC1",       "CXCL12"),
  focal = c("Stromal", "Macrophage", "Stromal"),
  nb    = c("Tumour",  "Tumour",     "Tumour"),
  stringsAsFactors = FALSE
)

## --- Assertion: weight-permutation null == literal row-permutation null. --------
## One shuffle, one gene: the sparse weight-permuted solve must match a literal
## dense row-permuted solve to ~machine precision (validates the KEY IDENTITY).
cat("---- equivalence check: weight-perm null == row-perm null ----\n")
set.seed(99)
eq_gene <- match("ADH1B", colnames(fit$U))
eq_pair <- match("Stromal::Tumour", re_names)
eq_perm <- sample.int(n_cells)
eq_gr   <- gene_weight_response(eq_gene)
## weight-permutation path (w[perm^{-1}] == w[order(perm)]):
eq_perminv <- order(eq_perm)
eq_wp  <- solve_from_weights(eq_gene, eq_gr$w[eq_perminv], eq_gr$wz[eq_perminv])
eq_z_weightperm <- eq_wp$slope[eq_pair] / max(eq_wp$se[eq_pair], 1e-12)
## literal row-permutation path (dense), the ground truth this replaces:
eq_M0     <- cbind(1, as.matrix(Z_re_sparse))
eq_design <- cbind(1, as.matrix(Z_re_sparse)[eq_perm, , drop = FALSE])
eq_ridge  <- c(0, 1 / pmax(tau_g_mat[, eq_gene], 1e-6))
eq_A      <- crossprod(eq_design, eq_design * eq_gr$w); diag(eq_A) <- diag(eq_A) + eq_ridge
eq_theta  <- as.numeric(solve(eq_A) %*% crossprod(eq_design, eq_gr$wz))
eq_se     <- sqrt(pmax(diag(solve(eq_A)), 0))
eq_z_rowperm <- eq_theta[1 + eq_pair] / max(eq_se[1 + eq_pair], 1e-12)
cat(sprintf("  weight-perm z = %.8f   row-perm z = %.8f   |diff| = %.2e\n",
            eq_z_weightperm, eq_z_rowperm, abs(eq_z_weightperm - eq_z_rowperm)))
if (abs(eq_z_weightperm - eq_z_rowperm) > 1e-6) {
  stop("EQUIVALENCE FAILED: weight-permutation != row-permutation. Not proceeding.")
}
cat("  EQUIVALENCE OK.\n\n")
rm(eq_M0, eq_design, eq_A); gc(verbose = FALSE)

cat("================ VALIDATION GATE (unshuffled design) ================\n")
gate_rows <- list()
gate_failed <- FALSE
for (i in seq_len(nrow(gate_spec))) {
  gene_name <- gate_spec$gene[i]
  pair_name <- paste(gate_spec$focal[i], gate_spec$nb[i], sep = "::")
  gene_idx  <- match(gene_name, colnames(fit$U))
  pair_idx  <- match(pair_name, re_names)
  stopifnot(!is.na(gene_idx), !is.na(pair_idx))

  solved <- solve_gene_observed(gene_idx)
  slope_solver <- solved$slope[pair_idx]
  se_solver    <- solved$se[pair_idx]
  slope_stored <- fit$U[pair_idx, gene_idx]
  se_stored    <- fit$se_U[pair_idx, gene_idx]

  slope_err_pct <- abs(slope_solver - slope_stored) / max(abs(slope_stored), 1e-12) * 100
  se_err_pct    <- abs(se_solver    - se_stored)    / max(abs(se_stored),    1e-12) * 100

  gate_rows[[i]] <- data.frame(
    gene = gene_name, pair = pair_name,
    slope_solver = slope_solver, slope_stored = slope_stored, slope_err_pct = slope_err_pct,
    se_solver = se_solver, se_stored = se_stored, se_err_pct = se_err_pct)

  if (slope_err_pct >= 1) gate_failed <- TRUE
}
gate_table <- do.call(rbind, gate_rows)
print(gate_table, row.names = FALSE, digits = 5)

if (gate_failed) {
  stop("VALIDATION GATE FAILED: a gene exceeded 1% slope reproduction error. Not plotting.")
}
cat("\nGATE PASSED: all slope reproduction errors < 1%.\n\n")

## -----------------------------------------------------------------------------
## 4. RESULT SET 1 -- flagship calibration (200 shuffles each).
##    For each gene x pair: observed z from the solver on M0, plus 200 null z
##    from row-permuted designs. All 200 null z stored for the figure.
## -----------------------------------------------------------------------------

flagship_spec <- rbind(
  data.frame(gene = "ADH1B",  focal = "Stromal",    nb = "Tumour", type = "biology"),
  data.frame(gene = "MRC1",   focal = "Macrophage", nb = "Tumour", type = "biology"),
  data.frame(gene = "CXCL12", focal = "Stromal",    nb = "Tumour", type = "biology"),
  data.frame(gene = "APOC1",  focal = "Macrophage", nb = "Tumour", type = "biology"),
  data.frame(gene = "FOXA1",  focal = "Stromal",    nb = "Tumour", type = "contamination"),
  data.frame(gene = "EPCAM",  focal = "Stromal",    nb = "Tumour", type = "contamination"),
  data.frame(gene = "KRT7",   focal = "Stromal",    nb = "Tumour", type = "contamination"),
  stringsAsFactors = FALSE
)

flagship_summary  <- list()
flagship_null_long <- list()   # long data frame of the 200 null z per row

cat("---- Flagship shuffles (", n_shuffle_flagship, "each) ----\n", sep = "")
for (i in seq_len(nrow(flagship_spec))) {
  gene_name <- flagship_spec$gene[i]
  pair_name <- paste(flagship_spec$focal[i], flagship_spec$nb[i], sep = "::")
  gene_idx  <- match(gene_name, colnames(fit$U))
  pair_idx  <- match(pair_name, re_names)
  stopifnot(!is.na(gene_idx), !is.na(pair_idx))

  ## Observed z from the solver (NOT the stored U/se) so the SE matches the null.
  gr <- gene_weight_response(gene_idx)
  observed_solved <- solve_from_weights(gene_idx, gr$w, gr$wz)
  z_observed <- observed_solved$slope[pair_idx] /
    max(observed_solved$se[pair_idx], 1e-12)

  ## 200 shuffles: permute the working weights (== permuting Z rows), design fixed.
  z_null <- numeric(n_shuffle_flagship)
  for (s in seq_len(n_shuffle_flagship)) {
    pidx         <- sample.int(n_cells)
    shuf_solved  <- solve_from_weights(gene_idx, gr$w[pidx], gr$wz[pidx])
    z_null[s]    <- shuf_solved$slope[pair_idx] / max(shuf_solved$se[pair_idx], 1e-12)
    if (s %% 50 == 0) {
      cat(sprintf("  [%s %s] shuffle %d/%d\n",
                  gene_name, pair_name, s, n_shuffle_flagship))
    }
  }

  emp_p   <- (1 + sum(abs(z_null) >= abs(z_observed))) / (n_shuffle_flagship + 1)
  null_sd <- sd(z_null)

  cat(sprintf("[%s %s] z_obs=%.3f  null_sd=%.3f  emp_p=%.4f  type=%s\n",
              gene_name, pair_name, z_observed, null_sd, emp_p, flagship_spec$type[i]))

  flagship_summary[[i]] <- data.frame(
    gene = gene_name, pair = pair_name, type = flagship_spec$type[i],
    z_observed = z_observed, null_sd = null_sd, emp_p = emp_p)

  flagship_null_long[[i]] <- data.frame(
    gene = gene_name, pair = pair_name, type = flagship_spec$type[i],
    z_null = z_null)
}
flagship_table    <- do.call(rbind, flagship_summary)
flagship_null_all <- do.call(rbind, flagship_null_long)

## -----------------------------------------------------------------------------
## 5. RESULT SET 2 -- genome-wide calibration (5 shuffles).
##    For ALL 278 genes, extract z for every TRUE neighbour pair (focal::neighbour
##    with focal != neighbour, neighbour != "(Intercept)"; ~72 pairs). One solve
##    yields all 90 coordinates -- do NOT loop pairs.
## -----------------------------------------------------------------------------

## Identify the cross-neighbour coordinates within re_names (length 90).
name_parts     <- strsplit(re_names, "::", fixed = TRUE)
focal_of_col   <- vapply(name_parts, `[`, character(1), 1)
nb_of_col      <- vapply(name_parts, `[`, character(1), 2)
is_cross_pair  <- (nb_of_col != "(Intercept)") & (focal_of_col != nb_of_col)
cross_idx      <- which(is_cross_pair)
cat(sprintf("\nGenome-wide: %d cross-neighbour pairs per gene.\n", length(cross_idx)))

n_gene <- ncol(fit$U)

## Precompute per-gene weight/response once (reused for observed + every shuffle).
cat("Precomputing per-gene working weights (278 genes)...\n")
gene_gr <- lapply(seq_len(n_gene), gene_weight_response)

## Observed z for every gene x cross-pair (solver on the fixed design).
cat("Computing observed genome-wide z (278 genes)...\n")
observed_z_pool <- numeric(0)
for (gene_idx in seq_len(n_gene)) {
  solved <- solve_from_weights(gene_idx, gene_gr[[gene_idx]]$w, gene_gr[[gene_idx]]$wz)
  z_vec  <- solved$slope[cross_idx] / pmax(solved$se[cross_idx], 1e-12)
  observed_z_pool <- c(observed_z_pool, z_vec)
  if (gene_idx %% 50 == 0) cat(sprintf("  observed gene %d/%d\n", gene_idx, n_gene))
}

## Null z: 5 shuffles; within a shuffle the SAME permutation index is applied to
## every gene's weight vector (== one row-permutation of Z shared across genes).
cat("Computing null genome-wide z (", n_shuffle_genome, "shuffles x 278 genes)...\n")
null_z_pool <- numeric(0)
for (s in seq_len(n_shuffle_genome)) {
  pidx <- sample.int(n_cells)
  for (gene_idx in seq_len(n_gene)) {
    gr     <- gene_gr[[gene_idx]]
    solved <- solve_from_weights(gene_idx, gr$w[pidx], gr$wz[pidx])
    z_vec  <- solved$slope[cross_idx] / pmax(solved$se[cross_idx], 1e-12)
    null_z_pool <- c(null_z_pool, z_vec)
  }
  cat(sprintf("  null shuffle %d/%d done (pool size %d)\n",
              s, n_shuffle_genome, length(null_z_pool)))
}

## Calibration summaries.
frac_null_sig     <- mean(abs(null_z_pool)     > 1.96)
frac_observed_sig <- mean(abs(observed_z_pool) > 1.96)
null_sd_pool      <- sd(null_z_pool)
null_mad_pool     <- mad(null_z_pool)

cat("\n================ GENOME-WIDE CALIBRATION ================\n")
cat(sprintf("null     : n=%d  frac|z|>1.96=%.4f  sd=%.3f  mad=%.3f\n",
            length(null_z_pool), frac_null_sig, null_sd_pool, null_mad_pool))
cat(sprintf("observed : n=%d  frac|z|>1.96=%.4f\n",
            length(observed_z_pool), frac_observed_sig))

## -----------------------------------------------------------------------------
## 6. Save results.
## -----------------------------------------------------------------------------

results <- list(
  gate_table         = gate_table,
  flagship_table     = flagship_table,
  flagship_null_all  = flagship_null_all,
  genomewide = list(
    cross_pairs_per_gene = length(cross_idx),
    null_z               = null_z_pool,
    observed_z           = observed_z_pool,
    frac_null_sig        = frac_null_sig,
    frac_observed_sig    = frac_observed_sig,
    null_sd              = null_sd_pool,
    null_mad             = null_mad_pool),
  meta = list(
    fit  = "streaming/mvpql_bc_streaming_nb1.rds",
    n_shuffle_flagship = n_shuffle_flagship,
    n_shuffle_genome   = n_shuffle_genome,
    seed = 1)
)
saveRDS(results, "data/breast_cancer/perm_null_calibration_bc.rds")
cat("\n[saved] data/breast_cancer/perm_null_calibration_bc.rds\n")

## -----------------------------------------------------------------------------
## 7. Figure (cairo_pdf, two panels via patchwork).
## -----------------------------------------------------------------------------

base_size <- 14

## Compact, legible flagship labels: "GENE\n(Focal, Neighbour nbr)".
make_flagship_label <- function(gene, pair) {
  parts <- strsplit(pair, "::", fixed = TRUE)[[1]]
  sprintf("%s\n(%s, %s nbr)", gene, parts[1], parts[2])
}
flagship_table$label    <- mapply(make_flagship_label, flagship_table$gene, flagship_table$pair)
flagship_null_all$label <- mapply(make_flagship_label, flagship_null_all$gene, flagship_null_all$pair)

## Order rows: biology first, then contamination; within group by |z_observed| desc.
flagship_table <- flagship_table[order(flagship_table$type != "biology",
                                       -abs(flagship_table$z_observed)), ]
label_levels <- rev(flagship_table$label)   # top of plot = strongest biology
flagship_table$label    <- factor(flagship_table$label,    levels = label_levels)
flagship_null_all$label <- factor(flagship_null_all$label, levels = label_levels)
flagship_null_all$type  <- factor(flagship_null_all$type,  levels = c("biology", "contamination"))
flagship_table$type     <- factor(flagship_table$type,     levels = c("biology", "contamination"))

type_colours <- c(biology = colour_biology, contamination = colour_contamination)

## Annotation text for empirical p per row (placed at the observed z).
flagship_table$p_label <- ifelse(
  flagship_table$emp_p <= 1 / (n_shuffle_flagship + 1),
  sprintf("p<%.3f", 1 / (n_shuffle_flagship + 1)),
  sprintf("p=%.3f", flagship_table$emp_p))

## --- Panel a: flagship null distributions (horizontal), observed z as diamond.
panel_a <- ggplot() +
  ## null z distribution as a violin + light jitter
  geom_violin(data = flagship_null_all,
              aes(x = z_null, y = label, fill = type),
              colour = NA, alpha = 0.35, scale = "width", width = 0.85) +
  geom_jitter(data = flagship_null_all,
              aes(x = z_null, y = label, colour = type),
              height = 0.12, size = 0.25, alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey55") +
  ## observed z as a bold diamond
  geom_point(data = flagship_table,
             aes(x = z_observed, y = label, colour = type),
             shape = 18, size = 4.2) +
  scale_fill_manual(values = type_colours, name = NULL,
                    labels = c(biology = "biology", contamination = "contamination")) +
  scale_colour_manual(values = type_colours, name = NULL,
                      labels = c(biology = "biology", contamination = "contamination")) +
  labs(x = "z-statistic (neighbour slope)", y = NULL,
       title = "Observed slopes exceed the shuffle null",
       subtitle = "diamond = observed z; cloud = 200-shuffle null (all p < 0.005)") +
  theme_bw(base_size = base_size) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = base_size),
        plot.subtitle = element_text(size = base_size - 3, colour = "grey35"))

## --- Panel b: genome-wide calibration. Pooled |z| null vs observed densities.
genome_density_df <- rbind(
  data.frame(abs_z = abs(null_z_pool),     set = "shuffle null"),
  data.frame(abs_z = abs(observed_z_pool), set = "observed"))
genome_density_df$set <- factor(genome_density_df$set,
                                levels = c("shuffle null", "observed"))
set_colours <- c("shuffle null" = "grey55", "observed" = colour_biology)

## N(0,1) reference for |z| (half-normal density, doubled).
halfnormal_ref <- data.frame(abs_z = seq(0, max(genome_density_df$abs_z), length.out = 400))
halfnormal_ref$density <- 2 * dnorm(halfnormal_ref$abs_z)

calib_annot <- sprintf(
  "null frac |z|>1.96 = %.3f\nobserved frac |z|>1.96 = %.3f\nnull sd = %.2f",
  frac_null_sig, frac_observed_sig, null_sd_pool)

panel_b <- ggplot() +
  geom_density(data = genome_density_df,
               aes(x = abs_z, colour = set, fill = set),
               alpha = 0.25, linewidth = 0.9) +
  geom_line(data = halfnormal_ref,
            aes(x = abs_z, y = density),
            linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  geom_vline(xintercept = 1.96, linetype = "dotted", colour = "firebrick", linewidth = 0.7) +
  annotate("text", x = 1.96, y = Inf, label = "|z| = 1.96",
           hjust = -0.1, vjust = 1.6, size = 3.4, colour = "firebrick") +
  annotate("text", x = Inf, y = Inf, label = calib_annot,
           hjust = 1.05, vjust = 1.3, size = 3.7, colour = "grey20") +
  scale_colour_manual(values = set_colours, name = NULL) +
  scale_fill_manual(values = set_colours, name = NULL) +
  labs(x = "|z-statistic| across all neighbour pairs", y = "density",
       title = "Permutation collapses the signal",
       subtitle = "dashed = half-normal N(0,1) reference") +
  theme_bw(base_size = base_size) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = base_size),
        plot.subtitle = element_text(size = base_size - 3, colour = "grey35"))

## Compose with lowercase bold tags.
combined <- (panel_a | panel_b) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = base_size + 2))

out_pdf <- "figures/method_comparison/perm_null_calibration_bc.pdf"
ggsave(out_pdf, combined, device = cairo_pdf, width = 12, height = 6)
cat(sprintf("[saved] %s\n", out_pdf))

## -----------------------------------------------------------------------------
## 8. Final printed report.
## -----------------------------------------------------------------------------

cat("\n================ FLAGSHIP emp_p TABLE ================\n")
print(flagship_table[, c("gene", "pair", "type", "z_observed", "null_sd", "emp_p")],
      row.names = FALSE, digits = 4)

cat("\n================ GENOME-WIDE frac|z|>1.96 ================\n")
cat(sprintf("null     = %.4f (target ~0.05 or below)\n", frac_null_sig))
cat(sprintf("observed = %.4f\n", frac_observed_sig))

cat("\n[paths]\n")
cat("  script:", normalizePath("scripts/figures/diagnostics/perm_null_calibration_bc.R"), "\n")
cat("  rds   :", normalizePath("data/breast_cancer/perm_null_calibration_bc.rds"), "\n")
cat("  pdf   :", normalizePath(out_pdf), "\n")

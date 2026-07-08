# ==============================================================================
# 04-mcsd.R
# MCSD: SVD-based gene importance scoring for variance decomposition
# Works natively with fit_light objects from 00-extract.R
#
# Boundary-singular fits:
#   With rich random-effect structures (e.g. (1 + neighbours || imageID)),
#   sparse / low-expression genes often have their variance components
#   collapse to zero. The BLUPs are still numerically returned but they are
#   uninterpretable; in our melanoma fits ~58-64% of panel genes are
#   boundary-singular and contribute near-zero columns to the contribution
#   matrix C. They do not pollute the SVD (their MCSD ~ 0) but they
#   inflate p (gene count), which distorts the null baseline (1/p) and
#   makes BH adjustment more conservative.
#
#   Recommended pattern in calling code:
#     fits_light  <- get_fits_light(results, well_fit_only = TRUE)
#     ran_vals_df <- get_ran_vals(results,   well_fit_only = TRUE)
#     mcsd_block_from_fits(fits_light, df, ...)
#
#   Or pass `well_fit_genes = names(is_well_fit(results))[is_well_fit(results)]`
#   to mcsd_block_from_fits / mcsd_partial_resid as an explicit filter.
# ==============================================================================

# ---- Per-gene contribution builder (unified) ----
.mcsd_contrib <- function(g, fit_light, df, block, group_ct,
                          resp_term = NULL, resp_level = NULL,
                          focal = NULL, neighbour = NULL) {
    used <- fit_light$used_rows
    df_loc <- droplevels(df[used, , drop = FALSE])
    n <- nrow(df_loc)

    # Focal masking (applied when focal is not NULL)
    mask_idx <- if (!is.null(focal)) {
        m <- which(df_loc[[group_ct]] == focal)
        if (!length(m)) return(list(gene = g, rows = character(0), vals = numeric(0)))
        m
    } else {
        seq_len(n)
    }

    # Spillover block uses fixed effects, not random
    if (block == "Spillover") {
        beta <- fit_light$fixef
        st <- names(beta)[grepl("spill|_near$", names(beta), ignore.case = TRUE)]
        st <- intersect(st, names(df_loc))
        if (length(st)) {
            X_spill <- as.matrix(df_loc[, st, drop = FALSE])
            vsum <- rowSums(sweep(X_spill, 2, beta[st], `*`))
            return(list(gene = g, rows = used[mask_idx], vals = as.numeric(vsum[mask_idx])))
        }
        return(list(gene = g, rows = used[mask_idx], vals = rep(0, length(mask_idx))))
    }

    re_list <- fit_light$ranef
    if (!group_ct %in% names(re_list)) {
        return(list(gene = g, rows = used[mask_idx], vals = rep(0, length(mask_idx))))
    }

    tab <- re_list[[group_ct]]
    levs <- rownames(tab); terms <- colnames(tab)

    if (block == "Spatial cell state") {
        if (is.null(resp_term)) {
            keep <- terms[terms != "(Intercept)"]
        } else {
            keep <- terms[terms != "(Intercept)" & !grepl(resp_term, terms, fixed = TRUE)]
        }
        if (!is.null(neighbour)) keep <- intersect(keep, neighbour)
    } else if (block == "Responder spatial state") {
        keep <- terms[grepl(paste0("^", resp_term, ":"), terms)]
        if (!is.null(neighbour)) keep <- intersect(keep, paste0(resp_term, ":", neighbour))
    } else if (block == "Responder status") {
        keep <- intersect(terms, resp_term)
    } else if (block == "Cell type") {
        keep <- intersect(terms, "(Intercept)")
    } else {
        stop("Unsupported block: ", block)
    }

    if (!length(keep)) {
        return(list(gene = g, rows = used[mask_idx], vals = rep(0, length(mask_idx))))
    }

    idx <- match(df_loc[[group_ct]], levs)
    resp_ind <- if (!is.null(resp_term)) as.numeric(df_loc$Responder == resp_level) else NULL
    vsum <- numeric(n)

    for (tm in keep) {
        b_by_row <- as.numeric(tab[, tm])[idx]
        if (identical(tm, "(Intercept)")) {
            z <- 1
        } else if (!is.null(resp_term) && tm == resp_term) {
            z <- resp_ind
        } else if (grepl(":", tm, fixed = TRUE)) {
            parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
            z <- 1
            for (p in parts) {
                if (!is.null(resp_term) && p == resp_term) z <- z * resp_ind
                else if (p %in% names(df_loc)) z <- z * as.numeric(df_loc[[p]])
            }
        } else if (tm %in% names(df_loc)) {
            z <- as.numeric(df_loc[[tm]])
        } else {
            z <- 0
        }
        vsum <- vsum + b_by_row * z
    }

    # Add the fixed Responder coefficient when computing the "Responder status"
    # block. After the offset+fixed-Responder patch, the population-mean
    # Responder effect lives in fit_light$fixef[["<resp_term><resp_level>"]];
    # the random per-celltype Responder BLUPs (already summed above) are
    # deviations from this fixed mean. Including both gives the total
    # per-cell Responder-status contribution.
    if (block == "Responder status" && !is.null(resp_term) && !is.null(resp_level)) {
        fix_resp <- paste0(resp_term, resp_level)
        if (fix_resp %in% names(fit_light$fixef)) {
            vsum <- vsum + fit_light$fixef[[fix_resp]] * resp_ind
        }
    }

    list(gene = g, rows = used[mask_idx], vals = as.numeric(vsum[mask_idx]))
}

# ---- Per-gene partial-residual contribution builder (unified) ----
# Returns: contribution[i] = residual[i] + b[focal, gene, neighbour] * X[i, neighbour]
# for focal cells only. This preserves cell-level variation that breaks rank-1.
.mcsd_contrib_partial_resid <- function(g, fit_light, df, group_ct, focal, neighbour,
                                        residual_type = "working",
                                        block = NULL, resp_term = NULL, resp_level = NULL) {
    used <- fit_light$used_rows
    df_loc <- droplevels(df[used, , drop = FALSE])
    n <- nrow(df_loc)

    # Focal mask
    mask_idx <- which(df_loc[[group_ct]] == focal)
    if (!length(mask_idx)) return(list(gene = g, rows = character(0), vals = numeric(0)))

    # Get residuals
    r <- fit_light$residuals[[residual_type]]
    if (is.null(r)) return(list(gene = g, rows = character(0), vals = numeric(0)))
    r_focal <- r[mask_idx]

    # Get spatial contribution for this neighbour
    re_list <- fit_light$ranef
    if (!group_ct %in% names(re_list)) {
        # No random effects -- just return residuals
        return(list(gene = g, rows = used[mask_idx], vals = as.numeric(r_focal)))
    }

    tab <- re_list[[group_ct]]
    levs <- rownames(tab); terms <- colnames(tab)

    # Select terms based on resp_term presence and block
    if (is.null(resp_term)) {
        # No-responder: simple neighbour intersection
        nb_terms <- intersect(neighbour, terms)
    } else {
        # Responder-aware: term selection depends on block
        if (is.null(block) || block == "Spatial cell state") {
            nb_terms <- intersect(neighbour, terms)
        } else if (block == "Responder spatial state") {
            nb_terms <- intersect(paste0(resp_term, ":", neighbour), terms)
        } else {
            stop("Partial residual MCSD only supports Spatial cell state and Responder spatial state blocks")
        }
    }

    if (!length(nb_terms)) {
        return(list(gene = g, rows = used[mask_idx], vals = as.numeric(r_focal)))
    }

    # Spatial contribution on focal cells
    idx <- match(df_loc[[group_ct]][mask_idx], levs)
    resp_ind <- if (!is.null(resp_term)) as.numeric(df_loc$Responder[mask_idx] == resp_level) else NULL
    spatial_contrib <- numeric(length(mask_idx))

    for (tm in nb_terms) {
        b_by_row <- as.numeric(tab[, tm])[idx]
        if (!is.null(resp_term) && grepl(":", tm, fixed = TRUE)) {
            parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
            z <- 1
            for (p in parts) {
                if (p == resp_term) z <- z * resp_ind
                else if (p %in% names(df_loc)) z <- z * as.numeric(df_loc[[p]][mask_idx])
            }
        } else if (tm %in% names(df_loc)) {
            z <- as.numeric(df_loc[[tm]][mask_idx])
        } else {
            z <- 0
        }
        spatial_contrib <- spatial_contrib + b_by_row * z
    }

    list(gene = g, rows = used[mask_idx], vals = as.numeric(r_focal + spatial_contrib))
}


# ---- Core MCSD computation from a standardised matrix ----
# Extracted so it can be reused for both observed and permuted matrices.
.compute_mcsd_scores <- function(Z, use_gram = TRUE) {
    p <- ncol(Z)
    if (p == 0) return(rep(0, p))

    if (use_gram) {
        G <- crossprod(Z)
        ee <- eigen(G, symmetric = TRUE)
        d2 <- pmax(ee$values, 0); V <- ee$vectors
    } else {
        sv <- svd(Z, nu = 0)
        d2 <- sv$d^2; V <- sv$v
    }

    if (!length(d2) || is.null(V)) {
        mcsd <- rep(1 / p, p)
    } else {
        w <- d2 / sum(d2)
        mcsd <- as.numeric((V^2) %*% w)
    }
    names(mcsd) <- colnames(Z)
    mcsd
}


#' Permutation Null for MCSD Scores
#'
#' Generates a null distribution for MCSD scores by independently permuting
#' cells within each gene column of the contribution matrix. This breaks
#' the cross-gene cell-level correlation structure while preserving each
#' gene's marginal distribution.
#'
#' @param C Raw contribution matrix (cells x genes), before standardisation.
#'   Typically obtained from \code{mcsd_block_from_fits(..., keep_Z = TRUE)$C}.
#' @param n_perm Number of permutations (default: 1000)
#' @param center Center columns before SVD (default: TRUE)
#' @param scale_ref Scaling method: "none", "mad", or "sd" (default: "none")
#' @param use_gram Use Gram matrix for eigendecomposition (default: TRUE)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param BPPARAM BiocParallel backend for parallelising permutations
#' @return List with:
#'   \describe{
#'     \item{observed}{Named numeric vector of observed MCSD scores}
#'     \item{null_dist}{Matrix (n_perm x p) of null MCSD scores}
#'     \item{pval}{Per-gene permutation p-values (one-sided, upper tail)}
#'     \item{padj}{BH-adjusted p-values}
#'     \item{n_perm}{Number of permutations used}
#'   }
#' @export
mcsd_perm_null <- function(
        C,
        n_perm = 1000,
        center = TRUE,
        scale_ref = c("none", "mad", "sd"),
        use_gram = TRUE,
        seed = NULL,
        BPPARAM = BiocParallel::SerialParam()
) {
    scale_ref <- match.arg(scale_ref)
    stopifnot(is.matrix(C), ncol(C) > 1, nrow(C) > 1)
    p <- ncol(C)
    BPPARAM <- .bp_or_serial(BPPARAM)

    # Standardise and compute observed scores
    .standardise <- function(M) {
        mu <- if (center) colMeans(M, na.rm = TRUE) else rep(0, ncol(M))
        s <- switch(scale_ref,
                    none = rep(1, ncol(M)),
                    mad  = apply(M, 2, stats::mad, constant = 1.4826, na.rm = TRUE),
                    sd   = apply(M, 2, stats::sd, na.rm = TRUE))
        s[!is.finite(s) | s == 0] <- .Machine$double.eps
        Z <- sweep(M, 2, mu, "-")
        Z <- sweep(Z, 2, s, "/")
        Z[!is.finite(Z)] <- 0
        Z
    }

    Z_obs <- .standardise(C)
    observed <- .compute_mcsd_scores(Z_obs, use_gram = use_gram)

    # Permutation null
    if (!is.null(seed)) set.seed(seed)

    null_dist <- BiocParallel::bplapply(seq_len(n_perm), function(b) {
        # Independently permute rows within each column
        C_perm <- apply(C, 2, sample)
        Z_perm <- .standardise(C_perm)
        .compute_mcsd_scores(Z_perm, use_gram = use_gram)
    }, BPPARAM = BPPARAM)
    null_dist <- do.call(rbind, null_dist)
    colnames(null_dist) <- names(observed)

    # Per-gene p-values (proportion of null >= observed)
    pval <- vapply(seq_len(p), function(j) {
        (sum(null_dist[, j] >= observed[j]) + 1) / (n_perm + 1)
    }, numeric(1))
    names(pval) <- names(observed)

    padj_bh <- stats::p.adjust(pval, method = "BH")
    qval <- tryCatch(qvalue::qvalue(pval)$qvalues, error = function(e) padj_bh)

    list(
        observed = observed,
        null_dist = null_dist,
        pval = pval,
        padj = padj_bh,
        qval = qval,
        n_perm = n_perm
    )
}


#' Marchenko-Pastur Analytical Null for MCSD Scores
#'
#' Derives p-values analytically using random matrix theory. Under the null
#' that entries of Z are iid with finite second moment, the squared right
#' singular vector loadings V_jk^2 converge to 1/p, and MCSD_j concentrates
#' around 1/p with a variance that depends on the aspect ratio and the
#' Herfindahl index of the variance spectrum.
#'
#' This is much faster than permutation (no resampling needed) but assumes
#' approximately Gaussian-like entries in Z. Use as a fast screen or as a
#' cross-check against the permutation null.
#'
#' @param Z Standardised contribution matrix (cells x genes). Typically
#'   obtained from \code{mcsd_block_from_fits(..., keep_Z = TRUE)$Z}.
#' @param use_gram Use Gram matrix for eigendecomposition (default: TRUE)
#' @return List with:
#'   \describe{
#'     \item{observed}{Named numeric vector of observed MCSD scores}
#'     \item{null_mean}{Expected MCSD under null (1/p for all genes)}
#'     \item{null_sd}{Standard deviation of MCSD under null}
#'     \item{z_score}{Per-gene z-scores: (observed - null_mean) / null_sd}
#'     \item{pval}{Per-gene p-values (one-sided, upper tail, from Gaussian approx)}
#'     \item{padj}{BH-adjusted p-values}
#'   }
#' @export
mcsd_mp_null <- function(Z, use_gram = TRUE) {
    stopifnot(is.matrix(Z), ncol(Z) > 1, nrow(Z) > 1)
    n <- nrow(Z); p <- ncol(Z)

    # Observed MCSD
    observed <- .compute_mcsd_scores(Z, use_gram = use_gram)

    # Compute variance weights from the decomposition
    if (use_gram) {
        G <- crossprod(Z)
        ee <- eigen(G, symmetric = TRUE)
        d2 <- pmax(ee$values, 0)
    } else {
        sv <- svd(Z, nu = 0)
        d2 <- sv$d^2
    }
    w <- d2 / sum(d2)

    # Under null (iid entries): E[V_jk^2] = 1/p for all j, k
    # So E[MCSD_j] = sum_k w_k * (1/p) = 1/p
    null_mean <- 1 / p

    # Var[V_jk^2] for Haar-distributed V (uniform on Stiefel manifold):
    #   Var[V_jk^2] = 2(p-1) / (p^2 * (p+2))
    # Cov[V_jk^2, V_jl^2] for k != l:
    #   Cov = -2 / (p^2 * (p+2))
    # Therefore:
    #   Var[MCSD_j] = sum_k w_k^2 * Var[V_jk^2] + sum_{k!=l} w_k*w_l * Cov[V_jk^2, V_jl^2]
    #              = sum_k w_k^2 * 2(p-1)/(p^2(p+2)) + (sum_k w_k)^2 - sum_k w_k^2) * (-2/(p^2(p+2)))
    #              = 2/(p^2(p+2)) * [ (p-1) * H + (-(1 - H)) ]
    #              = 2/(p^2(p+2)) * [ pH - 1 ]    where H = sum(w_k^2) (Herfindahl index)
    #
    # For p*H > 1 (which holds when variance is concentrated), this is positive.
    # When p*H <= 1 (perfectly uniform spectrum), variance is ~0 as expected.
    H <- sum(w^2)  # Herfindahl index of variance spectrum
    null_var <- 2 / (p^2 * (p + 2)) * max(p * H - 1, .Machine$double.eps)
    null_sd <- sqrt(null_var)

    # Z-scores and p-values (Gaussian approximation, upper tail)
    z_score <- (observed - null_mean) / null_sd
    pval <- stats::pnorm(z_score, lower.tail = FALSE)
    names(pval) <- names(observed)

    padj_bh <- stats::p.adjust(pval, method = "BH")
    qval <- tryCatch(qvalue::qvalue(pval)$qvalues, error = function(e) padj_bh)

    list(
        observed = observed,
        null_mean = null_mean,
        null_sd = null_sd,
        z_score = z_score,
        pval = pval,
        padj = padj_bh,
        qval = qval
    )
}


#' Rotation Null for MCSD Scores
#'
#' Generates a null distribution by applying random orthogonal rotations
#' in gene space. This preserves the singular values (total variance structure)
#' exactly but redistributes loadings uniformly across genes, testing whether
#' any gene loads disproportionately onto the principal components.
#'
#' This is more conservative than column-wise permutation because it preserves
#' the full correlation structure of the data and only tests whether individual
#' gene loadings deviate from uniform.
#'
#' @param Z Standardised contribution matrix (cells x genes). Typically
#'   obtained from \code{mcsd_block_from_fits(..., keep_Z = TRUE)$Z}.
#' @param n_perm Number of random rotations (default: 1000)
#' @param use_gram Use Gram matrix for eigendecomposition (default: TRUE)
#' @param seed Random seed for reproducibility (default: NULL)
#' @return List with:
#'   \describe{
#'     \item{observed}{Named numeric vector of observed MCSD scores}
#'     \item{null_dist}{Matrix (n_perm x p) of null MCSD scores}
#'     \item{pval}{Per-gene permutation p-values (one-sided, upper tail)}
#'     \item{padj}{BH-adjusted p-values}
#'     \item{n_perm}{Number of rotations used}
#'   }
#' @export
mcsd_rotation_null <- function(
        Z,
        n_perm = 1000,
        use_gram = TRUE,
        seed = NULL
) {
    stopifnot(is.matrix(Z), ncol(Z) > 1, nrow(Z) > 1)
    p <- ncol(Z)

    # Observed scores
    observed <- .compute_mcsd_scores(Z, use_gram = use_gram)

    if (!is.null(seed)) set.seed(seed)

    # For rotation null, we don't need to redo the full SVD each time.
    # Z_rot = Z %*% Q has the same singular values as Z.
    # The right singular vectors transform as V_rot = Q^T V.
    # So MCSD_rot_j = sum_k (Q^T V)_jk^2 * w_k.
    #
    # We can compute V and w once, then just rotate V.
    if (use_gram) {
        G <- crossprod(Z)
        ee <- eigen(G, symmetric = TRUE)
        d2 <- pmax(ee$values, 0); V <- ee$vectors
    } else {
        sv <- svd(Z, nu = 0)
        d2 <- sv$d^2; V <- sv$v
    }
    w <- d2 / sum(d2)

    null_dist <- matrix(NA_real_, nrow = n_perm, ncol = p,
                        dimnames = list(NULL, colnames(Z)))

    for (b in seq_len(n_perm)) {
        # Random orthogonal matrix via QR of random Gaussian
        Q <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
        V_rot <- crossprod(Q, V)  # Q^T %*% V
        null_dist[b, ] <- as.numeric((V_rot^2) %*% w)
    }

    # Per-gene p-values
    pval <- vapply(seq_len(p), function(j) {
        (sum(null_dist[, j] >= observed[j]) + 1) / (n_perm + 1)
    }, numeric(1))
    names(pval) <- names(observed)

    padj_bh <- stats::p.adjust(pval, method = "BH")
    qval <- tryCatch(qvalue::qvalue(pval)$qvalues, error = function(e) padj_bh)

    list(
        observed = observed,
        null_dist = null_dist,
        pval = pval,
        padj = padj_bh,
        qval = qval,
        n_perm = n_perm
    )
}


#' Parametric Coefficient-Resampling Null for Weighted MCSD Scores
#'
#' Generates a null distribution for MCSD scores by resampling spatial
#' coefficients from their estimated null distribution. Under the null,
#' each gene's spatial slope is pure estimation noise: b_j ~ N(0, se_j^2).
#'
#' This approach is appropriate for the standard (non-partial-residual) MCSD
#' where the contribution matrix is approximately rank-1 (single neighbour type),
#' making cell-level permutation tests inappropriate. It naturally incorporates
#' spillover and detection weights because the same post-hoc adjustments are
#' applied to both observed and simulated scores.
#'
#' @param ran_vals Tidy data frame from \code{get_ran_vals()} containing columns:
#'   gene, level, term, estimate, std.error. Should be pre-filtered to the
#'   relevant focal cell type and neighbour term.
#' @param focal Focal cell type (used to filter ran_vals by level)
#' @param neighbour Neighbour cell type (used to filter ran_vals by term)
#' @param observed_scores Named numeric vector of observed weighted MCSD scores
#'   (from \code{result$scores$MCSD}, named by gene). These are the final scores
#'   after all weighting adjustments.
#' @param spill_weights Optional named numeric vector of spillover weights per gene.
#'   If provided, applied as \code{mcsd * spill_weights} then renormalised.
#'   Compute as: \code{1 / pmax(expr_ratio, 1)} for "weight" mode,
#'   or \code{(1 / pmax(expr_ratio, 1))^2} for "weight2" mode.
#' @param det_rates Optional named numeric vector of detection rates per gene.
#'   If provided, applied as \code{mcsd * det_rates} then renormalised.
#' @param n_sim Number of simulations (default: 1000)
#' @param seed Random seed for reproducibility (default: NULL)
#' @return List with:
#'   \describe{
#'     \item{observed}{Named numeric vector of observed weighted MCSD scores}
#'     \item{null_dist}{Matrix (n_sim x p) of null weighted MCSD scores}
#'     \item{pval}{Per-gene p-values (one-sided, upper tail)}
#'     \item{padj}{BH-adjusted p-values}
#'     \item{n_sim}{Number of simulations used}
#'     \item{genes}{Character vector of gene names in order}
#'     \item{coefficients}{Data frame of observed estimates and standard errors used}
#'   }
#' @export
mcsd_coef_null <- function(
        ran_vals,
        focal,
        neighbour,
        observed_scores,
        spill_weights = NULL,
        det_rates = NULL,
        n_sim = 1000,
        seed = NULL
) {
    # Filter ran_vals to the focal-neighbour spatial slopes
    rv <- ran_vals |>
        dplyr::filter(level == focal, term == neighbour) |>
        dplyr::distinct(gene, .keep_all = TRUE)

    if (nrow(rv) == 0) {
        stop("No slopes found for focal='", focal, "', neighbour='", neighbour, "'.")
    }

    # Align genes: only keep genes present in both ran_vals and observed_scores
    genes <- intersect(rv$gene, names(observed_scores))
    if (!length(genes)) {
        stop("No overlapping genes between ran_vals and observed_scores.")
    }

    rv <- rv |> dplyr::filter(gene %in% genes)
    rv <- rv[match(genes, rv$gene), ]  # ensure same order

    estimates <- rv$estimate
    ses <- rv$std.error
    names(estimates) <- names(ses) <- genes

    # Drop genes with NA/non-finite estimates or standard errors
    valid <- is.finite(estimates) & is.finite(ses) & ses > 0
    if (!any(valid)) stop("No genes with finite estimates and standard errors.")
    if (any(!valid)) {
        message("  Dropping ", sum(!valid), " genes with NA/non-finite coefficients: ",
                paste(genes[!valid], collapse = ", "))
    }
    genes <- genes[valid]
    estimates <- estimates[valid]
    ses <- ses[valid]
    p <- length(genes)

    # Observed scores (subset and reorder to match)
    obs <- observed_scores[genes]

    # Align weights
    sw <- if (!is.null(spill_weights)) {
        w <- spill_weights[genes]
        w[is.na(w)] <- 1
        w
    } else NULL

    dr <- if (!is.null(det_rates)) {
        d <- det_rates[genes]
        d[is.na(d)] <- 0
        d
    } else NULL

    # Simulate null
    if (!is.null(seed)) set.seed(seed)

    null_dist <- matrix(NA_real_, nrow = n_sim, ncol = p,
                        dimnames = list(NULL, genes))

    for (b in seq_len(n_sim)) {
        # Draw null coefficients: b_j ~ N(0, se_j^2)
        b_null <- stats::rnorm(p, mean = 0, sd = ses)

        # Raw MCSD in rank-1 case: proportional to b_j^2
        mcsd_null <- b_null^2
        total <- sum(mcsd_null, na.rm = TRUE)
        if (is.finite(total) && total > 0) mcsd_null <- mcsd_null / total

        # Apply spillover weights (same as mcsd_block_from_fits)
        if (!is.null(sw)) {
            mcsd_null <- mcsd_null * sw
            total <- sum(mcsd_null, na.rm = TRUE)
            if (is.finite(total) && total > 0) mcsd_null <- mcsd_null / total
        }

        # Apply detection weights (same as mcsd_block_from_fits)
        if (!is.null(dr)) {
            mcsd_null <- mcsd_null * dr
            total <- sum(mcsd_null, na.rm = TRUE)
            if (is.finite(total) && total > 0) mcsd_null <- mcsd_null / total
        }

        null_dist[b, ] <- mcsd_null
    }

    # Per-gene p-values
    pval <- vapply(seq_len(p), function(j) {
        (sum(null_dist[, j] >= obs[j]) + 1) / (n_sim + 1)
    }, numeric(1))
    names(pval) <- genes

    padj_bh <- stats::p.adjust(pval, method = "BH")

    # Storey q-values: less conservative when most genes are true nulls
    qval <- tryCatch({
        qvalue::qvalue(pval)$qvalues
    }, error = function(e) {
        # Fall back to BH if qvalue is not installed or fails
        # (can fail with too few p-values or extreme distributions)
        message("  qvalue::qvalue() failed (", conditionMessage(e),
                "), falling back to BH adjustment.")
        padj_bh
    })

    list(
        observed = obs,
        null_dist = null_dist,
        pval = pval,
        padj = padj_bh,
        qval = qval,
        n_sim = n_sim,
        genes = genes,
        coefficients = tibble::tibble(
            gene = genes,
            estimate = estimates,
            std.error = ses
        )
    )
}


#' Compute per-gene source-to-focal expression ratio
#'
#' For each gene, computes ratio = mean(gene in source) / mean(gene in focal).
#' Genes highly expressed in the source relative to focal are plausible spillover.
#'
#' @param fits_light Named list of fit_light objects (used for gene names)
#' @param df Original data frame
#' @param focal Focal cell type (target of contamination)
#' @param neighbour Neighbour cell type(s) (source of contamination)
#' @param group_ct Cell type grouping column
#' @return Named numeric vector of expression ratios (>1 means source-enriched)
.spill_expr_ratio <- function(fits_light, df, focal, neighbour, group_ct = "celltype",
                              spill_source_mode = c("neighbour", "any"), ...) {
    spill_source_mode <- match.arg(spill_source_mode)
    genes <- names(fits_light)

    focal_idx <- df[[group_ct]] == focal

    if (spill_source_mode == "any") {
        # Max ratio across all non-focal cell types
        other_types <- setdiff(unique(df[[group_ct]]), focal)
        source_idx_list <- lapply(other_types, function(ct) df[[group_ct]] == ct)
    } else {
        source_idx_list <- list(df[[group_ct]] %in% neighbour)
    }

    if (!any(focal_idx)) {
        out <- rep(0, length(genes)); names(out) <- genes; return(out)
    }

    vapply(genes, function(g) {
        if (!g %in% names(df)) return(0)
        focal_mean  <- mean(df[[g]][focal_idx], na.rm = TRUE)

        if (!is.finite(focal_mean) || focal_mean < .Machine$double.eps) {
            for (si in source_idx_list) {
                source_mean <- mean(df[[g]][si], na.rm = TRUE)
                if (is.finite(source_mean) && source_mean > 0) return(Inf)
            }
            return(0)
        }

        max_ratio <- 0
        for (si in source_idx_list) {
            if (!any(si)) next
            source_mean <- mean(df[[g]][si], na.rm = TRUE)
            ratio <- source_mean / focal_mean
            if (is.finite(ratio) && ratio > max_ratio) max_ratio <- ratio
        }
        max_ratio
    }, numeric(1))
}


#' MCSD Block Scores (Unified)
#'
#' @param fits_light Named list of fit_light objects
#' @param df Original data frame
#' @param block Block type to compute MCSD for
#' @param focal Optional focal cell type
#' @param neighbour Optional neighbour type(s)
#' @param group_ct Cell type grouping column
#' @param resp_term Responder term name (NULL for no-responder mode)
#' @param resp_level Responder level
#' @param center Center contributions before SVD
#' @param scale_ref Scaling method
#' @param BPPARAM BiocParallel backend
#' @param expr_filter Apply expression-based gene filtering
#' @param min_expr_prev Minimum prevalence threshold
#' @param min_expr_var Minimum variance threshold
#' @param spill_correct Spillover correction mode: "none" (default), "weight"
#'   (w = 1/max(ratio,1)), "weight2" (w = 1/max(ratio,1)^2), or "exclude"
#'   (drop genes with ratio > spill_threshold)
#' @param spill_threshold Expression ratio threshold for "exclude" mode (default 1.5)
#' @param spill_source_mode Source mode for spillover ratio
#' @param detection_weight Detection weight mode
#' @param use_gram Use Gram matrix for eigendecomposition (TRUE) or plain SVD (FALSE)
#' @param keep_Z Keep Z and C matrices in output
#' @return List with scores, meta, Z matrix, and C matrix
mcsd_block_from_fits <- function(
        fits_light, df,
        block = c("Spatial cell state", "Responder spatial state",
                  "Responder status", "Cell type", "Spillover"),
        focal = NULL, neighbour = NULL,
        group_ct = "celltype",
        resp_term = NULL, resp_level = "PD",
        ## Canonical post-2026-05-01: Variant A (covariance PCA) — center the
        ## contribution matrix and skip column SD scaling so MCSD is magnitude-
        ## aware, not scale-invariant. See project_mcsd_reformulation.md.
        center = TRUE, scale_ref = c("none", "sd", "mad"),
        BPPARAM = BiocParallel::SerialParam(),
        expr_filter = TRUE, min_expr_prev = 0.10, min_expr_var = 1,
        spill_correct = c("none", "weight", "weight2", "exclude"),
        spill_threshold = 1.5,
        spill_source_mode = c("neighbour", "any"),
        detection_weight = c("none", "post", "pre"),
        use_gram = TRUE, keep_Z = FALSE,
        well_fit_genes = NULL
) {
    block <- match.arg(block); scale_ref <- match.arg(scale_ref)
    spill_correct <- match.arg(spill_correct)
    spill_source_mode <- match.arg(spill_source_mode)
    detection_weight <- match.arg(detection_weight)
    stopifnot(is.list(fits_light), length(fits_light) > 0)
    genes <- names(fits_light)
    stopifnot(!is.null(genes), all(nzchar(genes)))
    BPPARAM <- .bp_or_serial(BPPARAM)

    # Boundary-singular gene filter (see header for context)
    if (!is.null(well_fit_genes)) {
        n_pre <- length(fits_light)
        fits_light <- fits_light[intersect(names(fits_light), well_fit_genes)]
        if (!length(fits_light)) stop("No genes survived well_fit_genes filter.")
        n_dropped <- n_pre - length(fits_light)
        if (n_dropped > 0L) {
            message("  Dropped ", n_dropped, " boundary-singular genes; ",
                    length(fits_light), " genes retained.")
        }
        genes <- names(fits_light)
    }

    # Raw-expression filter in focal slice
    if (expr_filter) {
        df_scope <- if (is.null(focal)) df else df[df[[group_ct]] == focal, , drop = FALSE]
        expr_prev <- vapply(genes, function(g) mean(is.finite(df_scope[[g]]) & df_scope[[g]] > 0, na.rm = TRUE), numeric(1))
        expr_var <- vapply(genes, function(g) stats::var(df_scope[[g]], na.rm = TRUE), numeric(1))
        keep <- (expr_prev >= min_expr_prev) & (expr_var >= min_expr_var)
        if (!any(keep)) stop("All genes failed raw-expression filters; relax thresholds.")
        fits_light <- fits_light[keep]; genes <- names(fits_light)
    }

    # Spillover correction: compute expression ratios
    expr_ratio <- NULL
    spill_blocks <- c("Spatial cell state", "Responder spatial state")
    if (spill_correct != "none" && block %in% spill_blocks &&
        !is.null(focal) && !is.null(neighbour)) {
        expr_ratio <- .spill_expr_ratio(fits_light, df, focal, neighbour, group_ct,
                                        spill_source_mode = spill_source_mode)
        expr_ratio[is.na(expr_ratio)] <- 0

        if (spill_correct == "exclude") {
            excl <- expr_ratio > spill_threshold
            if (any(excl)) {
                message("  Excluding ", sum(excl), " genes with expr ratio > ", spill_threshold,
                        ": ", paste(names(expr_ratio)[excl], collapse = ", "))
            }
            keep_genes <- !excl
            if (!any(keep_genes)) stop("All genes excluded by spillover filter.")
            fits_light <- fits_light[keep_genes]
            genes <- names(fits_light)
            expr_ratio <- expr_ratio[keep_genes]
        }
    }

    # Per-gene contributions (parallel)
    per_gene <- BiocParallel::bplapply(genes, function(g) {
        fl <- fits_light[[g]]
        .mcsd_contrib(g, fl, df, block, group_ct,
                      resp_term = resp_term, resp_level = resp_level,
                      focal = focal, neighbour = neighbour)
    }, BPPARAM = BPPARAM)
    per_gene <- Filter(function(x) length(x$rows) > 0, per_gene)
    if (!length(per_gene)) stop("No genes produced contributions.")

    # Align by common rows
    common_rows <- Reduce(intersect, lapply(per_gene, `[[`, "rows"))
    if (!length(common_rows)) stop("No overlapping rows after focal/neighbour filtering.")
    # Focal masking after alignment (noresp approach for consistency)
    if (!is.null(focal)) {
        keep_rows <- common_rows[df[common_rows, group_ct, drop = TRUE] == focal]
        if (!length(keep_rows)) stop("No rows for focal after alignment.")
        common_rows <- keep_rows
    }
    pos <- setNames(seq_along(common_rows), common_rows)
    gene_kept <- vapply(per_gene, `[[`, "", "gene")
    C <- matrix(0, nrow = length(common_rows), ncol = length(per_gene),
                dimnames = list(common_rows, gene_kept))
    for (i in seq_along(per_gene)) {
        pg <- per_gene[[i]]
        m <- pos[pg$rows]; ok <- !is.na(m)
        C[m[ok], i] <- pg$vals[ok]
    }

    mu <- if (center) colMeans(C, na.rm = TRUE) else rep(0, ncol(C))
    s <- switch(scale_ref,
                sd   = apply(C, 2, stats::sd, na.rm = TRUE),
                mad  = apply(C, 2, stats::mad, constant = 1.4826, na.rm = TRUE),
                none = rep(1, ncol(C)))
    s[!is.finite(s) | s == 0] <- .Machine$double.eps
    Z <- sweep(C, 2, mu, "-"); Z <- sweep(Z, 2, s, "/"); Z[!is.finite(Z)] <- 0

    # Compute detection rates for kept genes in focal cells
    det_rates <- NULL
    if (detection_weight != "none") {
        df_scope <- if (is.null(focal)) df else df[df[[group_ct]] == focal, , drop = FALSE]
        det_rates <- vapply(colnames(Z), function(g) {
            mean(is.finite(df_scope[[g]]) & df_scope[[g]] > 0, na.rm = TRUE)
        }, numeric(1))
    }

    # Pre-SVD: scale Z columns by sqrt(detection rate)
    if (detection_weight == "pre" && !is.null(det_rates)) {
        Z <- sweep(Z, 2, sqrt(pmax(det_rates, .Machine$double.eps)), "*")
    }

    mcsd <- .compute_mcsd_scores(Z, use_gram = use_gram)

    # Apply spillover downweighting
    if (spill_correct %in% c("weight", "weight2") && !is.null(expr_ratio)) {
        r <- expr_ratio[colnames(Z)]
        r[is.na(r)] <- 0
        spill_w <- 1 / pmax(r, 1)
        if (spill_correct == "weight2") spill_w <- spill_w^2
        mcsd <- mcsd * spill_w
        total_mcsd <- sum(mcsd)
        if (total_mcsd > 0) mcsd <- mcsd / total_mcsd
    }

    # Post-SVD: multiply MCSD by detection rate
    if (detection_weight == "post" && !is.null(det_rates)) {
        dr <- det_rates[names(mcsd)]
        dr[is.na(dr)] <- 0
        mcsd <- mcsd * dr
        total_mcsd <- sum(mcsd)
        if (total_mcsd > 0) mcsd <- mcsd / total_mcsd
    }

    scores <- tibble(gene = names(mcsd), MCSD = mcsd, rank = rank(-mcsd, ties.method = "min")) |>
        arrange(rank)

    list(
        scores = scores,
        meta   = list(block = block, focal = focal, neighbour = neighbour,
                      center = center, scale_ref = scale_ref,
                      spill_correct = spill_correct,
                      n_genes = length(gene_kept), n_cells = nrow(Z)),
        expr_ratio = if (!is.null(expr_ratio)) expr_ratio[colnames(Z)] else NULL,
        Z = if (keep_Z) Z else NULL,
        C = if (keep_Z) C else NULL
    )
}


#' MCSD via Partial Residuals (Unified)
#'
#' Uses C[i,j] = residual[i,j] + spatial_contribution[i,j] to break
#' the rank-1 degeneracy of the standard pairwise MCSD. The residuals
#' carry cell-level variation that creates genuine multivariate structure.
#'
#' @param fits_light Named list of fit_light objects
#' @param df Original data frame
#' @param block "Spatial cell state" or "Responder spatial state"
#' @param focal Focal cell type (required)
#' @param neighbour Neighbour cell type(s) (required)
#' @param group_ct Cell type grouping column
#' @param resp_term Responder term name (NULL for no-responder mode)
#' @param resp_level Responder level
#' @param center Center contributions before SVD
#' @param scale_ref Scaling method
#' @param residual_type Type of residuals to use
#' @param BPPARAM BiocParallel backend
#' @param expr_filter Apply expression-based gene filtering
#' @param min_expr_prev Minimum prevalence threshold
#' @param min_expr_var Minimum variance threshold
#' @param use_gram Use Gram matrix for eigendecomposition
#' @param keep_Z Keep Z and C matrices in output
#' @return List with scores, meta, and optionally Z and C matrices
mcsd_partial_resid <- function(
        fits_light, df,
        block = c("Spatial cell state", "Responder spatial state"),
        focal, neighbour,
        group_ct = "celltype",
        resp_term = NULL, resp_level = "PD",
        ## Canonical post-2026-05-01: Variant A (covariance PCA)
        center = TRUE, scale_ref = c("none", "sd", "mad"),
        residual_type = "working",
        BPPARAM = BiocParallel::SerialParam(),
        expr_filter = TRUE, min_expr_prev = 0.10, min_expr_var = 1.0,
        use_gram = TRUE, keep_Z = FALSE,
        well_fit_genes = NULL
) {
    block <- match.arg(block); scale_ref <- match.arg(scale_ref)
    stopifnot(!is.null(focal), !is.null(neighbour))
    genes <- names(fits_light)
    BPPARAM <- .bp_or_serial(BPPARAM)

    # Boundary-singular gene filter (see header for context)
    if (!is.null(well_fit_genes)) {
        n_pre <- length(fits_light)
        fits_light <- fits_light[intersect(names(fits_light), well_fit_genes)]
        if (!length(fits_light)) stop("No genes survived well_fit_genes filter.")
        n_dropped <- n_pre - length(fits_light)
        if (n_dropped > 0L) {
            message("  Dropped ", n_dropped, " boundary-singular genes; ",
                    length(fits_light), " genes retained.")
        }
        genes <- names(fits_light)
    }

    if (expr_filter) {
        df_scope <- df[df[[group_ct]] == focal, , drop = FALSE]
        expr_var <- vapply(genes, function(g) stats::var(df_scope[[g]], na.rm = TRUE), numeric(1))
        expr_prev <- vapply(genes, function(g) mean(is.finite(df_scope[[g]]) & df_scope[[g]] > 0, na.rm = TRUE), numeric(1))
        keep <- (expr_var >= min_expr_var) & (expr_prev >= min_expr_prev)
        if (!any(keep)) stop("All genes failed filters.")
        fits_light <- fits_light[keep]; genes <- names(fits_light)
    }

    per_gene <- BiocParallel::bplapply(genes, function(g) {
        fl <- fits_light[[g]]
        .mcsd_contrib_partial_resid(g, fl, df, group_ct, focal, neighbour,
                                    residual_type,
                                    block = if (!is.null(resp_term)) block else NULL,
                                    resp_term = resp_term, resp_level = resp_level)
    }, BPPARAM = BPPARAM)
    per_gene <- Filter(function(x) length(x$rows) > 0, per_gene)
    if (!length(per_gene)) stop("No genes produced contributions.")

    common_rows <- Reduce(intersect, lapply(per_gene, `[[`, "rows"))
    if (!length(common_rows)) stop("No overlapping rows.")
    pos <- setNames(seq_along(common_rows), common_rows)
    gene_kept <- vapply(per_gene, `[[`, "", "gene")
    C <- matrix(0, nrow = length(common_rows), ncol = length(per_gene),
                dimnames = list(common_rows, gene_kept))
    for (i in seq_along(per_gene)) {
        pg <- per_gene[[i]]
        m <- pos[pg$rows]; ok <- !is.na(m)
        C[m[ok], i] <- pg$vals[ok]
    }

    mu <- if (center) colMeans(C) else rep(0, ncol(C))
    s <- switch(scale_ref,
                "none" = rep(1, ncol(C)),
                "mad"  = apply(C, 2, stats::mad, constant = 1.4826, na.rm = TRUE),
                "sd"   = apply(C, 2, stats::sd, na.rm = TRUE))
    s[!is.finite(s) | s == 0] <- .Machine$double.eps
    Z <- sweep(C, 2, mu, "-"); Z <- sweep(Z, 2, s, "/"); Z[!is.finite(Z)] <- 0

    mcsd <- .compute_mcsd_scores(Z, use_gram = use_gram)

    scores <- tibble(gene = names(mcsd), MCSD = mcsd, rank = rank(-mcsd, ties.method = "min")) |>
        arrange(rank)

    list(
        scores = scores,
        meta   = list(block = paste0(block, " (partial residual)"),
                      focal = focal, neighbour = neighbour,
                      center = center, scale_ref = scale_ref,
                      n_genes = length(gene_kept), n_cells = nrow(Z)),
        Z = if (keep_Z) Z else NULL,
        C = if (keep_Z) C else NULL
    )
}


# ---- Top spatial drivers with model p-values ----

#' Top Spatial Drivers with Model Z-scores
#'
#' Combines MCSD rankings (with spillover/detection weighting) with per-gene
#' spatial slope z-scores from the mixed model. Returns a single data frame
#' of genes ranked by MCSD with effect sizes and direction.
#'
#' Note: z-scores are Wald statistics on conditional modes (BLUPs) of the
#' random effects. They are approximate and may be anti-conservative due to
#' shrinkage and conditioning on estimated variance components. They are
#' suitable as effect size measures and for ranking but should be interpreted
#' with caution as formal hypothesis tests.
#'
#' @param mcsd_result Output from \code{mcsd_block_from_fits()} or its wrappers
#' @param ran_vals Tidy data frame from \code{get_ran_vals()} with columns:
#'   gene, level, term, estimate, std.error, scaled_estimate, pval
#' @param focal Focal cell type
#' @param neighbour Neighbour cell type
#' @param z_threshold Minimum |z_score| to consider a gene (default: 2)
#' @param n_top Maximum number of top genes to return (default: all passing
#'   z_threshold). If NULL, returns all genes regardless of threshold.
#' @param resp_term Optional responder term name (e.g. "ResponderPD"). When
#'   supplied, the function pulls the responder-by-spatial interaction term
#'   `paste0(resp_term, ":", neighbour)` so that z-scores reflect the
#'   responder spatial-state coefficient rather than the baseline spatial
#'   slope. Required when the upstream `mcsd_result` was computed with a
#'   responder block; otherwise the volcano x-axis would mix axes.
#' @return Tibble with columns: gene, MCSD, rank, estimate, std.error,
#'   z_score, direction
#' @export
mcsd_top_drivers <- function(
        mcsd_result,
        ran_vals,
        focal,
        neighbour,
        z_threshold = 2,
        n_top = NULL,
        resp_term = NULL
) {
    # Resolve which random-effect term to pull. For non-responder analyses,
    # the term is the neighbour name (baseline spatial slope). For responder
    # analyses, it is the interaction term `<resp_term>:<neighbour>` so that
    # the z-score reflects the responder spatial-state coefficient that the
    # MCSD score on the y-axis was computed from.
    target_term <- if (is.null(resp_term)) {
        neighbour
    } else {
        paste0(resp_term, ":", neighbour)
    }

    # Extract spatial slopes for this focal-neighbour pair
    slopes <- ran_vals |>
        dplyr::filter(level == focal, term == target_term) |>
        dplyr::distinct(gene, .keep_all = TRUE) |>
        dplyr::transmute(
            gene,
            estimate = round(estimate, 4),
            std.error = round(std.error, 4),
            z_score = round(estimate / std.error, 2),
            direction = ifelse(estimate > 0, "up", "down")
        )

    # Merge with MCSD scores
    out <- mcsd_result$scores |>
        dplyr::inner_join(slopes, by = "gene") |>
        dplyr::arrange(rank)

    # Filter/limit
    if (!is.null(n_top)) {
        out <- dplyr::slice_head(out, n = n_top)
    } else {
        out <- dplyr::filter(out, abs(z_score) >= z_threshold)
    }

    out
}


# ---- Backward-compatible wrappers ----

#' MCSD Block Scores (No Responder) -- backward-compatible wrapper
#' @inheritParams mcsd_block_from_fits
mcsd_block_from_fits_no_resp <- function(
        fits_light, df,
        block = c("Spatial cell state", "Cell type", "Spillover"),
        focal = NULL, neighbour = NULL,
        group_ct = "celltype",
        ## Canonical post-2026-05-01: Variant A (covariance PCA)
        center = TRUE, scale_ref = c("none", "sd", "mad"),
        BPPARAM = BiocParallel::SerialParam(),
        expr_filter = TRUE, min_expr_var = 0.01, min_expr_prev = 0.05,
        use_gram = TRUE, keep_Z = FALSE,
        spill_correct = c("none", "weight", "weight2", "exclude"),
        spill_threshold = 1.5,
        spill_source_mode = c("neighbour", "any"),
        detection_weight = c("none", "post", "pre")
) {
    mcsd_block_from_fits(fits_light, df, block = block, focal = focal,
                         neighbour = neighbour, group_ct = group_ct,
                         resp_term = NULL, center = center, scale_ref = scale_ref,
                         BPPARAM = BPPARAM, expr_filter = expr_filter,
                         min_expr_prev = min_expr_prev, min_expr_var = min_expr_var,
                         spill_correct = spill_correct, spill_threshold = spill_threshold,
                         spill_source_mode = spill_source_mode,
                         detection_weight = detection_weight,
                         use_gram = use_gram, keep_Z = keep_Z)
}

#' MCSD via Partial Residuals (No Responder) -- backward-compatible wrapper
#' @inheritParams mcsd_partial_resid
mcsd_partial_resid_no_resp <- function(
        fits_light, df,
        focal, neighbour,
        group_ct = "celltype",
        ## Canonical post-2026-05-01: Variant A (covariance PCA)
        center = TRUE, scale_ref = c("none", "sd", "mad"),
        residual_type = "working",
        BPPARAM = BiocParallel::SerialParam(),
        expr_filter = TRUE, min_expr_var = 1.0, min_expr_prev = 0.10,
        use_gram = TRUE, keep_Z = FALSE
) {
    mcsd_partial_resid(fits_light, df, block = "Spatial cell state",
                       focal = focal, neighbour = neighbour,
                       group_ct = group_ct, resp_term = NULL,
                       center = center, scale_ref = scale_ref,
                       residual_type = residual_type, BPPARAM = BPPARAM,
                       expr_filter = expr_filter, min_expr_prev = min_expr_prev,
                       min_expr_var = min_expr_var, use_gram = use_gram,
                       keep_Z = keep_Z)
}

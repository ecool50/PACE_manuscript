#' ============================================================================
#' DenoIST-style Poisson Mixture Model for Transcript Contamination Correction
#' ============================================================================
#'
#' Fully vectorized R implementation — no cell×gene inner loops.
#' EM steps use matrix operations throughout.
#'
#' Dependencies: Matrix, FNN, data.table
#' ============================================================================

required_pkgs <- c("Matrix", "FNN", "data.table")
for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop(sprintf("Package '%s' required. Install: install.packages('%s')", pkg, pkg))
}

library(Matrix)
library(FNN)
library(data.table)

# =============================================================================
# STEP 1: Compute contamination profiles (vectorized with sparse neighbor mat)
# =============================================================================

compute_contamination_profiles <- function(counts, coords, cell_types,
                                           radius = 30, lambda_decay = 10) {
    
    n_cells <- nrow(counts)
    message(sprintf("Computing contamination profiles for %d cells...", n_cells))
    
    k_max <- min(50, n_cells - 1)
    nn <- FNN::get.knnx(coords, coords, k = k_max)
    message("  kNN search complete. Building sparse weight matrix...")
    
    verbose_internal <- TRUE
    
    # Build sparse weight matrix: only different-type neighbors within radius
    # Accumulate triplets for sparse construction
    from_list <- vector("list", n_cells)
    to_list   <- vector("list", n_cells)
    w_list    <- vector("list", n_cells)
    
    for (i in seq_len(n_cells)) {
        if (verbose_internal && i %% 10000 == 0)
            message(sprintf("    Building neighbors: %d / %d cells...", i, n_cells))
        ids   <- nn$nn.index[i, ]
        dists <- nn$nn.dist[i, ]
        keep  <- dists > 0 & dists <= radius & cell_types[ids] != cell_types[i]
        if (!any(keep)) next
        
        ids_k   <- ids[keep]
        dists_k <- dists[keep]
        wts     <- exp(-dists_k / lambda_decay)
        wts     <- wts / sum(wts)
        
        from_list[[i]] <- rep(i, length(ids_k))
        to_list[[i]]   <- ids_k
        w_list[[i]]    <- wts
    }
    
    ii <- unlist(from_list)
    jj <- unlist(to_list)
    ww <- unlist(w_list)
    
    # Sparse weight matrix × count matrix = contamination profiles
    W_sparse <- Matrix::sparseMatrix(i = ii, j = jj, x = ww,
                                     dims = c(n_cells, n_cells),
                                     repr = "C")  # force dgCMatrix
    contam_profiles <- as.matrix(W_sparse %*% as.matrix(counts))
    colnames(contam_profiles) <- colnames(counts)
    
    message("  Done.")
    contam_profiles
}

# =============================================================================
# STEP 2: Cell-type profiles
# =============================================================================

compute_celltype_profiles <- function(counts, cell_types) {
    types <- unique(cell_types)
    profiles <- matrix(0, length(types), ncol(counts))
    rownames(profiles) <- types
    colnames(profiles) <- colnames(counts)
    
    for (ct in types) {
        idx <- which(cell_types == ct)
        profiles[ct, ] <- if (length(idx) == 1) as.numeric(counts[idx, ])
        else colMeans(as.matrix(counts[idx, , drop = FALSE]))
    }
    profiles
}

# =============================================================================
# STEP 3: Fully vectorized EM
# =============================================================================

fit_denoist_em <- function(counts, cell_types, contam_profiles,
                           max_iter = 30, tol = 1e-4, verbose = TRUE) {
    
    n_cells <- nrow(counts)
    n_genes <- ncol(counts)
    counts_mat <- as.matrix(counts)
    
    ct_profiles <- compute_celltype_profiles(counts, cell_types)
    
    # Size factors
    lib <- rowSums(counts_mat); lib[lib == 0] <- 1
    sf <- lib / median(lib)
    
    # Initialise lambda_true from cell-type profiles × size factor
    lambda_true <- ct_profiles[cell_types, , drop = FALSE] * sf
    lambda_true <- pmax(lambda_true, 1e-10)
    
    # Initialise lambda_contam
    lambda_contam <- pmax(contam_profiles * 0.1, 1e-10)
    
    # Per-cell mixing proportion (start 90% true)
    pi_true <- rep(0.9, n_cells)
    
    prev_ll <- -Inf
    
    for (iter in seq_len(max_iter)) {
        
        # === E-step (fully vectorized, numerically stable) ===
        # Clamp lambdas to prevent log(0) and overflow
        lt <- pmin(pmax(lambda_true, 1e-10), 1e6)
        lc <- pmin(pmax(lambda_contam, 1e-10), 1e6)
        
        log_p_true   <- counts_mat * log(lt) - lt
        log_p_contam <- counts_mat * log(lc) - lc
        
        log_wt_true   <- log_p_true   + log(pi_true + 1e-10)
        log_wt_contam <- log_p_contam + log(1 - pi_true + 1e-10)
        
        # Stable softmax
        log_max <- pmax(log_wt_true, log_wt_contam)
        # Guard against -Inf - (-Inf) = NaN
        finite_mask <- is.finite(log_max)
        exp_true   <- exp(log_wt_true   - log_max)
        exp_contam <- exp(log_wt_contam - log_max)
        denom_post <- exp_true + exp_contam
        denom_post[denom_post == 0] <- 1  # avoid 0/0
        
        posterior_true <- exp_true / denom_post
        posterior_true[!finite_mask] <- 0.5
        posterior_true[is.nan(posterior_true)] <- 0.5
        posterior_true <- pmax(pmin(posterior_true, 1 - 1e-10), 1e-10)
        
        # === M-step (vectorized per cell type) ===
        pi_true <- rowMeans(posterior_true)
        pi_true <- pmax(pmin(pi_true, 1 - 1e-6), 1e-6)
        
        # Update cell-type profiles
        for (ct in unique(cell_types)) {
            idx <- which(cell_types == ct)
            if (length(idx) < 2) next
            
            wt_counts <- counts_mat[idx, , drop = FALSE] * posterior_true[idx, , drop = FALSE]
            denom <- colSums(posterior_true[idx, , drop = FALSE] * sf[idx])
            denom[denom == 0] <- 1e-10
            updated <- colSums(wt_counts) / denom
            
            lambda_true[idx, ] <- outer(sf[idx], updated)
            ct_profiles[ct, ] <- updated
        }
        
        # Update contamination rates (vectorized)
        contam_wt_sum <- rowSums((1 - posterior_true) * counts_mat)
        contam_profile_sum <- rowSums(contam_profiles * sf + 1e-10)
        scale_fac <- contam_wt_sum / contam_profile_sum
        scale_fac[!is.finite(scale_fac)] <- 0
        lambda_contam <- pmax(contam_profiles * sf * scale_fac, 1e-10)
        
        lambda_true   <- pmax(lambda_true, 1e-10)
        lambda_contam <- pmax(lambda_contam, 1e-10)
        
        # Convergence — use only finite values
        ll_terms <- log_wt_true * posterior_true + log_wt_contam * (1 - posterior_true)
        ll <- sum(ll_terms[is.finite(ll_terms)])
        
        if (verbose)
            message(sprintf("  EM iter %d: log-lik = %.2f, mean(π) = %.4f", iter, ll, mean(pi_true)))
        
        if (is.finite(ll) && is.finite(prev_ll) &&
            abs(ll - prev_ll) / (abs(prev_ll) + 1) < tol) {
            if (verbose) message(sprintf("  Converged at iteration %d", iter))
            break
        }
        prev_ll <- ll
    }
    
    list(pi = pi_true, lambda_true = lambda_true, lambda_contam = lambda_contam,
         posterior_true = posterior_true, ct_profiles = ct_profiles)
}

# =============================================================================
# STEP 4: Corrected counts (vectorized)
# =============================================================================

generate_corrected_counts <- function(counts, posterior_true,
                                      method = "expected", seed = 42) {
    counts_mat <- as.matrix(counts)
    
    corrected <- if (method == "expected") {
        round(counts_mat * posterior_true)
    } else if (method == "binomial") {
        set.seed(seed)
        matrix(rbinom(length(counts_mat), size = as.integer(counts_mat),
                      prob = as.numeric(posterior_true)),
               nrow = nrow(counts_mat), ncol = ncol(counts_mat),
               dimnames = dimnames(counts_mat))
    } else stop("method must be 'expected' or 'binomial'")
    
    corrected <- pmax(corrected, 0L)
    storage.mode(corrected) <- "integer"
    corrected
}

# =============================================================================
# MAIN WRAPPER
# =============================================================================

run_denoist_correction <- function(counts, coords, cell_types,
                                   radius = 30, lambda_decay = 10,
                                   max_iter = 30, method = "expected",
                                   verbose = TRUE) {
    
    if (verbose) message("=== DenoIST Poisson Mixture Correction ===")
    
    contam_profiles <- compute_contamination_profiles(
        counts, coords, cell_types, radius = radius, lambda_decay = lambda_decay)
    
    if (verbose) message("Fitting Poisson mixture model via EM...")
    em <- fit_denoist_em(counts, cell_types, contam_profiles,
                         max_iter = max_iter, verbose = verbose)
    
    if (verbose) message("Generating corrected integer counts...")
    corrected <- generate_corrected_counts(counts, em$posterior_true, method = method)
    
    total_orig <- sum(counts)
    total_corr <- sum(corrected)
    
    stats <- do.call(rbind, lapply(unique(cell_types), function(ct) {
        idx <- which(cell_types == ct)
        data.frame(cell_type = ct, n_cells = length(idx),
                   mean_pi = mean(em$pi[idx]),
                   original_counts = sum(counts[idx, ]),
                   corrected_counts = sum(corrected[idx, ]),
                   pct_removed = 100 * (1 - sum(corrected[idx, ]) / sum(counts[idx, ])),
                   stringsAsFactors = FALSE)
    }))
    
    if (verbose) {
        message(sprintf("\nOverall: removed %.1f%% of transcripts",
                        100 * (1 - total_corr / total_orig)))
        message("\nPer-cell-type summary:"); print(stats)
    }
    
    list(corrected_counts = corrected, pi = em$pi,
         posterior_true = em$posterior_true, contam_profiles = contam_profiles,
         ct_profiles = em$ct_profiles, stats = stats)
}

# =============================================================================
# SPE Helper
# =============================================================================

correct_spe_denoist <- function(spe, celltype_col = "celltype",
                                radius = 30, lambda_decay = 10) {
    counts_mat <- t(as.matrix(assay(spe, "counts")))
    coords_mat <- as.matrix(spatialCoords(spe))
    cell_types <- colData(spe)[[celltype_col]]
    
    result <- run_denoist_correction(counts = counts_mat, coords = coords_mat,
                                     cell_types = cell_types,
                                     radius = radius, lambda_decay = lambda_decay)
    
    assay(spe, "corrected") <- t(result$corrected_counts)
    metadata(spe)$denoist_stats <- result$stats
    metadata(spe)$denoist_pi    <- result$pi
    spe
}
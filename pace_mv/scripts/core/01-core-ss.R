# ==============================================================================
# 01-core-ss.R
# Utilities and core sum-of-squares computation for variance decomposition
# Works natively with fit_light objects from 00-extract.R
# ==============================================================================

## ============================== UTILITIES ===================================
## Helper to ensure we never spawn nested pools
.bp_or_serial <- function(BPPARAM) {
    if (is.null(BPPARAM)) return(BiocParallel::SerialParam())
    BPPARAM
}

## Check if object is a fit_light
.is_fit_light <- function(x) {
    inherits(x, "fit_light") || 
        (is.list(x) && all(c("used_rows", "fixef", "ranef") %in% names(x)))
}

## ====================== COEF SS: CORE COMPUTATIONS ==========================

# ---- Fixed-effects SS: vectorized ----
# resp_term/resp_level let this function recognise factor-coded fixed effects
# such as `ResponderPD` (i.e., when df has column `Responder` and the model
# has a fixed effect named `<resp_term><resp_level>`). Without this, the new
# fixed Responder term in the responder-aware model would be silently dropped
# from variance decomposition.
.fixed_ss <- function(beta, mask_idx, df_loc, resp_term = NULL, resp_level = NULL) {
    keep <- intersect(names(beta), c("(Intercept)", names(df_loc)))

    resp_coef <- if (!is.null(resp_term) && !is.null(resp_level)) {
        paste0(resp_term, resp_level)
    } else NULL
    if (!is.null(resp_coef) && resp_coef %in% names(beta) &&
        !resp_coef %in% keep) {
        keep <- c(keep, resp_coef)
    }

    if (!length(keep)) return(tibble())
    n <- nrow(df_loc)
    X <- matrix(NA_real_, nrow = n, ncol = length(keep), dimnames = list(NULL, keep))
    for (k in keep) {
        X[, k] <- if (k == "(Intercept)") {
            1
        } else if (!is.null(resp_coef) && identical(k, resp_coef)) {
            as.numeric(df_loc[[resp_term]] == resp_level)
        } else {
            as.numeric(df_loc[[k]])
        }
    }
    Xk <- X[mask_idx, , drop = FALSE]
    b2 <- (beta[keep])^2
    SS <- as.numeric(colSums(Xk * Xk, na.rm = TRUE)) * b2
    tibble(
        type = "fixed", group = NA_character_, level = NA_character_,
        term = keep, coef_id = paste0("fixed:", keep), SS_link = SS
    )
}

# ---- Random-effects SS by z^2 group sum (unified) ----
# When resp_term is NULL, uses no-responder logic.
# When resp_term is not NULL, computes responder indicator and uses resp-aware logic.
.ranef_ss <- function(df_loc, mask_idx, re_list, group_terms, resp_term = NULL, resp_level = NULL) {
    out <- list()
    if (!is.null(resp_term)) {
        resp_ind_full <- as.numeric(df_loc$Responder == resp_level)
    }
    for (grp in intersect(names(re_list), group_terms)) {
        tab   <- re_list[[grp]]
        levs  <- rownames(tab); terms <- colnames(tab)
        idx_all <- match(df_loc[[grp]], levs)
        idx_m   <- idx_all[mask_idx]
        nlev    <- length(levs)

        z2_by_level <- function(tm) {
            if (identical(tm, "(Intercept)")) {
                cnt <- suppressWarnings(rowsum(rep.int(1, length(idx_m)), idx_m, reorder = FALSE))
                S <- numeric(nlev); S[as.integer(rownames(cnt))] <- as.numeric(cnt[,1]); return(S)
            }
            if (!is.null(resp_term)) {
                # Responder-aware path
                z_m <- if (tm == resp_term) {
                    resp_ind_full[mask_idx]
                } else if (grepl(":", tm, fixed = TRUE)) {
                    parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
                    z <- rep.int(1, length(idx_m))
                    for (p in parts) {
                        if (p == resp_term) {
                            z <- z * resp_ind_full[mask_idx]
                        } else if (p %in% names(df_loc)) {
                            z <- z * as.numeric(df_loc[[p]][mask_idx])
                        } else { z <- 0; break }
                    }
                    z
                } else if (tm %in% names(df_loc)) {
                    as.numeric(df_loc[[tm]][mask_idx])
                } else 0
            } else {
                # No-responder path
                z_m <- if (tm %in% names(df_loc)) {
                    as.numeric(df_loc[[tm]][mask_idx])
                } else if (grepl(":", tm, fixed = TRUE)) {
                    parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
                    z <- rep.int(1, length(idx_m))
                    for (p in parts) {
                        if (p %in% names(df_loc)) z <- z * as.numeric(df_loc[[p]][mask_idx]) else { z <- 0; break }
                    }
                    z
                } else 0
            }
            z2 <- if (is.numeric(z_m) && length(z_m) == length(idx_m)) z_m * z_m else rep(0, length(idx_m))
            srs <- suppressWarnings(rowsum(z2, idx_m, reorder = FALSE))
            S <- numeric(nlev); if (length(srs)) S[as.integer(rownames(srs))] <- as.numeric(srs[,1]); S
        }

        if (!length(terms)) next
        S_mat <- do.call(cbind, lapply(terms, z2_by_level))
        SS_mat <- as.matrix(tab)^2 * S_mat
        keep <- SS_mat > 0 & is.finite(SS_mat)
        if (!any(keep)) next
        kk <- which(keep, arr.ind = TRUE)
        out[[length(out)+1L]] <- tibble(
            type="ranef", group=grp,
            level=levs[kk[,1]], term=terms[kk[,2]],
            coef_id=paste0("ranef:",grp,"[",levs[kk[,1]],"]:",terms[kk[,2]]),
            SS_link=SS_mat[keep]
        )
    }
    if (length(out)) bind_rows(out) else tibble()
}

# ---- Residual SS (masked) - for fit_light ----
.resid_ss_light <- function(residuals_list, residual, mask_idx) {
    if (identical(residual, "none")) return(tibble())
    r <- residuals_list[[residual]]
    if (is.null(r)) stop("Residual type '", residual, "' not pre-computed in fit_light.")
    tibble(
        type="residual", group=NA_character_, level=NA_character_,
        term=residual, coef_id=paste0("residual:", residual),
        SS_link=sum(r[mask_idx]^2, na.rm=TRUE)
    )
}

## ====================== PUBLIC INTERFACE ====================================

#' Coefficient Sum-of-Squares
#'
#' Unified interface for both responder-aware and no-responder analysis.
#' When \code{resp_term} is NULL, no responder logic is applied.
#' When \code{resp_term} is not NULL, computes responder-aware random effects.
#'
#' @param fit_light A fit_light object from extractRandomEffectsNew
#' @param df Original data frame
#' @param focal Optional focal cell type to filter
#' @param group_terms Random effect grouping terms
#' @param resp_term Responder term name (NULL for no-responder analysis)
#' @param resp_level Responder level to use (ignored when resp_term is NULL)
#' @param residual Residual type
#' @return Tibble of SS contributions by coefficient
coef_ss_link <- function(
        fit_light, df,
        focal = NULL,
        group_terms = c("celltype", "imageID"),
        resp_term   = NULL,
        resp_level  = "PD",
        residual    = c("working", "none", "response", "pearson", "deviance")
) {
    residual <- match.arg(residual)
    stopifnot(.is_fit_light(fit_light))

    used_rows <- fit_light$used_rows
    df_loc <- droplevels(df[used_rows, , drop = FALSE])
    n <- nrow(df_loc)

    mask_idx <- if (!is.null(focal)) {
        m <- which(df_loc$celltype == focal)
        if (!length(m)) stop("No rows for focal = '", focal, "'.")
        m
    } else seq_len(n)

    # Fixed effects
    beta <- fit_light$fixef
    out_fix <- .fixed_ss(beta, mask_idx, df_loc,
                         resp_term = resp_term, resp_level = resp_level)

    # Random effects
    re_list <- fit_light$ranef
    out_re <- .ranef_ss(df_loc, mask_idx, re_list, group_terms,
                        resp_term = resp_term, resp_level = resp_level)

    # Residuals
    out_res <- .resid_ss_light(fit_light$residuals, residual, mask_idx)

    bind_rows(out_fix, out_re, out_res) |>
        mutate(n_obs = length(mask_idx)) |>
        arrange(desc(.data$SS_link))
}

#' Coefficient Sum-of-Squares (No Responder) - backward-compatible wrapper
#'
#' @param fit_light A fit_light object from extractRandomEffectsNew
#' @param df Original data frame
#' @param focal Optional focal cell type to filter
#' @param group_terms Random effect grouping terms
#' @param residual Residual type
#' @return Tibble of SS contributions by coefficient
coef_ss_link_no_resp <- function(fit_light, df, focal = NULL,
                                  group_terms = c("celltype", "imageID"),
                                  residual = c("working", "none", "response", "pearson", "deviance")) {
    coef_ss_link(fit_light, df, focal = focal, group_terms = group_terms,
                 resp_term = NULL, residual = residual)
}
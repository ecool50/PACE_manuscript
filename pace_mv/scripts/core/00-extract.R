# ==============================================================================
# 00-extract.R
# Model fitting with lightweight extraction for variance decomposition pipeline
#
# Canonical post-sensitivity model (locked 2026-05-01):
#   family = nbinom1()                 — wins per-gene AIC vs NB2 on flagship
#                                         drivers (FN1, SPP1, COL1A1, COL6A2;
#                                         see test_nb2_top_drivers.R) and by
#                                         ~12k AIC vs Poisson (sensitivity_poisson.R)
#   RE structure = `||` (uncorrelated) — full `|` correlation is unidentifiable
#                                         with K=6 cell-type levels (115 RE
#                                         params, all 6/6 fits Hessian non-PD;
#                                         see test_nb1_correlated_re.R)
# A cell-complexity term scale(log1p(n_features_minus_g)) was tested in
# test_nb1_nfeatures.R; it improves per-gene AIC by 4-126 on flagship drivers
# without changing their BLUPs, but is intentionally NOT included in the
# canonical model — the offset(log(nCount)) library-size correction was kept
# as the only per-cell normalisation.
# Do not change the family or RE operator without re-running the corresponding
# sensitivity scripts.
# ==============================================================================

#' Extract Random Effects with Lightweight Fit Objects (Unified)
#'
#' Fits glmmTMB models for each gene and returns both tidy random effects
#' and lightweight fit objects suitable for variance decomposition.
#'
#' When \code{resp_var} is NULL (single-sample mode):
#'   gene ~ 1 + offset(log(nCount)) + spill_terms + (1 + vars || celltype)
#'
#' When \code{resp_var} is not NULL (multi-sample responder mode):
#'   gene ~ 1 + offset(log(nCount)) + resp_var + spill_terms +
#'          (1 + resp_var * (vars) || celltype) + (1 || imageID)
#'
#' The fixed \code{resp_var} term stabilises the random per-celltype Responder
#' effects. Without it, the random Responder BLUPs are partially confounded
#' with imageID random intercepts (Responder is constant within each image),
#' which inflates the BLUPs and makes them uninterpretable in isolation.
#' The fixed term absorbs the population-mean Responder effect; the random
#' per-celltype Responder slopes then represent well-defined deviations from it.
#'
#' The \code{offset(log(nCount))} term provides per-cell library-size correction
#' so that the conditional mean is modelled as a rate (counts per total
#' transcript output) rather than a raw count, the standard approach for NB
#' GLMs on single-cell count data. Setting \code{use_offset = FALSE} reverts
#' to the unadjusted model.
#'
#' @param genes Character vector of gene names
#' @param df Data frame with cell-level covariates. Must contain gene columns,
#'   celltype factor, neighbourhood abundance columns (vars), and spillover
#'   columns (vars_spill or auto-derived as paste0(vars, "_spill")).
#'   When resp_var is not NULL, must also contain resp_var and imageID columns.
#'   When use_offset is TRUE, must also contain an nCount column with per-cell
#'   total transcript counts.
#' @param vars Character vector of covariate names for random slopes
#'   (typically cell type names matching neighbourhood columns)
#' @param vars_spill Character vector of spillover column names.
#'   If NULL, derived as paste0(vars, "_spill")
#' @param include_spillover Logical; include spillover terms in model
#' @param resp_var Name of the responder column in df, or NULL for
#'   single-sample mode (default: NULL)
#' @param resp_level Active level of responder (default: "PD")
#' @param family glmmTMB family (default: nbinom1)
#' @param BPPARAM BiocParallel backend
#' @param residual_types Character vector of residual types to pre-compute
#' @param keep_fit Logical; also return full fit object (memory intensive)
#' @param use_offset Logical; include offset(log(nCount)) for library-size
#'   correction (default: TRUE). When TRUE, df must contain an nCount column.
#' @param area_var Optional column name in df giving a per-cell physical size
#'   measurement (e.g. segmentation area in µm^2). When non-NULL, adds
#'   scale(log(<area_var>)) as a fixed effect, which provides a true cell-size
#'   correction beyond the library-size offset. Default NULL (no area term).
#'
#' @return Named list (by gene) of lists containing:
#'   - ran_vals: tidy data frame of random effects
#'   - fit_light: lightweight list for variance decomposition
#'   - gene: gene name
#'   - error: logical
extractRandomEffects <- function(
        genes, df, vars,
        vars_spill = NULL,
        include_spillover = TRUE,
        resp_var = NULL,
        resp_level = "PD",
        family = glmmTMB::nbinom1(),
        BPPARAM = BiocParallel::MulticoreParam(workers = 4),
        residual_types = "working",
        keep_fit = FALSE,
        use_offset = TRUE,
        area_var = NULL,
        image_slope_vars = NULL
) {

    if (is.null(vars_spill)) vars_spill <- paste0(vars, "_spill")

    has_responder <- !is.null(resp_var)
    has_area <- !is.null(area_var)
    if (has_area) {
        if (!area_var %in% colnames(df)) {
            stop(sprintf("area_var = '%s' not found in df columns.", area_var))
        }
        if (any(df[[area_var]] <= 0, na.rm = TRUE)) {
            stop(sprintf("%s must be strictly positive for log scaling.", area_var))
        }
        df$.area_z <- as.numeric(scale(log(df[[area_var]])))
    }

    # Resolve image_slope_vars:
    #   NULL or FALSE -> no per-image random slopes (intercept only)
    #   TRUE          -> per-image random slopes for ALL neighbour vars
    #   character     -> per-image random slopes for the named subset
    if (isTRUE(image_slope_vars)) {
        image_slope_vars <- vars
    } else if (isFALSE(image_slope_vars)) {
        image_slope_vars <- NULL
    }

    if (has_responder) {
        stopifnot(resp_var %in% names(df))
        stopifnot("imageID" %in% names(df))
        stopifnot("celltype" %in% names(df))
        if (!is.null(image_slope_vars)) {
            missing_vars <- setdiff(image_slope_vars, names(df))
            if (length(missing_vars))
                stop("image_slope_vars not found in df: ",
                     paste(missing_vars, collapse = ", "))
        }
    }
    if (use_offset) {
        if (!"nCount" %in% names(df)) {
            stop("use_offset = TRUE requires an `nCount` column in df ",
                 "(per-cell total transcript counts). Either compute it ",
                 "before calling, or pass use_offset = FALSE.")
        }
        if (any(df$nCount <= 0, na.rm = TRUE)) {
            stop("nCount must be strictly positive for log-offset. ",
                 "Filter zero-count cells first.")
        }
    }

    BiocParallel::bplapply(genes, function(gene) {

        # Fixed effects: intercept + optional library-size offset +
        #                fixed Responder (multi-sample only) + spillover
        # Fixed Responder stabilises the random Responder BLUPs:
        # without it, the random Responder per-celltype effects are partially
        # confounded with imageID random intercepts (Responder is constant
        # within image), leading to inflated BLUPs. Adding it as fixed gives
        # the random per-celltype Responder effects a well-defined zero mean.
        parts <- c("1")
        if (use_offset) parts <- c(parts, "offset(log(nCount))")
        if (has_responder) parts <- c(parts, resp_var)
        if (include_spillover) {
            parts <- c(parts, paste(vars_spill, collapse = " + "))
        }
        if (has_area) parts <- c(parts, ".area_z")
        cond_rhs <- paste(parts, collapse = " + ")

        # Random effects
        # image_slope_vars: NULL/FALSE -> intercept-only image RE (default);
        #   TRUE -> per-image random slopes for ALL neighbour vars (recommended
        #   for cross-patient robustness, since the per-celltype spatial
        #   differentials should not depend on a few outlier patients);
        #   character vector -> per-image slopes for the named subset.
        # Adding image random slopes produces large AIC improvements in
        # melanoma flagship genes (-150 to -600 per gene) and prevents
        # per-patient slope variation from being absorbed into the fixed/
        # random Responder interactions. Uses `||` (uncorrelated) for
        # stability when many slope terms are added.
        if (has_responder) {
            vars_str <- paste(vars, collapse = " + ")
            rand_ct  <- glue::glue("(1 + {resp_var} * ({vars_str}) || celltype)")
            if (!is.null(image_slope_vars) && length(image_slope_vars)) {
                slope_str <- paste(image_slope_vars, collapse = " + ")
                rand_img <- glue::glue("(1 + {slope_str} || imageID)")
            } else {
                rand_img <- "(1 || imageID)"
            }
            cond_formula <- stats::as.formula(
                glue::glue("`{gene}` ~ {cond_rhs} + {rand_ct} + {rand_img}")
            )
        } else {
            rand_ct <- glue::glue("(1 + {paste(vars, collapse = ' + ')} || celltype)")
            cond_formula <- stats::as.formula(
                glue::glue("`{gene}` ~ {cond_rhs} + {rand_ct}")
            )
        }

        # Fit model
        fit_args <- list(
            formula = cond_formula,
            data    = df,
            family  = family,
            verbose = FALSE
        )
        if (has_responder) {
            fit_args$control <- glmmTMB::glmmTMBControl(parallel = 1L)
        }

        fit <- tryCatch(
            do.call(glmmTMB::glmmTMB, fit_args),
            error = function(e) NULL
        )

        if (is.null(fit)) {
            return(list(ran_vals = NULL, fit_light = NULL, gene = gene, error = TRUE))
        }

        ran_vals <- broom.mixed::tidy(fit, effects = "ran_vals") |>
            dplyr::mutate(
                gene            = gene,
                lower           = estimate - 2 * std.error,
                upper           = estimate + 2 * std.error,
                scaled_estimate = estimate / std.error,
                pval            = 2 * stats::pnorm(-abs(scaled_estimate))
            )

        mf <- stats::model.frame(fit)
        used_rows <- rownames(mf)
        if (is.null(used_rows)) used_rows <- as.character(seq_len(nrow(mf)))

        # pdHess captures whether the random-effect Hessian is positive-definite.
        # When FALSE, the variance components hit the boundary and per-BLUP
        # standard errors are unreliable (broom.mixed::tidy returns NA).
        # Stored here so downstream code (MCSD, plotting) can filter
        # boundary-singular genes without re-deriving from ran_vals SE.
        pd_hess <- tryCatch(isTRUE(fit$sdr$pdHess), error = function(e) NA)

        # Components needed for Nakagawa-Johnson-Schielzeth variance
        # decomposition (per pace_variance_decomposition() in helpers/):
        #   sigma_disp   = NB1 dispersion (alpha) for the residual term
        #                  log(1 + 1/mu_bar + alpha_NB1)
        #   varcorr      = random-effect covariance matrices, named by
        #                  group, used to read off sigma^2_alpha and the
        #                  per-slope sigma^2_gammak directly (not from
        #                  shrunken BLUPs)
        #   mu_summary   = small summary of fitted(type="response") used to
        #                  compute the residual term without storing the
        #                  full per-cell vector
        sigma_disp  <- tryCatch(glmmTMB::sigma(fit), error = function(e) NA_real_)
        varcorr_cond <- tryCatch(glmmTMB::VarCorr(fit)$cond,
                                  error = function(e) NULL)
        mu_summary <- tryCatch({
          mu <- as.numeric(stats::fitted(fit, type = "response"))
          ct_in_mf <- if ("celltype" %in% colnames(mf)) {
            as.character(mf[["celltype"]])
          } else NULL
          list(mean = mean(mu, na.rm = TRUE),
               by_celltype = if (!is.null(ct_in_mf))
                              tapply(mu, ct_in_mf, mean, na.rm = TRUE)
                              else NULL)
        }, error = function(e) NULL)

        fit_light <- list(
            gene         = gene,
            used_rows    = used_rows,
            n_obs        = length(used_rows),
            fixef        = glmmTMB::fixef(fit)$cond,
            ranef        = lme4::ranef(fit)$cond,
            pdHess       = pd_hess,
            sigma_disp   = sigma_disp,
            varcorr      = varcorr_cond,
            mu_summary   = mu_summary,
            residuals    = lapply(
                stats::setNames(residual_types, residual_types),
                function(rt) as.numeric(stats::residuals(fit, type = rt))
            )
        )
        class(fit_light) <- c("fit_light", "list")

        result <- list(
            ran_vals  = ran_vals,
            fit_light = fit_light,
            gene      = gene,
            error     = FALSE
        )

        if (keep_fit) result$fit <- fit

        result

    }, BPPARAM = BPPARAM) |>
        stats::setNames(genes)
}


# ==============================================================================
# Backward-compatible wrappers
# ==============================================================================

#' @describeIn extractRandomEffects Single-sample wrapper (resp_var = NULL)
extractRandomEffectsNew <- function(genes, df, vars, ...) {
    extractRandomEffects(genes, df, vars, resp_var = NULL, ...)
}

#' @describeIn extractRandomEffects Multi-sample responder wrapper
extractRandomEffectsResponder <- function(
        genes, df, vars,
        vars_spill = NULL,
        include_spillover = TRUE,
        resp_var = "Responder",
        resp_level = "PD",
        family = glmmTMB::nbinom1(),
        BPPARAM = BiocParallel::MulticoreParam(workers = 4),
        residual_types = "working",
        keep_fit = FALSE,
        use_offset = TRUE,
        area_var = NULL,
        image_slope_vars = NULL
) {
    extractRandomEffects(genes, df, vars,
                         vars_spill = vars_spill,
                         include_spillover = include_spillover,
                         resp_var = resp_var,
                         resp_level = resp_level,
                         family = family,
                         BPPARAM = BPPARAM,
                         residual_types = residual_types,
                         keep_fit = keep_fit,
                         use_offset = use_offset,
                         area_var = area_var,
                         image_slope_vars = image_slope_vars)
}


# ==============================================================================
# Helper extractors
# ==============================================================================

#' Identify well-fit (non-boundary-singular) genes
#'
#' Returns a logical vector indicating which genes have a numerically stable
#' fit. A gene is considered well-fit if either:
#'   - its `fit_light$pdHess` is TRUE (variance-component Hessian is
#'     positive-definite; preferred signal, available for fits made after
#'     2026-04-30), OR
#'   - at least one of its `ran_vals$std.error` values is finite (fallback
#'     for older saved results that pre-date the `pdHess` field).
#'
#' Genes with all-NA standard errors and `pdHess = FALSE` are boundary-
#' singular: their variance components collapsed to zero, the BLUPs are
#' uninterpretable, and they should be excluded from downstream aggregations
#' (MCSD, variance pies, summary tables) where they would inflate gene
#' counts without contributing real signal.
#'
#' @param results Output from extractRandomEffects (or wrappers)
#' @return Named logical vector, TRUE for well-fit genes
#' @export
is_well_fit <- function(results) {
    vapply(results, function(r) {
        if (isTRUE(r$error)) return(FALSE)
        # Prefer pdHess if recorded
        pd <- r$fit_light$pdHess
        if (!is.null(pd) && !is.na(pd)) return(isTRUE(pd))
        # Fallback: any finite std.error in ran_vals
        if (is.null(r$ran_vals) || nrow(r$ran_vals) == 0) return(FALSE)
        any(is.finite(r$ran_vals$std.error))
    }, logical(1))
}

#' Extract fit_light objects from extraction results
#'
#' @param results Output from extractRandomEffects (or wrappers)
#' @param well_fit_only Logical; if TRUE, drop genes flagged as boundary-
#'   singular by \code{is_well_fit()}. Default FALSE for backward
#'   compatibility, but TRUE is recommended for any downstream aggregation
#'   that pools across genes (MCSD, variance decomposition, etc.).
#' @return Named list of fit_light objects
get_fits_light <- function(results, well_fit_only = FALSE) {
    if (well_fit_only) {
        ok <- is_well_fit(results)
        results <- results[ok]
    }
    out <- lapply(results, `[[`, "fit_light")
    Filter(Negate(is.null), out)
}

#' Extract ran_vals from extraction results
#'
#' @param results Output from extractRandomEffects (or wrappers)
#' @param bind Logical; bind into single data frame
#' @param well_fit_only Logical; if TRUE, drop genes flagged as boundary-
#'   singular by \code{is_well_fit()}. Default FALSE for backward
#'   compatibility.
#' @return List of data frames or single bound data frame
get_ran_vals <- function(results, bind = TRUE, well_fit_only = FALSE) {
    if (well_fit_only) {
        ok <- is_well_fit(results)
        results <- results[ok]
    }
    rv <- lapply(results, `[[`, "ran_vals")
    rv <- Filter(Negate(is.null), rv)
    if (bind) dplyr::bind_rows(rv) else rv
}

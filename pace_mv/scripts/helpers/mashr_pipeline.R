## Pipeline-step mashr shrinkage of per-celltype neighbour-slope BLUPs.
##
## Two entry points:
##   apply_mashr_shrinkage(results, focals, neighbours, resp_term = NULL)
##     -> long-form tibble (gene, focal, neighbour, estimate, std.error,
##                          estimate_shrunk, sd_shrunk, lfsr)
##     One mashr fit per neighbour slice, focals as conditions.
##
##   merge_shrunken_into_ranvals(ran_vals, shrunken_long, resp_term = NULL)
##     -> ran_vals + extra columns (estimate_shrunk, sd_shrunk, lfsr) joined
##     on (gene, level=focal, term=neighbour). Original estimate / std.error
##     untouched so downstream code can choose either.
##
## In melanoma (multi-sample) pass `resp_term = "ResponderPD"` and the
## function shrinks the ResponderPD:<neighbour> interaction BLUPs. In breast
## cancer (single-sample) leave `resp_term = NULL` and it shrinks the bare
## <neighbour> proximity-slope BLUPs.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(mashr); library(ashr)
})

#' Run mashr per-neighbour and return shrunken BLUPs in long form
#'
#' @param results PACE fit results from extractRandomEffects(*)
#' @param focals  Character vector of focal cell-type levels (the columns
#'                of each per-neighbour Bhat / Shat matrix).
#' @param neighbours Character vector of neighbour cell types whose slopes
#'                   are to be shrunk. The actual `term` looked up in
#'                   ran_vals is `paste0(resp_term, ":", neighbour)` if
#'                   resp_term is non-NULL, otherwise just `neighbour`.
#' @param resp_term Responder term name for the multi-sample model
#'                  (e.g., "ResponderPD"); NULL for single-sample.
#' @param min_significant Minimum number of "significant" genes required
#'                        before adding data-driven covariances; otherwise
#'                        canonical-only.
#' @return Long-form tibble.
apply_mashr_shrinkage <- function(results, focals, neighbours,
                                  resp_term = NULL,
                                  min_significant = 10,
                                  data_driven = TRUE) {
  ## data_driven=FALSE -> canonical covariances only (skip the slow cov_pca/cov_ed
  ## extreme-deconvolution step). Much faster when there are many neighbour terms
  ## (e.g. 16-type fits); shrinkage is slightly less adaptive but driver rankings
  ## are nearly identical. Default TRUE preserves the canonical manuscript behaviour.
  rv_all <- get_ran_vals(results, well_fit_only = TRUE)
  build_term <- function(nb) {
    if (is.null(resp_term)) nb else paste0(resp_term, ":", nb)
  }

  out <- list()
  for (nb in neighbours) {
    target_term <- build_term(nb)
    sub <- rv_all |>
      dplyr::filter(term == target_term, group == "celltype",
                    level %in% focals)
    if (nrow(sub) == 0) {
      message("  [mashr] no ran_vals for term '", target_term, "', skipping")
      next
    }

    Bhat <- sub |>
      dplyr::select(gene, focal = level, estimate) |>
      tidyr::pivot_wider(names_from = focal, values_from = estimate)
    Shat <- sub |>
      dplyr::select(gene, focal = level, std.error) |>
      tidyr::pivot_wider(names_from = focal, values_from = std.error)

    # Drop genes with NA / non-positive SE in any focal column
    cols_present <- intersect(focals, colnames(Bhat))
    if (length(cols_present) < 2) {
      message("  [mashr] term '", target_term,
              "' has <2 focal columns, skipping")
      next
    }
    genes <- Bhat$gene
    Bhat <- as.matrix(Bhat[, cols_present, drop = FALSE]); rownames(Bhat) <- genes
    Shat <- as.matrix(Shat[, cols_present, drop = FALSE]); rownames(Shat) <- genes
    ## Strict filter: require finite Bhat / Shat, positive Shat, and Shat
    ## not pathologically large (boundary-singular fits sometimes return
    ## absurd SEs that destabilise mash's likelihood matrix).
    se_cap <- 50  # any SE above this is a non-converged BLUP, drop
    se_floor <- 1e-4  # below this both Bhat and Shat are effectively zero
    keep <- rowSums(!is.finite(Bhat)) == 0 &
            rowSums(!is.finite(Shat)) == 0 &
            rowSums(Shat <= 0)        == 0 &
            rowSums(Shat > se_cap)    == 0 &
            ## drop rows where every focal's BLUP collapsed to ~zero with
            ## ~zero SE — these are degenerate boundary fits that crash
            ## mashr's likelihood-matrix check
            !(rowSums(abs(Bhat) < se_floor & Shat < se_floor) == ncol(Bhat))
    Bhat <- Bhat[keep, , drop = FALSE]
    Shat <- Shat[keep, , drop = FALSE]
    ## Floor any remaining tiny SEs to a small positive number to keep
    ## mashr's likelihood matrix well-conditioned.
    Shat <- pmax(Shat, se_floor)
    if (nrow(Bhat) < 5) {
      message("  [mashr] term '", target_term,
              "' has <5 well-fit genes after filter, skipping")
      next
    }

    data <- mash_set_data(Bhat, Shat)
    U_c  <- cov_canonical(data)
    ## Be defensive about data-driven covariances: cov_pca / cov_ed can
    ## produce non-PSD components on small / degenerate slices, which then
    ## break mash's likelihood-matrix check. Try canonical + ED; if it
    ## fails, fall back to canonical only.
    m <- tryCatch({
      m_c <- mash(data, U_c)
      top_idx <- get_significant_results(m_c)
      if (data_driven && length(top_idx) >= min_significant) {
        U_pca <- cov_pca(data, npc = min(5, length(top_idx) - 1),
                         subset = top_idx)
        U_ed  <- cov_ed(data, U_pca, subset = top_idx)
        mash(data, c(U_c, U_ed))
      } else {
        m_c
      }
    }, error = function(e) {
      message("  [mashr] data-driven cov failed for term '", target_term,
              "' (", conditionMessage(e), "); falling back to canonical only")
      tryCatch(mash(data, U_c),
               error = function(e2) {
                 message("  [mashr] canonical-only also failed (",
                         conditionMessage(e2), "); skipping")
                 NULL
               })
    })
    if (is.null(m)) next

    pm   <- get_pm(m)
    psd  <- get_psd(m)
    lfsr <- get_lfsr(m)

    out[[nb]] <- tibble::tibble(
      gene = rownames(pm)[row(pm)],
      focal = colnames(pm)[col(pm)],
      neighbour = nb,
      term      = target_term,
      estimate  = as.numeric(Bhat),
      std.error = as.numeric(Shat),
      estimate_shrunk = as.numeric(pm),
      sd_shrunk       = as.numeric(psd),
      lfsr            = as.numeric(lfsr)
    )
    message(sprintf("  [mashr] %s: %d genes shrunk; sig (lfsr<0.05) = %d",
                    target_term, nrow(Bhat),
                    length(get_significant_results(m, thresh = 0.05))))
  }

  if (length(out) == 0) {
    return(tibble::tibble(gene = character(), focal = character(),
                          neighbour = character(), term = character(),
                          estimate = numeric(), std.error = numeric(),
                          estimate_shrunk = numeric(), sd_shrunk = numeric(),
                          lfsr = numeric()))
  }
  dplyr::bind_rows(out)
}

#' Left-join shrunken values into ran_vals
#'
#' @param ran_vals  Output of get_ran_vals(results, ...).
#' @param shrunken_long Output of apply_mashr_shrinkage(...).
#' @return ran_vals plus columns estimate_shrunk, sd_shrunk, lfsr (NA where
#'   the (gene, level, term) combination was not shrunk).
merge_shrunken_into_ranvals <- function(ran_vals, shrunken_long) {
  if (nrow(shrunken_long) == 0) {
    return(ran_vals |> dplyr::mutate(estimate_shrunk = NA_real_,
                                     sd_shrunk       = NA_real_,
                                     lfsr            = NA_real_))
  }
  ran_vals |>
    dplyr::left_join(shrunken_long |>
                       dplyr::select(gene, level = focal, term,
                                     estimate_shrunk, sd_shrunk, lfsr),
                     by = c("gene", "level", "term"))
}

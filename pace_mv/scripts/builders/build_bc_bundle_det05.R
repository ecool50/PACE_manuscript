## build_bc_bundle_det05.R — reproducible BC fitting bundle with an EXPLICIT gene filter.
##
## Problem this fixes: the legacy `Y_df_for_mcsd.rds` carried 278 of the 313 panel
## genes, but the selection rule was not reconstructable from any script in the tree
## (it was not a clean detection threshold). This script rebuilds the bundle from the
## SAME cells and labels with a documented filter:
##
##   KEEP gene g  iff  max_c  mean( count_{ig} > 0  |  cellType(i) = c )  >= DET_MIN
##
## i.e. max-per-celltype detection >= 5% (DET_MIN, env R_DET_MIN, default 0.05),
## matching the Mel canonical builder's gene filter. Computed over the 9 scClassify
## cell types actually fit.
##
## Source of counts: the legacy bundle's `df` already holds ALL 313 panel gene
## columns as integer counts (verified identical to its `Y` for shared genes), so we
## rebuild from it directly — same 126,432 cells, same scClassify `cellType` labels.
## Output: data/breast_cancer/Y_df_for_mcsd_det05.rds  (drop-in for R_YDF).

setwd("/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv")
suppressPackageStartupMessages(library(Matrix))
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), sprintf(...)))

DET_MIN <- as.numeric(Sys.getenv("R_DET_MIN", unset = "0.05"))
IN  <- "data/breast_cancer/Y_df_for_mcsd.rds"
RAW <- "data/breast_cancer/Y_df_janesick_raw.rds"
OUT <- Sys.getenv("R_OUT", unset = "data/breast_cancer/Y_df_for_mcsd_det05.rds")

say("Loading legacy bundle %s", IN)
mc  <- readRDS(IN)
df  <- mc$df
panel <- colnames(readRDS(RAW)$Y)                      # 313-gene full panel
present <- intersect(panel, colnames(df))
if (length(present) != length(panel))
  stop(sprintf("only %d/%d panel genes present in df", length(present), length(panel)))
G <- as.matrix(df[, present, drop = FALSE]); storage.mode(G) <- "double"
stopifnot(all(G == round(G)))                          # integer counts
say("Full panel matrix: %d cells x %d genes", nrow(G), ncol(G))

## ---- max-per-celltype detection filter ----
ct <- as.character(df$cellType)
lev <- sort(unique(ct))
say("cell types (%d): %s", length(lev), paste(lev, collapse = ", "))
det_by_ct <- vapply(lev, function(c) colMeans(G[ct == c, , drop = FALSE] > 0),
                    numeric(ncol(G)))                  # genes x celltypes
max_det <- apply(det_by_ct, 1, max)
keep <- max_det >= DET_MIN
say("DET_MIN = %.3f (max-per-celltype); KEEP %d / %d genes (add back %d, drop %d)",
    DET_MIN, sum(keep), length(keep),
    sum(keep) - length(intersect(colnames(mc$Y), present)),
    sum(!keep))
added   <- setdiff(present[keep], colnames(mc$Y))
removed <- setdiff(colnames(mc$Y), present[keep])
say("  added back vs legacy 278: %s", if (length(added)) paste(added, collapse=", ") else "none")
say("  newly removed vs legacy 278: %s", if (length(removed)) paste(removed, collapse=", ") else "none")

Y <- G[, keep, drop = FALSE]; storage.mode(Y) <- "integer"

## ---- write bundle (same df; Y = filtered panel). nCount recomputed by the fit builder. ----
out <- list(Y = Y, df = df,
            gene_filter = list(rule = "max-per-celltype detection >= DET_MIN",
                               DET_MIN = DET_MIN, n_panel = length(panel),
                               n_kept = sum(keep),
                               max_det = max_det, kept = colnames(Y),
                               added_vs_legacy = added, removed_vs_legacy = removed,
                               built = format(Sys.time())))
saveRDS(out, OUT)
say("Saved %s  (%d cells x %d genes)", OUT, nrow(Y), ncol(Y))

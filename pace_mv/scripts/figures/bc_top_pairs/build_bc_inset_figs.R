## Reproduce legacy BC inset figures using the EXISTING helper compose_inset_figure().
## Two figures, same two tumour-ring ROIs:
##   1) Macrophage <- Tumour, genes MRC1 + APOC1
##   2) Stromal    <- Tumour, genes ADH1B + CCDC80
suppressMessages({
  library(ggplot2); library(patchwork); library(cowplot); library(ggrastr)
})

base <- "/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/pace_mv"
source(file.path(base, "scripts/helpers/compose_inset_figure.R"))
fig_dir <- file.path(base, "plots/bc_top_pairs")

## ---- build df ----
ymv <- readRDS(file.path(base, "data/breast_cancer/Y_df_for_mcsd.rds"))
Y <- ymv$Y; meta <- ymv$df
df <- data.frame(
  x = meta$x, y = meta$y,
  celltype = as.character(meta$cellType),
  stringsAsFactors = FALSE
)
genes_all <- c("ADH1B", "CCDC80", "MRC1", "APOC1", "MYLK", "KRT14")
stopifnot(all(genes_all %in% colnames(Y)))
for (g in genes_all) df[[g]] <- log1p(1e4 * Y[, g] / meta$nCount)

## ---- VALIDATION ----
cat("=== VALIDATION ===\n")
cat("nrow(df):", nrow(df), "\n")
cat("genes present in Y:", paste(genes_all, genes_all %in% colnames(Y), sep="="), "\n")
cat("celltype table:\n"); print(table(df$celltype))
expected <- c("Tumour","Myoepithelial","Stromal","Macrophage","T_Cell",
              "B_Cell","Endothelial","Dendritic_Cell","Mast")
cat("all 9 expected types present:", all(expected %in% df$celltype), "\n")
cat("expression ranges (log CP10K):\n")
for (g in genes_all) cat(sprintf("  %-7s [%.2f, %.2f]\n", g, min(df[[g]]), max(df[[g]])))

## ---- two tumour-ring ROIs (user-specified, right-side rings) ----
shrink_box <- function(xlim, ylim, f) {
  cx <- mean(xlim); cy <- mean(ylim)
  hx <- diff(xlim) / 2 * f
  hy <- diff(ylim) / 2 * f
  list(xlim = c(cx - hx, cx + hx), ylim = c(cy - hy, cy + hy))
}
shrink_to_nonoverlap <- function(xlim1, ylim1, xlim2, ylim2, gap = 0.02) {
  c1  <- c(mean(xlim1), mean(ylim1)); c2 <- c(mean(xlim2), mean(ylim2))
  hx1 <- diff(xlim1) / 2; hy1 <- diff(ylim1) / 2
  hx2 <- diff(xlim2) / 2; hy2 <- diff(ylim2) / 2
  dx  <- abs(c1[1] - c2[1]); dy <- abs(c1[2] - c2[2])
  fx <- dx / (hx1 + hx2); fy <- dy / (hy1 + hy2)
  f  <- min(1, max(fx, fy) - gap); f <- max(f, 0.05)
  list(factor = f, orange = shrink_box(xlim1, ylim1, f), blue = shrink_box(xlim2, ylim2, f))
}
xlim1 <- c(6300, 7300); ylim1 <- c(3700, 4700)
xlim2 <- c(6500, 7500); ylim2 <- c(2300, 3300)
res <- shrink_to_nonoverlap(xlim1, ylim1, xlim2, ylim2, gap = 0.01)
cat("\nshrink factor:", res$factor, "\n")
xlim1 <- res$orange$xlim; ylim1 <- res$orange$ylim   # box1 -> orange
xlim2 <- res$blue$xlim;   ylim2 <- res$blue$ylim     # box2 -> blue

roi_count <- function(ct, xl, yl) {
  s <- df[df$celltype == ct & df$x>=xl[1] & df$x<=xl[2] & df$y>=yl[1] & df$y<=yl[2], ]
  nrow(s)
}
cat("\n=== ROI focal-cell counts ===\n")
cat(sprintf("ROI1 x[%g,%g] y[%g,%g]: Stromal=%d Macrophage=%d Tumour=%d\n",
            xlim1[1],xlim1[2],ylim1[1],ylim1[2],
            roi_count("Stromal",xlim1,ylim1), roi_count("Macrophage",xlim1,ylim1),
            roi_count("Tumour",xlim1,ylim1)))
cat(sprintf("ROI2 x[%g,%g] y[%g,%g]: Stromal=%d Macrophage=%d Tumour=%d\n",
            xlim2[1],xlim2[2],ylim2[1],ylim2[2],
            roi_count("Stromal",xlim2,ylim2), roi_count("Macrophage",xlim2,ylim2),
            roi_count("Tumour",xlim2,ylim2)))

## ---- figure 1: Macrophage <- Tumour, MRC1 + APOC1 (reproduces reference) ----
f1 <- file.path(fig_dir, "BC_inset_Macrophage_MRC1_APOC1.pdf")
cairo_pdf(f1, width = 13.5, height = 7.5)
compose_inset_figure(df, genes = c("MRC1","APOC1"),
                     focal = "Macrophage", neighbour = "Tumour",
                     xlim1, ylim1, xlim2, ylim2,
                     roi_col1 = "darkorange", roi_col2 = "steelblue")
dev.off()
cat("\nSAVED:", f1, "\n")

## ---- figure 2: Stromal <- Tumour, ADH1B + CCDC80 ----
f2 <- file.path(fig_dir, "BC_inset_Stromal_ADH1B_CCDC80.pdf")
cairo_pdf(f2, width = 13.5, height = 7.5)
compose_inset_figure(df, genes = c("ADH1B","CCDC80"),
                     focal = "Stromal", neighbour = "Tumour",
                     xlim1, ylim1, xlim2, ylim2,
                     roi_col1 = "darkorange", roi_col2 = "steelblue")
dev.off()
cat("SAVED:", f2, "\n")

## ---- figure 3: Myoepithelial <- Tumour, MYLK + KRT14 (top-2 MCSD) ----
f3 <- file.path(fig_dir, "BC_inset_Myoepithelial_MYLK_KRT14.pdf")
cairo_pdf(f3, width = 13.5, height = 7.5)
compose_inset_figure(df, genes = c("MYLK","KRT14"),
                     focal = "Myoepithelial", neighbour = "Tumour",
                     xlim1, ylim1, xlim2, ylim2,
                     roi_col1 = "darkorange", roi_col2 = "steelblue")
dev.off()
cat("SAVED:", f3, "\n")

# =============================================================================
# build_wta_spe.R
# -----------------------------------------------------------------------------
# Build a SpatialExperiment for the 10x Atera WTA FFPE human breast cancer block
# (whole-transcriptome, 18,028 genes), mirroring the canonical BC SPE structure
# so it is a drop-in input for the PACE-MV BC pipeline.
#
# Inputs (10x Xenium Onboard Analysis v4 outs + Xenium Explorer exports):
#   cell_feature_matrix.h5                 counts (27,104 features; 18,028 Gene Expression)
#   cells.csv.gz                           per-cell centroids (um) + segmentation metrics
#   WTA_..._cell_groups.csv                10x cell-type annotation (cell_id, group, color)
#   WTA_..._gene_groups.csv                10x curated marker sets (stored in metadata)
#
# Cell-type mapping (user-approved):
#   - all DCIS subtypes + Apocrine            -> Tumour      (collapse)
#   - CAFs + CXCL14+ Fibroblasts + Pericytes  -> Stromal     (Pericyte folded in)
#   - Plasma & Mast Cell Mixture              -> DROPPED      (ambiguous)
#   -> final 9-type scheme matches the canonical BC SPE.
#
# Output: data/breast_cancer/spe_wta_ffpe_breast_cancer.rds
# NOTE: this only builds the SPE. It does NOT run any PACE fit.
# =============================================================================

suppressPackageStartupMessages({
  library(DropletUtils)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(dbscan)
})

data_dir   <- "/Users/elijahwillie/Documents/Academic/PhD/lateral_spillover/data/breast_cancer"
outs_dir   <- file.path(data_dir, "WTA_Preview_FFPE_Breast_Cancer_outs")
h5_path    <- file.path(outs_dir, "cell_feature_matrix.h5")
cells_path <- file.path(outs_dir, "cells.csv.gz")
grp_path   <- file.path(data_dir, "WTA_Preview_FFPE_Breast_Cancer_cell_groups.csv")
gg_path    <- file.path(data_dir, "WTA_Preview_FFPE_Breast_Cancer_gene_groups.csv")
out_path   <- file.path(data_dir, "spe_wta_ffpe_breast_cancer.rds")

say <- function(...) cat(sprintf(...), "\n")

# ---- 1. Counts: read full feature-barcode matrix, keep Gene Expression only ----
say("[1/6] reading cell_feature_matrix.h5 ...")
sce_all <- DropletUtils::read10xCounts(h5_path, col.names = TRUE)
say("    raw matrix: %d features x %d cells", nrow(sce_all), ncol(sce_all))

is_gene <- rowData(sce_all)$Type == "Gene Expression"
say("    Gene Expression features: %d (dropping %d control features)",
    sum(is_gene), sum(!is_gene))
sce_all <- sce_all[is_gene, ]

# Use gene symbols as rownames (canonical BC SPE convention), unique-ified.
gene_symbol <- rowData(sce_all)$Symbol
rownames(sce_all) <- make.unique(as.character(gene_symbol))
counts_mat <- as(assay(sce_all, "counts"), "CsparseMatrix")
row_data   <- rowData(sce_all)[, c("ID", "Symbol", "Type")]
barcodes   <- colnames(sce_all)
rm(sce_all); invisible(gc())

# ---- 2. Per-cell coordinates + segmentation metrics ----
say("[2/6] reading cells.csv.gz ...")
cells <- read.csv(gzfile(cells_path), stringsAsFactors = FALSE)
rownames(cells) <- cells$cell_id

# ---- 3. 10x cell-type annotation; drop the ambiguous mixture ----
say("[3/6] reading cell-group annotation ...")
grp <- read.csv(grp_path, stringsAsFactors = FALSE)   # cell_id, group, color
n_mixture <- sum(grp$group == "Plasma & Mast Cell Mixture")
grp <- grp[grp$group != "Plasma & Mast Cell Mixture", ]
say("    dropped ambiguous 'Plasma & Mast Cell Mixture': %d cells", n_mixture)

# Map the 10x groups to the canonical 9-type PACE scheme.
celltype_map <- c(
  "11q13 High Grade DCIS Cells"              = "Tumour",
  "11q13 High Grade DCIS Cells (Mitotic)"    = "Tumour",
  "11q13 High Grade DCIS Tumor Cells (G1/S)" = "Tumour",
  "Luminal-like Amorphous DCIS Cells"        = "Tumour",
  "Basal-like Structured DCIS Cells"         = "Tumour",
  "Apocrine Cells"                           = "Tumour",
  "CAFs, Low Grade DCIS Associated"          = "Stromal",
  "CAFs, High Grade DCIS Associated"         = "Stromal",
  "CXCL14+ Fibroblasts"                      = "Stromal",
  "Pericytes"                                = "Stromal",
  "Endothelial Cells"                        = "Endothelial",
  "Myoepithelial Cells"                      = "Myoepithelial",
  "T Lymphocytes"                            = "T_Cell",
  "Macrophages"                              = "Macrophage",
  "Myeloid Cells"                            = "Macrophage",
  "Dendritic Cells"                          = "Dendritic_Cell",
  "Plasma Cells"                             = "B_Cell",
  "Mast Cells"                               = "Mast"
)
stopifnot(all(grp$group %in% names(celltype_map)))   # fail loud on any unmapped group
grp$celltype <- celltype_map[grp$group]
rownames(grp) <- grp$cell_id

# ---- 4. Align the three sources on the SAME cells (matrix ∩ coords ∩ annotation) ----
say("[4/6] aligning matrix / coords / annotation ...")
keep <- Reduce(intersect, list(barcodes, cells$cell_id, grp$cell_id))
say("    cells kept (annotated, non-ambiguous, in matrix+coords): %d", length(keep))

counts_mat <- counts_mat[, keep, drop = FALSE]
cells_k    <- cells[keep, ]
grp_k      <- grp[keep, ]

# ---- 5. Validate coordinate units via within-section nearest-neighbour distance ----
say("[5/6] validating coordinate units (within-section 1-NN) ...")
coords <- as.matrix(cells_k[, c("x_centroid", "y_centroid")])
nn1    <- dbscan::kNN(coords, k = 1)$dist[, 1]
say("    median 1-NN = %.2f  (expect ~8-15 um for single-cell Xenium); IQR [%.2f, %.2f]",
    median(nn1), quantile(nn1, 0.25), quantile(nn1, 0.75))
if (median(nn1) < 4 || median(nn1) > 25)
  warning("median 1-NN outside the expected 8-15 um single-cell range; check coord units")

# ---- 6. Assemble the SpatialExperiment ----
say("[6/6] assembling SpatialExperiment ...")
nCount <- Matrix::colSums(counts_mat)            # library size over Gene Expression genes

col_data <- DataFrame(
  cell_id           = keep,
  celltype          = factor(grp_k$celltype),    # final 9-type PACE scheme
  annotation_10x    = grp_k$group,               # original 10x cell-group label
  annotation_color  = grp_k$color,
  slide             = "WTA_breast",              # single physical section (mirrors BC 'slide')
  sample_id         = "wta_ffpe_breast",
  nCount            = nCount,
  transcript_counts = cells_k$transcript_counts,
  total_counts      = cells_k$total_counts,
  cell_area         = cells_k$cell_area,
  nucleus_area      = cells_k$nucleus_area,
  nucleus_count     = cells_k$nucleus_count,
  segmentation      = cells_k$segmentation_method,
  row.names         = keep
)

spe <- SpatialExperiment(
  assays        = list(counts = counts_mat),
  colData       = col_data,
  rowData       = row_data,
  spatialCoords = as.matrix(setNames(
                    cells_k[, c("x_centroid", "y_centroid")],
                    c("cell_centroid_x", "cell_centroid_y")))
)
# curated 10x marker sets, kept alongside for downstream validation (not used in the fit)
metadata(spe)$gene_groups <- read.csv(gg_path, stringsAsFactors = FALSE)
metadata(spe)$source      <- "10x Atera WTA FFPE Human Breast Cancer (DCIS Grade 3); single section"
metadata(spe)$pixel_size  <- 0.2125

say("    final SPE: %d genes x %d cells", nrow(spe), ncol(spe))
say("    cell-type table:")
print(table(spe$celltype))

saveRDS(spe, out_path)
say("WROTE %s", out_path)

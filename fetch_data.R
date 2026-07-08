## fetch_data.R -- download the PACE manuscript analysis data from Zenodo and
## place it at the paths the notebooks expect. Run from the repository root:
##
##   Rscript fetch_data.R
##
## The data are not stored in git (raw inputs and fits total several GB). They
## live in a single Zenodo archive that unzips to `data/` and `pace_mv/data/`
## under the repository root, reproducing the layout the notebooks read from.
##
## >>> BEFORE FIRST USE: create the Zenodo deposit, upload the archive, and set
## >>> ZENODO_RECORD (and the archive name + md5) below. <<<

## ---- Zenodo record -- FILL THESE IN AFTER DEPOSITING -----------------------
ZENODO_RECORD <- "REPLACE_WITH_ZENODO_RECORD_ID"   # e.g. "14812345"
ARCHIVE_NAME  <- "pace_manuscript_data.zip"
ARCHIVE_MD5   <- "REPLACE_WITH_ARCHIVE_MD5"        # md5sum of the .zip

## Files the notebooks require, relative to the repository root, so we can verify
## the archive unpacked correctly. (Canonical fits are refit by the notebooks and
## are optional; the raw inputs below are required.)
REQUIRED <- c(
  "data/spe_10x_nuclei_withMetrics.rds",                       # breast cancer SPE
  "data/simvi_melanoma/14708000/Melanoma_5612.h5ad",           # melanoma h5ad
  "pace_mv/data/breast_cancer/Y_df_for_mcsd.rds",
  "pace_mv/data/breast_cancer/mvpql_percell_hc.rds",           # PACE fit for method comparison
  "pace_mv/data/simvi_melanoma/Y_df_for_mcsd.rds"
)

## ---- download + verify + unpack --------------------------------------------
repo_root <- normalizePath(".")
if (ZENODO_RECORD == "REPLACE_WITH_ZENODO_RECORD_ID")
  stop("Set ZENODO_RECORD (and ARCHIVE_MD5) in fetch_data.R first.", call. = FALSE)

url <- sprintf("https://zenodo.org/records/%s/files/%s?download=1",
               ZENODO_RECORD, ARCHIVE_NAME)
dest <- file.path(tempdir(), ARCHIVE_NAME)

message("Downloading ", ARCHIVE_NAME, " from Zenodo record ", ZENODO_RECORD, " ...")
utils::download.file(url, dest, mode = "wb")

if (ARCHIVE_MD5 != "REPLACE_WITH_ARCHIVE_MD5") {
  got <- tools::md5sum(dest)[[1]]
  if (!identical(unname(got), ARCHIVE_MD5))
    stop("md5 mismatch: expected ", ARCHIVE_MD5, " got ", got, call. = FALSE)
  message("md5 verified.")
}

message("Unpacking into ", repo_root, " ...")
utils::unzip(dest, exdir = repo_root)

## ---- sanity check -----------------------------------------------------------
missing <- REQUIRED[!file.exists(file.path(repo_root, REQUIRED))]
if (length(missing)) {
  warning("These expected files are still missing after unpacking:\n  ",
          paste(missing, collapse = "\n  "))
} else {
  cat("\nData ready. All required inputs are in place.\n")
}

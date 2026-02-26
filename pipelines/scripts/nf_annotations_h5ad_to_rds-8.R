library(anndata)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No script directory provided.")
}

script_dir <- args[1]
source(file.path(script_dir, "pipeline_utils.R"))

SEURATCOUNTFILE <- get_param_value("--seuratcountfile")
ANNDATAFILE     <- get_param_value("--annfile")
INTEGRATED      <- get_param_value("--integrate")
SAMPLEID        <- get_param_value("--sample_id")
ANNOFIELD1      <- get_param_value("--annofield1")

DATA_INTEGRATED <- toupper(INTEGRATED) != "NONE"
OUTFILE         <- sub("\\.h5ad$", ".rds", ANNDATAFILE)

cat("Seurat file:", SEURATCOUNTFILE, "\n")
cat("AnnData file:", ANNDATAFILE, "\n")
cat("Integrated:", DATA_INTEGRATED, "\n")
cat("Sample ID field:", SAMPLEID, "\n")
cat("Output:", OUTFILE, "\n")

# load adata
adata <- read_h5ad(ANNDATAFILE, backed = "r")
obs   <- as.data.frame(adata$obs)
rm(adata)
gc()

# load seurat object
if (!DATA_INTEGRATED) {
  sample <- as.character(obs[[SAMPLEID]][1])
  count_data_list <- readRDS(SEURATCOUNTFILE)
  if (!is.list(count_data_list)) {
    stop("Expected a list object for count_data_list.")
  }
  count_data <- count_data_list[[sample]]
} else {
  count_data <- readRDS(SEURATCOUNTFILE)
  sample <- "INTEGRATED"
}

# check metadata alignment
if (!identical(rownames(obs), rownames(count_data@meta.data))) {
  stop(sprintf(
    "Cell names are not identical or not in the same order for sample '%s'.",
    sample
  ))
}
cat(sprintf(
  "Cell names match for sample '%s'. Transferring metadata...\n",
  sample
))

# Always transfer primary annotation field
if (!ANNOFIELD1 %in% colnames(obs)) {
  stop(sprintf("Primary annotation field '%s' not found in AnnData obs.", ANNOFIELD1))
}
count_data@meta.data[[ANNOFIELD1]] <- obs[[ANNOFIELD1]]
# Transfer all scanvipred_ and celltypist_ columns automatically
cols_to_transfer <- colnames(obs)[
  startsWith(colnames(obs), "scanvipred_") |
    startsWith(colnames(obs), "celltypist_")
]
if (length(cols_to_transfer) > 0) {
  count_data@meta.data[, cols_to_transfer] <- obs[, cols_to_transfer, drop = FALSE]
}
cat("Transferred", length(cols_to_transfer), "dynamic annotation columns.\n")
saveRDS(count_data, OUTFILE)
cat("Done.\n")

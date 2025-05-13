library(reticulate)
use_condaenv("scseq", require = T)
sc <- reticulate::import("scanpy")
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

# Check if the script directory argument is provided
if (length(args) == 0) {
    stop("No script directory provided.")
}

# Store the script directory
script_dir <- args[1]

# Print the script directory for debugging
cat("Script directory:", script_dir, "\n")

# Load the helper function
source(file.path(script_dir, "pipeline_utils.R"))


SEURATCOUNTFILE <- get_param_value("--seuratcountfile")
ANNDATAFILE <- get_param_value("--annfile")
INTEGRATED <- get_param_value("--integrate")
SAMPLEID <- get_param_value("--sample_id")
ANNOFIELD1 <- get_param_value("--annofield1")

DATA_INTEGRATED <- ifelse(toupper(INTEGRATED) == "NONE", FALSE, TRUE)
OUTFILE <- gsub(".h5ad", ".rds", ANNDATAFILE)

cat(SEURATCOUNTFILE, "\n")
cat(ANNDATAFILE, "\n")
cat(DATA_INTEGRATED, "\n")
cat(SAMPLEID, "\n")
cat(OUTFILE, "\n")

# backprop_metadata <- function(adata, count_data) {
#   assertion_result <- identical(rownames(adata$obs), rownames(count_data@meta.data))
#   if (assertion_result) {
#     cat("The row names or cells are identical for files \n")
#     # merge metadata from h5ad to seurat object
#     count_data@meta.data <- merge(count_data@meta.data, adata$obs, by = "row.names", all.x = TRUE)
#     saveRDS(count_data, OUTFILE)
#   } else {
#     stop(paste0("Assertion failed: The row names (cells) are not identical for ", sample,"; Try to identify what is missing/out of order"))
#   }
# }


if (!DATA_INTEGRATED) {
  adata <- py_to_r(sc$read_h5ad(ANNDATAFILE))
  sample <- as.character(adata$obs[[SAMPLEID]][[1]])
  count_data_list <- readRDS(SEURATCOUNTFILE)
  # Check if count_data_list is a list
  if (is.list(count_data_list)) {
    count_data <- count_data_list[[sample]]
  } else {
    stop("Expected a list object for count_data_list, but got a different type.")
  }
} else {
  count_data <- readRDS(SEURATCOUNTFILE) # read raw count data file
  adata <- py_to_r(sc$read_h5ad(ANNDATAFILE)) # read annotated data file
}

# # backpropagate metadata
# backprop_metadata(adata, count_data)

assertion_result <- identical(rownames(adata$obs), rownames(count_data@meta.data))

if (assertion_result) {
  cat("The row names or cells are identical for files \n")
  # merge metadata from h5ad to seurat object
  count_data@meta.data[[ANNOFIELD1]] <- adata$obs[[ANNOFIELD1]]
  count_data@meta.data$`scanvipred_mouse-tabula-muris-senis-droplet` <- adata$obs$`scanvipred_mouse-tabula-muris-senis-droplet`
  count_data@meta.data$`celltypist_mouse-tabula-muris-senis-droplet.predicted_labels` <- adata$obs$`celltypist_mouse-tabula-muris-senis-droplet.predicted_labels`
  count_data@meta.data$`celltypist_mouse-tabula-muris-senis-droplet.over_clustering` <- adata$obs$`celltypist_mouse-tabula-muris-senis-droplet.over_clustering`
  count_data@meta.data$`celltypist_mouse-tabula-muris-senis-droplet.conf_score` <- adata$obs$`celltypist_mouse-tabula-muris-senis-droplet.conf_score`
  count_data@meta.data$`celltypist_mouse-tabula-muris-senis-droplet.majority_voting` <- adata$obs$`celltypist_mouse-tabula-muris-senis-droplet.majority_voting`
  saveRDS(count_data, OUTFILE)
} else {
  stop(paste0("Assertion failed: The row names (cells) are not identical for ", sample,"; Try to identify what is missing/out of order"))
}

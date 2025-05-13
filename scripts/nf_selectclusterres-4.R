library(Seurat)
library(dplyr)
options(bitmapType='cairo')
# Parse command-line arguments + import utility functions
# Get the script directory from the command line arguments
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

DATASET_LABEL <- get_param_value("--datlabel")
INFILE <- get_param_value("--infile")
OUTFILE <- get_param_value("--outfile")
INTEGRATION <- get_param_value("--integrate")

CLUST_RES <- as.numeric(get_param_value("--clusterreslist")) 
RES <- get_param_value("--clusterres")
SUBCLUSTER = get_param_value("--subcluster")
SUBCLUSTER = ifelse(SUBCLUSTER == 1, T, F)

# read data
seurat_data <- readRDS(INFILE)

# Finalize cluster resolution for downstream analysis ---------------------------
INTEGRATED = T
if(INTEGRATION == "NONE") {
  INTEGRATED = F
}
if (is.list(seurat_data)) {
  assaylab <- ifelse(INTEGRATED & !SUBCLUSTER,"integrated", DefaultAssay(seurat_data[[1]])) # assume list of objects underwent identical preprocessing
} else {
  assaylab <- ifelse(INTEGRATED & !SUBCLUSTER,"integrated", DefaultAssay(seurat_data))
}

# RNA_snn_res.0.5 <assay>_snn_res.<res>
clusterlab_list <- setNames(
  lapply(CLUST_RES, function(name) paste0(assaylab, "_snn_res.", name)),
  as.character(CLUST_RES)
)
print(clusterlab_list)
print("\n")
print(RES)
if (is.list(seurat_data)) {
  seurat_data <- lapply(seurat_data, function(seurat_obj) {
    Idents(seurat_obj) <- clusterlab_list[[RES]]
    seurat_obj$seurat_clusters <- seurat_obj[[clusterlab_list[[RES]]]] # extra setting - as I will typically refer to seurat cluster when using python to analyze
    return(seurat_obj)
  })
} else {
  Idents(seurat_data) <- clusterlab_list[[RES]]
  seurat_data$seurat_clusters <- seurat_data[[clusterlab_list[[RES]]]] # extra setting - as I will typically refer to seurat cluster when using python to analyze
}

saveRDS(seurat_data, file = OUTFILE)



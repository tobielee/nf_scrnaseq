#!/usr/local/bin R

# library(scran)
library(ggplot2)
library(Seurat)
library(dplyr)
options(bitmapType='cairo')
options(future.globals.maxSize = 16000 * 1024^2) # 16 GB memory limit
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
# Retrieve parameter values
DATASET_LABEL <- get_param_value("--datlabel")
INFILE <- get_param_value("--infile")
OUTFILE <- get_param_value("--outfile")
NORM <- get_param_value("--norm")
INTEGRATION <- get_param_value("--integrate")

DOUBLETDETECT = get_param_value("--detectdoub")
REMOVEDOUBLETS = get_param_value("--removedoub")
DOUBLETDETECT = ifelse(DOUBLETDETECT == 1, T, F)
REMOVEDOUBLETS = ifelse(REMOVEDOUBLETS == 1, T, F)

cat(NORM, INTEGRATION, REMOVEDOUBLETS)

seurat_data <- readRDS(INFILE)

# Normalization/integration
# scran normalization -----------------------------------------------------
# TODO assume log_transformation
ScranNorm <- function(seurat_obj, already_clustered = DOUBLETDETECT, log_trans = T) {
  sce <- as.SingleCellExperiment(seurat_obj)
  if (DOUBLETDETECT) {
    clusters <- seurat_obj@meta.data$seurat_clusters
  } else {
    clusters <- quickCluster(sce,
                             use.ranks = FALSE, # suggested by the authors
                             min.size = 100) # require at least 100 cells per cluster
  }
  print(table(clusters))
  sce <- computeSumFactors(sce,
                           clusters = clusters,
                           min.mean = 0.1)
  # size_factors = sizeFactors(
  #   computeSumFactors(
  #     SingleCellExperiment(
  #       list(counts=data_mat)),
  #     clusters = clusters,
  #     min.mean = 0.1,
  #     BPPARAM = MulticoreParam()
  #   )
  # )
  # apply size factors to generate log normalized data
  # print(summary(sizeFactors(sce)))
  # sce <- logNormCounts(sce)
  seurat_obj$size_factors <- sizeFactors(sce)
  original_counts <- seurat_obj@assays$RNA$counts # TODO this is not ideal way of fetching counts 
  seurat_obj@assays$RNA$counts <- original_counts/seurat_obj$size_factors
  # Check if log transformation is required
  if (log_trans) {
    cat("Note! Performing log1p-transformation after normalization.\n")
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj@assays$RNA$counts <- original_counts
  } else {
    seurat_obj@assays$RNA$data <- seurat_obj@assays$RNA$counts
    seurat_obj@assays$RNA$counts <- original_counts
    cat("No log-transformation performed after normalization.\n")
  }
  return(seurat_obj)
}
isListOfLists <- any(sapply(seurat_data, function(x) length(x) > 1))

if (isListOfLists) {
  for (sample_name in names(seurat_data)) {
    if (toupper(NORM) == "SCRAN") {
      seurat_data[[sample_name]] <- sapply(seurat_data[[sample_name]], ScranNorm)
    }
    else if (toupper(NORM) == "SCT") {
      seurat_data[[sample_name]] <- sapply(seurat_data[[sample_name]], SCTransform)
      # seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE) # TODO add other options for normalization SCT transform etc.
    }
    else {
      seurat_data[[sample_name]] <- sapply(seurat_data[[sample_name]], NormalizeData)
      # seurat_obj <- NormalizeData(object = seurat_obj, verbose = FALSE)
    }
    if (toupper(NORM) != "SCT") {
      seurat_data[[sample_name]] <- sapply(seurat_data[[sample_name]], FindVariableFeatures) # TODO adjust variable feature amount
      # seurat_obj <- FindVariableFeatures(object = seurat_obj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
    }
  }
} else {
  if (toupper(NORM) == "SCRAN") {
    seurat_data<- lapply(seurat_data, ScranNorm)
  } else if (toupper(NORM) == "SCT") {
    seurat_data <- lapply(seurat_data, SCTransform)
    # seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE) # TODO add other options for normalization SCT transform etc.
  } else {
    seurat_data <- lapply(seurat_data, NormalizeData)
    # seurat_obj <- NormalizeData(object = seurat_obj, verbose = FALSE)
  }
  if (toupper(NORM) != "SCT") {
    seurat_data <- lapply(seurat_data, FindVariableFeatures) # TODO adjust variable feature amount
    # seurat_obj <- FindVariableFeatures(object = seurat_obj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  }
}
# per sample_name integration or leave separate -------------------------------------------
scale_pca_integrate <- function(seurat_input) {
  if (is.list(seurat_input)) {
    # If the input is a list of Seurat objects
    features <- SelectIntegrationFeatures(object.list = seurat_input)
    seurat_output <- lapply(seurat_input, function(obj) {
      obj <- ScaleData(obj, features = features, verbose = FALSE) # TODO may need to avoid this if SCT
      obj <- RunPCA(obj, features = features, verbose = FALSE)
      return(obj)
    })
    cat("Integrating...")
    start <- Sys.time()
    anchors <- FindIntegrationAnchors(object.list = seurat_output, anchor.features = features, reduction = "rpca") # TODO adjust integration method
    seurat_output <- IntegrateData(anchorset = anchors)
    seurat_output[["RNA"]] <- JoinLayers(seurat_output[["RNA"]]) # immediately join layers for seurat v5
    print(Sys.time() - start )
  } else {
    # If the input is a single Seurat object
    # TODO This might not be needed - as it is done in next step for clustering 
    features <- VariableFeatures(object = seurat_input)
    seurat_output <- ScaleData(seurat_input, features = features, verbose = FALSE)
    seurat_output <- RunPCA(seurat_output, features = features, verbose = FALSE)
  }
  return(seurat_output)
}

if (toupper(INTEGRATION) == "GROUP") { # TODO need to fix this this still seems to list everything out though grouped doesn't reduce to each sample_name
  seurat_data <- lapply(seurat_data, function(sample_name) {
    lapply(sample_name, scale_pca_integrate)
  })
} else {
  if (isListOfLists) {
    seurat_data <- unlist(seurat_data, recursive = FALSE)
  }
  if (toupper(INTEGRATION) == "ALL") {
    seurat_data = scale_pca_integrate(seurat_data)
  } else {
    seurat_data <- lapply(seurat_data, scale_pca_integrate)
  }
}


saveRDS(seurat_data, file = OUTFILE)
# saveRDS(seurat_data, file = paste0("/Users/tlee/Documents/scrnaseq/outs/", DATASET_LABEL,filesuffix))


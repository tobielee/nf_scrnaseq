library(ggplot2)
library(Seurat)
library(tidyverse)
library(clustree)
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

NPCA_COMPS <- as.numeric(get_param_value("--npca"))
CLUST_RES <- as.numeric(get_param_value("--clusterreslist")) 
TOPN_VARFEAT <- as.numeric(get_param_value("--topnfeat"))
REDUCTION = get_param_value("--reduction")

SUBCLUSTER = get_param_value("--subcluster")
SUBCLUSTER = ifelse(SUBCLUSTER == 1, T, F)
# SUBCLUSTERNAME = "epi_stromal"
# subcluster_integration_modes <- list("HARMONY" = "harmony", "RPCA" = "integrated.rpca", "MNN" = "integrated.mnn")
# if (toupper(INTEGRATION) %in% names(subcluster_integration_modes)) {
#   REDUCTION <- subcluster_integration_modes[[toupper(INTEGRATION)]]
# }

PLOTSDIR <- file.path("cluster_deg") # TODO could expose this as a parameter
if (!dir.exists(PLOTSDIR)) {
  dir.create(PLOTSDIR)
}
if (SUBCLUSTER) {
  seurat_data <- readRDS(INFILE)
  seurat_data[["RNA"]] <- JoinLayers(seurat_data[["RNA"]])
} else {
  seurat_data <- readRDS(INFILE)
}

cat(CLUST_RES, "##!@#@!#@#") # TODO remove
# post integration clustering --------------------------------------------------------------
umap_cluster <- function(seurat_obj) {
  if (INTEGRATION != "NONE") {
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE) # need to scale again post integration 
  }
  seurat_obj <- RunPCA(seurat_obj, npcs = NPCA_COMPS, verbose = FALSE)
  png(paste0(DATASET_LABEL,"_", "elbowplot_check.png"), width = 750, height = 500)
  print(ElbowPlot(seurat_obj, ndims = NPCA_COMPS))
  dev.off()
  seurat_obj <- RunUMAP(seurat_obj, reduction = REDUCTION, dims = 1:NPCA_COMPS)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = REDUCTION, dims = 1:NPCA_COMPS)
  seurat_obj <- FindClusters(seurat_obj, resolution = CLUST_RES)
  return(seurat_obj)
}

INTEGRATED = T
if(INTEGRATION == "NONE") {
  INTEGRATED = F
}
if (is.list(seurat_data)) {
  assaylab <- ifelse(INTEGRATED & !SUBCLUSTER,"integrated", DefaultAssay(seurat_data[[1]])) # assume list of objects underwent identical preprocessing
  seurat_data <- lapply(seurat_data, umap_cluster)
} else {
  assaylab <- ifelse(INTEGRATED & !SUBCLUSTER,"integrated", DefaultAssay(seurat_data))
  seurat_data <- umap_cluster(seurat_data)
}


# RNA_snn_res.0.5 <assay>_snn_res.<res>
clusterlab_list <- setNames(
  lapply(CLUST_RES, function(name) paste0(assaylab, "_snn_res.", name)),
  as.character(CLUST_RES)
)

plot_cluster_res <- function(seurat_obj) {
  if (SUBCLUSTER) {
    sample_name <- SUBCLUSTERNAME
  } else {
    sample_name <- ifelse(INTEGRATION == "ALL", DATASET_LABEL, seurat_obj$sample[[1]]) 
  }
  sample_subdir <- file.path(PLOTSDIR, sample_name)
  if (!dir.exists(sample_subdir)) {
    dir.create(sample_subdir)
  }
  clusterprefix <- paste0(DefaultAssay(seurat_obj), "_snn_res.")
  # Clustree
  clustree_filename <- file.path(
    sample_subdir,
    ifelse(INTEGRATED & !SUBCLUSTER, sprintf("%s_clustree.png", DATASET_LABEL), sprintf("%s_%s_clustree.png", DATASET_LABEL, sample_name))
  )
  png(filename = clustree_filename, width = 10, height = 10, res = 300, units = 'in')
  print(clustree(seurat_obj, prefix = clusterprefix))
  dev.off()
  # Umaps to check batch effect 
  # TODO set batch labels
  batch_labels <- c("sample") #c("model", "model_biorep")
  if (INTEGRATED) {
    for (batch in batch_labels) {
      p <- DimPlot(seurat_data, reduction = "umap", group.by = batch) +
        ggtitle(sprintf("%s", DATASET_LABEL)) 
      
      dimplot_filename <- file.path(
        sample_subdir,
        ifelse(INTEGRATED & !SUBCLUSTER, sprintf("%s_%s_umap.png", DATASET_LABEL, batch), sprintf("%s_%s_%s_umap.png", DATASET_LABEL, sample_name, batch))
      )
      png(filename = dimplot_filename, width = 10, height = 5, res = 300, units = 'in')
      print(p)
      dev.off()
    }
  }
  for (res in names(clusterlab_list)) {
    clustering <- clusterlab_list[[res]]
    cluster_vec <- seurat_obj[[clustering]] %>%
      rownames_to_column() %>%
      deframe()
    Idents(seurat_obj) <- cluster_vec
    
    # UMAP DimPlot
    p <- DimPlot(seurat_obj, label = TRUE, reduction = "umap") +
      ggtitle(sprintf("%s: (res = %s)", DATASET_LABEL, res)) # TODO may want to include sample name here
    
    dimplot_filename <- file.path(
      sample_subdir,
      ifelse(INTEGRATED & !SUBCLUSTER, sprintf("%s_%s_umap.png", DATASET_LABEL, res), sprintf("%s_%s_%s_umap.png", DATASET_LABEL, sample_name, res))
    )
    png(filename = dimplot_filename, width = 10, height = 5, res = 300, units = 'in')
    print(p)
    dev.off()
    
    # Variable features Heatmap
    topFeatures <- head(VariableFeatures(seurat_obj), TOPN_VARFEAT)
    heatmap_filename <- file.path(
      sample_subdir,
      ifelse(INTEGRATED & !SUBCLUSTER, sprintf("%s_top%s_varmarkers_heatmap_%s.png", DATASET_LABEL, TOPN_VARFEAT, res),
             sprintf("%s_%s_top%s_varmarkers_heatmap_%s.png", DATASET_LABEL, sample_name, TOPN_VARFEAT, res))
    )
    png(filename = heatmap_filename, width = 10, height = 8, res = 300, units = 'in')
    p2 <- DoHeatmap(seurat_obj, features = topFeatures)
    print(p2)
    dev.off()
  }
  
  # Clustree Overlay
  overlay_list <- clustree_overlay(seurat_obj, red_dim = "umap", x_value = "umap1", y_value = "umap2", plot_sides = TRUE)
  
  # ClustreeX
  clustreex_filename <- file.path(
    sample_subdir,
    ifelse(INTEGRATED & !SUBCLUSTER, sprintf("%s_clustreeX.png", DATASET_LABEL), sprintf("%s_%s_clustreeX.png", DATASET_LABEL, sample_name))
  )
  png(filename = clustreex_filename, width = 10, height = 10, res = 300, units = 'in')
  print(overlay_list$x_side)
  dev.off()
  
  # ClustreeY
  clustreey_filename <- file.path(
    sample_subdir,
    ifelse(INTEGRATED & !SUBCLUSTER, sprintf("%s_clustreeY.png", DATASET_LABEL), sprintf("%s_%s_clustreeY.png", DATASET_LABEL, sample_name))
  )
  png(filename = clustreey_filename, width = 10, height = 10, res = 300, units = 'in')
  print(overlay_list$y_side)
  dev.off()
}


# Assuming 'seurat_data' is your Seurat object or list of Seurat objects
if (is.list(seurat_data)) {
  lapply(seurat_data, plot_cluster_res)
} else {
  plot_cluster_res(seurat_data)
}

if (SUBCLUSTER) {
  saveRDS(seurat_data, file = OUTFILE)
} else {
  saveRDS(seurat_data, file = OUTFILE)
}

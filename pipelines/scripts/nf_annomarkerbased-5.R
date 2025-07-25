library(Seurat)
library(tidyverse)
library(ggplot2)
library(openxlsx)
library(HGNChelper)
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
OUTFILE_PRE = get_param_value("--outfilepre")
NORM <- get_param_value("--norm")
INTEGRATION <- get_param_value("--integrate")
SAMPLE_LABEL <- get_param_value("--sample_id")
INTEGRATED = T
if(toupper(INTEGRATION) == "NONE") {
  INTEGRATED = F
}
RES <- as.numeric(get_param_value("--clusterres"))
SUBCLUSTER = get_param_value("--subcluster")
SUBCLUSTERNAME = get_param_value("--subclustername")
SUBCLUSTER = ifelse(SUBCLUSTER == 1, T, F)


PLOTSDIR <- file.path("sctype_anno") # could expose this name as a parameter
if (!dir.exists(PLOTSDIR)) {
  dir.create(PLOTSDIR)
}

MARKERDIR <- file.path("marker_diffgenes") # could expose this name as a parameter
if (!dir.exists(MARKERDIR)) {
  dir.create(MARKERDIR)
}

seurat_data <- readRDS(INFILE)

# get DEGs from clusters --------------------------------------------------
# # install.packages("devtools")
# devtools::install_github("immunogenomics/presto") # for efficiency
getTopDEGs <- function(seurat_object, outfile_prefix = OUTFILE_PRE) {
  if (SUBCLUSTER) {
    sample_name <- SUBCLUSTERNAME
  } else {
    sample_name <- ifelse(INTEGRATION == "ALL", DATASET_LABEL, seurat_object[[SAMPLE_LABEL]][[1]]) 
  }
  cat("Processing: ",sample_name, "\n")
  if (toupper(NORM) =="SCT") {
    seurat_object <- PrepSCTFindMarkers(object = seurat_object)
  }
  sc_markers <- FindAllMarkers(object = seurat_object, only.pos = TRUE, logfc.threshold = 0.25)
  # # conserved markers
  # cluster0_conserved_markers <- FindConservedMarkers(seurat_object,
  #                               ident.1 = 0,
  #                      	      grouping.var = "sample_id",
  #                               only.pos = TRUE,
  # 		              logfc.threshold = 0.25)
  # 
  # avg_log2fc_cols <- grep("_avg_log2FC$", colnames(cluster0_conserved_markers), value = TRUE)
  # 
  # top10 <- cluster0_conserved_markers %>% 
  #   select(all_of(avg_log2fc_cols)) %>%
  #   mutate(avg_fc = rowMeans(.)) %>%
  #   top_n(n = 10, wt = avg_fc)
  # top10
  saveRDS(sc_markers,file.path(MARKERDIR,sprintf("%s_%s_diffgenes.rds", outfile_prefix, sample_name)))
  markercounts <- c(10, 20, 50)
  for (markercount in markercounts) {
    topmarkers <- sc_markers %>% group_by(cluster) %>% top_n(markercount, avg_log2FC)
    write.csv(topmarkers, file.path(MARKERDIR, sprintf("%s_top%s_markers_%s.csv", outfile_prefix, markercount, sample_name)), row.names = FALSE)
  }
}

if (is.list(seurat_data)) {
  lapply(seurat_data, getTopDEGs)
} else {
  getTopDEGs(seurat_data)
}

# SCType Cell annotation ------------------------------------------------------------------

# load gene set preparation function and cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

marker_dir = get_param_value("--markerdir")
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
sctype <- "ScTypeDB_full.xlsx"
cellmarker2 <- "Mouse_cellmarker2_seq_sctype.xlsx" # use tissue ="Breast"
panglao_breast <- "Mm_panglao_sctype_breast.xlsx"
topp_immune <- "Mouse_toppFengshuo_seq_sctype.xlsx" # use tissue ="Immune"
db_local = file.path(marker_dir, cellmarker2) 
db_local2 = file.path(marker_dir, panglao_breast) 
db_local3 = file.path(marker_dir, topp_immune) 
db_local4 = file.path(marker_dir, sctype) 

immune_SCTYPE = gene_sets_prepare(db_local4, "Immune system")
cat(db_local4, "\n")

breast_CELLMARKER2_GENESETS = gene_sets_prepare(db_local, "Breast") # can do db_local
cat(db_local, "\n")

breast_PANGLAO_GENESETS = gene_sets_prepare(db_local2, "Breast") # can do db_local # nolint: line_length_linter.
cat(db_local2, "\n")

fengshuo_TOPP_IMMUNE_GENESETS = gene_sets_prepare(db_local3, "Immune") # this is all tissues
cat(db_local3, "\n")

# get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = seurat_object[["RNA"]]@scale.data, scaled = TRUE, 
#                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
#'gs_list formatting
#'
run_sctype <- function(seurat_object, gs_list = CELL_TYPE_GENESETS, norm_method, is_integrated = F, sample_indicator, short_marker_label = short_marker_label) {
  # TODO method of fetching scaled data is not ideal
  # Determine which assay to pull data from
  data_assay <- if (norm_method == "SCT") {
    "SCT"
  } else if (is_integrated) {
    "integrated"
  } else {
    "RNA"
  }
  
  # Extract appropriate expression matrix
  sc_data <- if (norm_method == "SCT") {
    mat <- seurat_object[[data_assay]]$scale.data
    if (nrow(mat) == 0) stop("SCT data slot is empty. Did you run SCTransform()?")
    mat
  } else {
    mat <- seurat_object[[data_assay]]$scale.data
    if (nrow(mat) == 0) stop(paste0(data_assay, " scale.data slot is empty. Did you run ScaleData()?"))
    mat
  }

  # Run sctype scoring
  es.max <- sctype_score(
    scRNAseqData = sc_data,
    scaled = TRUE,  # Both SCT@data and scale.data are standardized
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
  # In case Seurat is used, it is either seurat_object[["RNA"]]@scale.data (default), seurat_object[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or seurat_object[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
  
  # merge by cluster
  cL_results = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  seurat_object@meta.data[[short_marker_label]] = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat_object@meta.data[[short_marker_label]][seurat_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  # sample_name <- paste(seurat_object@meta.data$model_condition[[1]], seurat_object@meta.data$sample_id[[1]], sep = "-")

  sample_name <- ifelse(INTEGRATION == "ALL", DATASET_LABEL, seurat_object@meta.data[[sample_indicator]][[1]]) 
  if (SUBCLUSTER) {
    sample_name <- SUBCLUSTERNAME
  }
  sample_subdir <- file.path(PLOTSDIR, sample_name) # PLOTTING DIRECTORY separate plots per sample 
  if (!dir.exists(sample_subdir)) {
    dir.create(sample_subdir)
  }
  
  print(DimPlot(seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = short_marker_label) + 
          ggtitle(sprintf("%s %s (res = %s)", sample_name, data_assay,  RES)))
  
  # p1 <- DimPlot(seurat_object, label = T, repel = T, group.by = "seurat_clusters") + ggtitle("Unsupervised clustering")
  p2 <- DimPlot(seurat_object, label = T, repel = T, group.by = short_marker_label) + 
    ggtitle(sprintf("%s %s (res = %s) SCType \n %s", sample_name, data_assay,  RES, short_marker_label))
  
  png(filename = file.path(sample_subdir,paste0(OUTFILE_PRE,"_", sample_name,"_",
                        data_assay,"_sc_umap_labeled_sctype",
                        short_marker_label,".png")), 
      width = 750, height = 500) 
  print(p2)  
  dev.off()
  return(seurat_object)
}

gene_sets <- list(
  SCType      = immune_SCTYPE,
  CellMarker2 = breast_CELLMARKER2_GENESETS,
  Panglao     = breast_PANGLAO_GENESETS,
  ToppCell    = fengshuo_TOPP_IMMUNE_GENESETS
)

if (is.list(seurat_data)) {
  for (label in names(gene_sets)) {
    seurat_data <- lapply(seurat_data, run_sctype,
                          norm_method = NORM,
                          is_integrated = INTEGRATED,
                          gs_list = gene_sets[[label]],
                          sample_indicator = SAMPLE_LABEL,
                          short_marker_label = label)
  }
} else {
  for (label in names(gene_sets)) {
    seurat_data <- run_sctype(seurat_data,
                              norm_method = NORM,
                              is_integrated = INTEGRATED,
                              gs_list = gene_sets[[label]],
                              sample_indicator = SAMPLE_LABEL,
                              short_marker_label = label)
  }
}

saveRDS(seurat_data, file = OUTFILE)


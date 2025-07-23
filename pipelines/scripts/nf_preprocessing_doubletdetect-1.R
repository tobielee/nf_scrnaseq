#!/usr/local/bin R
library(ggplot2)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
# library(stringr)
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
# Retrieve parameter values
DATASET_LABEL <- get_param_value("--datlabel")
OUTSDIR <- get_param_value("--outdir")
OUTSFILE1 <- paste0(DATASET_LABEL,"_filtered_list.rds")
OUTSFILE2 <- paste0(DATASET_LABEL,"_", "detectedDoublets", "_filtered_list.rds")
OUTSFILE3 <- paste0(DATASET_LABEL,"_", "removedDoublets", "_filtered_list.rds")


FILTER_CELLS_MIN <- as.numeric(get_param_value("--fcellmin"))
FILTER_FEAT_MIN <- as.numeric(get_param_value("--ffeatmin"))
FILTER_FEAT_QUANTILE_MIN <- as.numeric(get_param_value("--fquantmin"))
FILTER_FEAT_QUANTILE_MAX <- as.numeric(get_param_value("--fquantmax"))
FILTER_MITO_MAX <- as.numeric(get_param_value("--fmito"))
INDIR <- get_param_value("--datadir")
SPECIES <- get_param_value("--species")
TISSUE <- get_param_value("--tissue")
INPUT_SEURAT <- get_param_value("--inputseurat")
SAMPLE_ID <- get_param_value("--sample_id")

DOUBLETDETECT = get_param_value("--detectdoub")
REMOVEDOUBLETS = get_param_value("--removedoub")
DOUBLETDETECT = ifelse(DOUBLETDETECT == 1, T, F)
REMOVEDOUBLETS = ifelse(REMOVEDOUBLETS == 1, T, F)

PARSING_REGEX = get_param_value("--parseregex")
PARSING_SAMPLESPLIT = get_param_value("--parsesplit")

cat(SPECIES, DOUBLETDETECT, REMOVEDOUBLETS, PARSING_REGEX)


all_files <- list.files(INDIR, full.names = F)
relevant_files <- grep("\\.(mtx|mtx.gz|tsv.gz|tsv)$", all_files, value = TRUE)
parsed_names <- gsub(PARSING_REGEX, "\\1", relevant_files) 
# parsed_names <- str_match(relevant_files, PARSING_REGEX)[,2]

sample_list <- unique(parsed_names)
use_seurat_input <- INPUT_SEURAT != 'null' && nchar(trimws(INPUT_SEURAT)) > 0
if (use_seurat_input) {
  seurat_input <- readRDS(INPUT_SEURAT)
  sample_list <- unique(seurat_input[[SAMPLE_ID]][, 1])
}

mito_10x_pattern <- switch(
  SPECIES,
  "mouse" = "^mt-",
  "rat"   = "^Mt-",
  "human" = "^MT-",
  "^MT-"  # default/fallback
)

seurat_data <- list()
all_files_full <- list.files(INDIR, full.names = T)
relevant_files_full <- grep("\\.(mtx|mtx.gz|tsv.gz|tsv)$", all_files_full, value = TRUE)
# cycle through samples 
for (sample in sample_list) {
  message("Processing sample: ", sample)
  
  if (use_seurat_input) {
    cells_to_keep <- rownames(seurat_input@meta.data)[seurat_input@meta.data[[SAMPLE_ID]] == sample]
    seurat_obj <- subset(seurat_input, cells = cells_to_keep)
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > FILTER_FEAT_MIN)
    
    if (ncol(seurat_obj) <= 10) {
      warning("Skipping sample ", sample, ": too few cells (<= 10)")
      next
    }
  } else {
    # Get file paths for this sample
    filtered_files <- grep(sample, relevant_files_full, value = TRUE)
    message("Matched files: ", toString(filtered_files))
    
    barcodes <- grep("barcodes", filtered_files, value = TRUE)
    features <- grep("genes|features", filtered_files, value = TRUE)
    mtx_file <- grep("\\.mtx", filtered_files, value = TRUE)
    
    # Parse metadata from sample name
    split_result <- strsplit(sample, PARSING_SAMPLESPLIT)[[1]]
    model <- split_result[1]
    condition <- ifelse(length(split_result) > 1, split_result[2], "")
    
    # Read and create Seurat object
    data <- ReadMtx(
      mtx = mtx_file,
      cells = barcodes,
      features = features,
      feature.column = 2
    )
    
    seurat_obj <- CreateSeuratObject(
      counts = data,
      project = sample,
      min.cells = FILTER_CELLS_MIN,
      min.features = FILTER_FEAT_MIN
    )
    
    seurat_obj$sample <- sample
    seurat_obj$model <- model
    seurat_obj$condition <- condition
  }

  seurat_obj$dataset <- DATASET_LABEL
  seurat_obj$tissue <- TISSUE
  seurat_obj$species <- SPECIES
  
  # Calculate percent mitochondrial genes (from raw counts layer)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
    seurat_obj,
    pattern = mito_10x_pattern)
  
  # Compute quantile cutoffs for feature count
  lb <- quantile(seurat_obj$nFeature_RNA, probs = FILTER_FEAT_QUANTILE_MIN, na.rm = TRUE)
  ub <- quantile(seurat_obj$nFeature_RNA, probs = FILTER_FEAT_QUANTILE_MAX, na.rm = TRUE)
  
  # Filter low-quality cells
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > lb &
      nFeature_RNA < ub &
      percent.mt < FILTER_MITO_MAX
  )
  
  seurat_data[[sample]] <- seurat_obj
}
saveRDS(seurat_data, file = OUTSFILE1)

# plot cell counts after preprocessing ------------------------------------
cell_counts <- sapply(seurat_data, function(seurat_object) {
  return(length(Cells(seurat_object)))
})
data <- data.frame(Sample = names(seurat_data), Cell_Count = cell_counts)
count_plot <- ggplot(data, aes(x = Sample, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Number of cells per sample after prelim preprocessing", x = "Sample", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Cell_Count), vjust = -0.3, color = "black", size = 3)

ggsave(file.path(OUTSDIR,paste0(DATASET_LABEL, "_cellcounts.png")), plot = count_plot, width = 10, height = 6, units = "in")

# doublet detection -------------------------------------------------------
preprocessForDoublet <- function(seurat_obj) {
  seurat_obj <- seurat_obj %>%
    # sc[, sc[["nFeature_RNA"]] < ub] %>% # filter outliers with too many features
    NormalizeData %>%
    FindVariableFeatures %>%
    ScaleData %>%
    RunPCA %>%
    RunUMAP(dims = 1:30)
  ElbowPlot(seurat_obj, ndims = 30)
  
  # clustering
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # adjusted res
  
  # visualization
  print(DimPlot(seurat_obj, reduction = "umap"))
  
  return(seurat_obj)
}

# BUG with paramSweep(seu, PCs, sct=FALSE) sct doesn't seem recognized as arg
findDoublets <- function(seurat_obj, PCs = 1:30, pN = 0.25) {
  # Run paramSweep and summarizeSweep
  sweep.res.list <- paramSweep(seurat_obj, PCs = PCs)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  
  # Find optimal pK
  bcmvn <- find.pK(sweep.stats)
  
  # # Plot pK vs BCmetric
  # ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
  #   geom_point() +
  #   geom_line()
  
  # Select the pK that corresponds to max BCmetric
  selected_pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    pull(pK) %>%
    as.numeric()
  print(paste0("selected PK:", selected_pK))
  # Calculate homotypic doublet proportion estimate
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder(seurat_obj,
                              PCs = PCs,
                              pN = pN,
                              pK = selected_pK,
                              nExp = nExp_poi.adj,
                              # reuse.pANN = FALSE,
                              # sct = sct
  )
  return(seurat_obj)
}


if (DOUBLETDETECT) {
  seurat_data <- lapply(seurat_data, preprocessForDoublet)
  seurat_data <- lapply(seurat_data, findDoublets)
  for (sample_name in names(seurat_data)) {
    seurat_obj <- seurat_data[[sample_name]]
    df_classification <- grep("^DF\\.classifications", names(seurat_obj@meta.data), value = TRUE)[1]
    png(file.path(OUTSDIR,paste0(DATASET_LABEL,"_", sample_name, "_doubletplot.png")), width = 750, height = 500)
    print(DimPlot(seurat_obj, group.by = df_classification))
    dev.off()
  }
  saveRDS(seurat_data, file = OUTSFILE2)
  remove_doublets <- function(seurat_obj) {
    df_classification <- grep("^DF\\.classifications", names(seurat_obj@meta.data), value = TRUE)[1]
    seurat_obj@meta.data[["doublet_finder"]] <- seurat_obj@meta.data[[df_classification]] # add new col
    seurat_obj@meta.data[[df_classification]] <- NULL
    return(subset(seurat_obj, doublet_finder == "Singlet"))
  }
  if (REMOVEDOUBLETS) {
    seurat_data <- lapply(seurat_data, remove_doublets)
    saveRDS(seurat_data, file = OUTSFILE3)
  }
}

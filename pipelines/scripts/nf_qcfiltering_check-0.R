
library(ggplot2)
library(Seurat)
library(argparse)
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
FILTER_CELLS_MIN <- as.numeric(get_param_value("--fcellmin"))
FILTER_FEAT_MIN <- as.numeric(get_param_value("--ffeatmin"))
FILTER_FEAT_QUANTILE_MIN <- as.numeric(get_param_value("--fquantmin"))
FILTER_FEAT_QUANTILE_MAX <- as.numeric(get_param_value("--fquantmax"))
FILTER_MITO_MAX <- as.numeric(get_param_value("--fmito"))
INDIR <- get_param_value("--datadir")
SPECIES <- get_param_value("--species")
INPUT_SEURAT <- get_param_value("--inputseurat")
PARSING_REGEX = get_param_value("--parseregex")
sample <- get_param_value("--sample")
cat(sample, "\n")
# # TODO sample parsing
# all_files <- list.files(INDIR, full.names = F)
# relevant_files <- grep("\\.(mtx|mtx.gz|tsv.gz|tsv)$", all_files, value = TRUE)
# parsed_names <- gsub("^[^_]+_(.*)_.*$", "\\1", relevant_files) #captures everything between the first and last underscore in the filename
# sample_list <- unique(parsed_names)

if (SPECIES == "mouse") {
  mito_10x_pattern <- "^mt-"
} else if (SPECIES == "rat") {
  mito_10x_pattern <- "^Mt-"
} else if (SPECIES == "human") {
  mito_10x_pattern <- "^MT-"
} else {
  # Handle the case for other species
  mito_10x_pattern <- "^MT-"
}

if (INPUT_SEURAT != 'null' && nchar(trimws(INPUT_SEURAT)) > 0) {
  seurat_obj <- readRDS(INPUT_SEURAT)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > FILTER_FEAT_MIN)
} else {
  all_files_full <- list.files(INDIR, full.names = T)
  relevant_files_full <- grep("\\.(mtx|mtx.gz|tsv.gz|tsv)$", all_files_full, value = TRUE)
  # cycle through samples 
  print(sample)
  parsed_names <- gsub(PARSING_REGEX, "\\1", basename(relevant_files_full))
  filtered_files <- relevant_files_full[parsed_names == sample]
  print(filtered_files)
  barcodes <- grep("barcodes", filtered_files, value = TRUE)
  features <- grep("genes|features", filtered_files, value = TRUE)
  mtx_file <- grep("\\.mtx", filtered_files, value = TRUE)

  data <- ReadMtx(mtx = mtx_file, cells = barcodes, features = features, feature.column = 2)
  seurat_obj <- CreateSeuratObject(counts = data, project = sample, min.cells = FILTER_CELLS_MIN, min.features = FILTER_FEAT_MIN)
}
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mito_10x_pattern)

lb <- quantile(seurat_obj[["nFeature_RNA"]]$nFeature_RNA, probs = FILTER_FEAT_QUANTILE_MIN)
ub <- quantile(seurat_obj[["nFeature_RNA"]]$nFeature_RNA, probs = FILTER_FEAT_QUANTILE_MAX)
png(paste0(DATASET_LABEL,"_", sample, "_",FILTER_CELLS_MIN,"-",FILTER_FEAT_MIN,"-",FILTER_FEAT_QUANTILE_MIN,"-",FILTER_FEAT_QUANTILE_MAX,"-",FILTER_MITO_MAX,"_qcviolin.png"), width = 750, height = 500)
# Create the VlnPlot
plot1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA"), ncol = 1) +
  geom_hline(yintercept = lb, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = ub, linetype = "dashed", color = "blue") + 
  theme(legend.position = "none")  
plot2 <- VlnPlot(seurat_obj, features = c("nCount_RNA"), ncol = 1) + theme(legend.position = "none")  

plot3 <- VlnPlot(seurat_obj, features = c("percent.mt"), ncol = 1) +
  geom_hline(yintercept = FILTER_MITO_MAX, linetype = "dashed", color = "red") +  # Line at 15 for percent.mt
  theme(legend.position = "none")  
plot <- plot1 | plot2 | plot3
print(plot)
dev.off()





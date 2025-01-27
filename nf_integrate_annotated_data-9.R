# install.packages("remotes")
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# install.packages("SeuratObject")
# install.packages("tidyverse")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
library("tidyverse")
library("Seurat")
library(patchwork)
# library(biomaRt)
source('/Users/tlee/Documents/misc_maybeuseful/rpipeline_singlecell/nextflow/pipeline_utils.R') # includes get_param_value


INT_DIR <- get_param_value("--intoutdir")
INT_LABEL <- get_param_value("--intlabel")

int_files <- list.files(INT_DIR, pattern = "\\_annotated1.rds$", full.names = TRUE)


seurat_list <- list()
for (file in int_files) {
  seurat_obj <- readRDS(file)

  seurat_list <- c(seurat_list, seurat_obj)
}

# Combine all Seurat objects into a single flat list
combined_seurat <- unlist(seurat_list, recursive = FALSE)
print(combined_seurat)



if (INTEGRATION == "NONE") {
  in_file <- paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION,"-", RES,"_list_multi_annotated.rds")
  in_data <- readRDS(file.path(DATADIR, in_file))
  # in_data <- in_data[grep("_negcontrol", names(in_data))] # filter for only neg controls
  out_file <- paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_multi_annotated_",INTEGRATION2,".rds")
  out_path <- paste0(DATADIR, out_file)
} else {
  in_data <- list()
  for (sample in sample_list) {
    RES <- ifelse(sample %in% highres, "1.2","0.8")
    sample_file <-paste0(sample, "-", NORM, "-",INTEGRATION,"-", RES,"_list_multi_annotated.rds")
    sample_data <- readRDS(file.path(DATADIR, sample_file))
    # sample_data[["sample_condition"]] <- sample
    in_data[[sample]] <- sample_data
  }
  out_file <- paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_multi_annotated_",INTEGRATION2,".rds")
  out_path <- file.path(DATADIR, out_file)
  out_file2 <- paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_list_multi_annotated.rds")
  out_path2 <- file.path(DATADIR, out_file2)
}

features <- SelectIntegrationFeatures(object.list = in_data, nfeatures = 2000) # adjust nfeatures for cancer cell comparison?
start <- Sys.time()
anchors <- FindIntegrationAnchors(object.list = in_data, anchor.features = features, reduction = "rpca") # cca reduction is slower
out_data <- IntegrateData(anchorset = anchors, normalization.method =c("LogNormalize"))
print( Sys.time() - start )

# # original unmodified data still resides in the 'RNA' assay
DefaultAssay(out_data) <- "integrated"

# Run the standard workflow for visualization and clustering
out_data <- ScaleData(out_data, verbose = FALSE)
out_data <- RunPCA(out_data, npcs = 30, verbose = FALSE)
out_data <- RunUMAP(out_data, reduction = "pca", dims = 1:30)
out_data <- FindNeighbors(out_data, reduction = "pca", dims = 1:30)
out_data <- FindClusters(out_data, resolution = 2) # adjust

# Visualization
p1 <- DimPlot(out_data, reduction = "umap", group.by = "model")
p2 <- DimPlot(out_data, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(out_data, reduction = "umap", group.by = "broad_cell_type")
# identical(Idents(out_data), out_data$seurat_clusters) # This is true
p1 + p2 + p3

# TODO This must be manually checked
# refine cell type annotation ---------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(UCell)
# display.brewer.all()
# colorRampPalette(brewer.pal(9,"Blues"))(100)


# Calculate relative proportions
proportions <- table(out_data$seurat_clusters, out_data$broad_cell_type) / rowSums(table(out_data$seurat_clusters, out_data$broad_cell_type))
proportions_df <- as.data.frame(proportions)
mycolors <- brewer.pal(length(unique(proportions_df$Var2)), "Set3")
write.csv(proportions_df, file.path(DATADIR,paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_broadcell_proptable_res2.csv")), row.names = FALSE)

result_df <- proportions_df %>%
  group_by(Var1) %>%
  slice(which.max(Freq)) %>%
  ungroup()
write.csv(result_df, file.path(DATADIR,paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_broadcell_filteredproptable_res2.csv")), row.names = FALSE)


# Plot using ggplot2
ggplot(proportions_df, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Relative Proportion of Cell Types in Clusters", x = "Cluster", y = "Relative Proportion") +
  scale_fill_manual(values = mycolors)

# Check NK cell annotations
rat_nk_sig <- c("Klrb1a", "Ncr1", "Klrb1c", "Klrk1", "Ncr3", "Cd244", "Cd226", "Klri2", "Klrc1", "Klrb1b", "Klrb1", "Klre1", "Klri1", "Gzmb", "Prf1")
# check fibroblast mural cell annotations
muhl_fibroblas_sig <- c("Col1a1", "Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln1", "Fbln2", "Cd34", "Pdgfra")
muhl_mural_sig <- c("Des", "Mcam", "Tagln", "Notch3", "Pdgfrb", "Anpep")
# https://www.nature.com/articles/s41467-020-17740-1

sig_list <- list(rat_nk_sig = rat_nk_sig, muhl_fibroblas_sig = muhl_fibroblas_sig, muhl_mural_sig = muhl_mural_sig)
ucell_features <- paste(names(sig_list), "_UCell", sep = "")
for (sample_name in names(in_data)) {
  seurat_obj <- in_data[[sample_name]]
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features = sig_list, ncores = 16)
  for (feat in ucell_features) {
    # VlnPlot for "seurat_clusters"
    vlnplot_clusters <- VlnPlot(seurat_obj, features = feat, group.by = "seurat_clusters")
    ggsave(file.path(DATADIR,paste0(sample_name,"_",feat,"_vlnplot_clusters.png")), vlnplot_clusters, width = 8, height = 6)
    # VlnPlot for "broad_cell_type"
    vlnplot_cell_type <- VlnPlot(seurat_obj, features = feat, group.by = "broad_cell_type")
    ggsave(file.path(DATADIR,paste0(sample_name,"_",feat,"_vlnplot_cell_type.png")), vlnplot_cell_type, width = 8, height = 6)
  }
}

# pericyte check
VlnPlot(out_data, features = c("Pdgfrb", "Rgs5"), group.by = "seurat_clusters")


out_data[["RNA"]] <- JoinLayers(out_data[["RNA"]]) 
DefaultAssay(out_data) <- "RNA"
out_data <- AddModuleScore_UCell(out_data, features = sig_list, ncores = 16)
for (feat in ucell_features) {
  # VlnPlot for "seurat_clusters"
  vlnplot_clusters <- VlnPlot(out_data, features = feat, group.by = "seurat_clusters")
  ggsave(file.path(DATADIR,paste0(DATASET_LABEL,"_",feat,"_vlnplot_clusters.png")), vlnplot_clusters, width = 8, height = 6)
  # VlnPlot for "broad_cell_type"
  vlnplot_cell_type <- VlnPlot(out_data, features = feat, group.by = "broad_cell_type")
  ggsave(file.path(DATADIR,paste0(DATASET_LABEL,"_",feat,"_vlnplot_cell_type.png")), vlnplot_cell_type, width = 8, height = 6)
}
FeaturePlot(out_data, features = "rat_nk_sig_UCell") & scale_color_viridis_c()
DefaultAssay(out_data) <- "integrated"


cluster_labels <- data.frame(Cluster = result_df$Var1,
                             Label = result_df$Var2)
rownames(cluster_labels) <- cluster_labels$Cluster
# Update cluster annotations using Idents
out_data$broad_cell_type2 <- cluster_labels$Label[match(out_data$seurat_clusters, cluster_labels$Cluster)]
# out_data <- readRDS("/Users/tobie/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/xiangzhang_lab/data_to_analyze/rat_proj/rat_human_breast_integrated/xz_ratbreastcancer_models-SCRAN-ALL_multi_annotated_SeuratRpca.rds")
out_data$model_biorep <- paste(out_data$model, out_data$bio_rep, sep = "-")
proportions2 <- table(out_data$model_biorep, out_data$broad_cell_type2) / rowSums(table(out_data$model_biorep, out_data$broad_cell_type2))
proportions_df2 <- as.data.frame(proportions2)
mycolors <- brewer.pal(length(unique(proportions_df2$Var2)), "Set3")
proportions_df2$model <-sapply(strsplit(as.character(proportions_df2$Var1), "-"), function(x) x[1])
facet_props_plot <- ggplot(proportions_df2, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Relative Proportion of Cell Types in Samples", x = "Sample", y = "Relative Proportion") +
  scale_fill_manual(values = mycolors) + 
  theme_bw() +
  facet_grid(~model, scales = "free", switch = "x") + RotatedAxis()
ggsave(file.path(DATADIR,paste0(DATASET_LABEL,"_cellproportion_by_model.png")), facet_props_plot, width = 8, height = 6)

write.csv(proportions_df2, file.path(DATADIR,paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_broadcell2_proptable2_res2.csv")), row.names = FALSE)


broad_cell2_umap <- DimPlot(out_data, reduction = "umap", label = TRUE, repel = TRUE, group.by = "broad_cell_type2")
ggsave(file.path(DATADIR,paste0(DATASET_LABEL,"_broadcell2_umap.png")), broad_cell2_umap, width = 8, height = 6)

# FeaturePlot(out_data, features = "Epcam") & scale_color_viridis_c()

# save updated results ----------------------------------------------------
saveRDS(out_data, out_path)
if (INTEGRATION != "NONE") {
  saveRDS(in_data, out_path2)
}

# TODO save additional for finer annotation
# # IL7R, CCR7
# # IL7R, S100A4
# # CD8A
# # 
# # Cd4
# # VlnPlot(out_data, c("Cd8b", "Il7r"), group.by = "seurat_clusters")
# # 
# # # Get the variable features (genes) from your Seurat object
# # seurat_genes <- VariableFeatures(out_data)
# 
# # Check which genes from your list are present in the Seurat object
# # genes_present <- nk_activ[nk_activ %in% seurat_genes]
# # seurat_genes[grep("Il7", seurat_genes, ignore.case = TRUE)]
# t_cell_subset <- subset(out_data, subset = broad_cell_type == "T Cell" | broad_cell_type2 == "T Cell")
# saveRDS(t_cell_subset, file = file.path(DATADIR, paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_multianno_t_cells.rds")))
# 
# t_cell_subset <- subset(out_data, subset = broad_cell_type == "T Cell" | broad_cell_type2 == "T Cell")
# saveRDS(t_cell_subset, file = file.path(DATADIR, paste0(DATASET_LABEL, "-", NORM, "-",INTEGRATION, "_multianno_t_cells.rds")))
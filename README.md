# scrna-seq-pipeline

# Pipeline details

This is a scRNA-seq analysis annotation pipeline. Downstream analyses to be added later. 

**Current input expected are fastq.**
 
1. Alignment for read counts done using CellRanger (STAR).
```bash
# single sample
cellranger count --id=XZ19335 \
           --transcriptome=/project/zhang/yard/ref/mouse111/Mmusculus_GRCm39_111 \
           --fastqs=/project/zhang/yard/data/fastq/XZ19335 \
           --sample=XZ19335_01,XZ19335_02,XZ19335_03,XZ19335_04

# multi sample
cellranger multi --id=20231124_lung --csv=/project/zhang/yard/nkcell_depletion/20231124_lung/config.csv
```

2. Read counts are copied [nf_cellrangerouts_copy.py]

3. Read counts are QC'd [nf_qcflitering_check-0.R, nf_preprocessing_doubletdetect-1.R]

4. QC'd Read counts are normalized [nf_normintegrate-2.R] 

5. Cells are clustered [nf_clustering-3.R, nf_selectclusterres-4.R]

6. Cells are annotated with marker-based approach e.g. scType [nf_annomarkerbased-5.R] and reference-based approach e.g. scanvi, celltypist [nf_cellannotation_scanvi_celltypist_refbased-6.py]

7. Confirm cluster annotations [nf_plot_cellannotations_eda-7.py], annotate [nf_annotate-8.py], backfill original R objects with annotations [nf_annotations_h5ad_to_rds-8.R]

[below is optional if using multiple datasets] 

8. Integrate annotated data [nf_scanvi_integration-9.py]

9. Refine annotations [nf_intannotation_majority.ipynb]

Break points which require human intervention led me to break the Nextflow pipeline into sections.
0. cellranger run and copy done separately
1. rsinglecellqcfiltering.nf - determining QC parameters for data
2. rsinglecellpipeline.nf - determining clustering resolution to use
3. rsinglecellpipeline2.nf - determining appropriate annotations for individual dataset
4. rsinglecellpipeline3.nf - determining datasets to integrate and potentially reannotate
5. rsinglecellintegrate.nf - integrating datasets

***
## Usage

### Dependencies
* Nextflow
* R
* Python

#### Installing Nextflow

https://www.nextflow.io/docs/latest/install.html

#### Installing Miniforge (for mamba)

https://github.com/conda-forge/miniforge

#### Setting up environment

```bash

# clone repository
git clone <repo>

# install conda environment fresh
mamba create -n scseq python=3.11
conda activate scseq

mamba install -y \
  -c conda-forge \
  -c bioconda \
  -c genomedk \
  -c defaults \
  jupyter jupyterlab \
  scanpy scikit-misc python-igraph leidenalg \
  celltypist anndata2ri bioconductor-ucell bioconductor-singlecellexperiment \
  r-seurat r-clustree r-tidyverse r-openxlsx r-hgnchelper r-doubletfinder r-argparse

# Install pytorch with CUDA 12.1 (adjust as needed)
mamba install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia

# Install JAX with CUDA 12 support
# For CPU-only version, use: pip install -U jax
pip install -U "jax[cuda12]"

# Install scvi-tools with CUDA support
# For CPU-only version, use: pip install -U scvi-tools
pip install -U "scvi-tools[cuda]"
```

**Note**: The `r-doubletfinder` package install via mamba only works on Linux. You can install it manually using the following steps:

```r
install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

```

For GPU accelerated training when using scvi-tools, you will need to install [PyTorch](https://pytorch.org/get-started/locally/) and [JAX](https://jax.readthedocs.io/en/latest/installation.html)

##### Another option to install conda environment

```bash
# gpu accelerated (assumes cuda 12 on linux os)
conda env create -n scseq --file scseq-cuda.yml

```


### Config file setting

There are params which you will need to set. 

#### directory params

`datlabel` = label for identifying individual dataset \
`datadir` = input data directory - this will cellranger output directory \
`outdir` = directory where analysis outputs are stored (i.e. read count QC, clustering analysis, DEG analysis, and cell annotation) \
`intoutdir` = directory for integrative analysis (inputs are assumed directly in this directory and outputs will be stored in this directory) \

#### QC params

`fcellmin` = minimum number of cells expressing a feature (gene) for that feature to be retained for analysis \
`ffeatmin` = minimum number of features present in a cell for that cell to be retained for analysis \
`fquantmin` = quantile minimum to filter out cells whose number of genes (nFeature_RNA) is below this quantile threshold for a given sample \
`fquantmax` = quantile max to filter out cells whose number of genes (nFeature_RNA) is below this quantile threshold for a given sample \
`fmito` = percentage of mitochondrial genes expressed in cell; if cell is above this threshold, it will be removed from analysis \

#### sample metadata params
`species` = label to identify species being analyzed (required for pathway analysis) \
`tissue` = label to identify tissue being analyzed

#### file parsing params
`parseregex` = regular expression used to filter files and identify individual samples for analysis \
`parsesplit` = delimiter for splitting sample annotations if samples have other metadata encoded beyond just sample ID

#### analysis params
`detectdoub` = indicator; 1 means you want to detect doublets, 0 means you don't (required to remove doublets but costs time) \
`removedoub` = indicator;  1 means you want to remove doublets, 0 means you don't (may be useful filter but costs time and additional storage) \
`norm` = types of normalization;  #"SCRAN" # "SCT" # "LOGNORM" \
`integrate` = if you want to integrate  # "ALL", "GROUP" # "NONE" either integrate everything together, just by sample_name, or none \
`npca` = number of princple components to consider prior to clustering \
`clusterreslist` = arraylist represented as a single string containing clustering resolutions to evaluate \
`topn_feat` = number of features to highlight for DEG analysis and top variable gene heatmaps \
`reduction` = dimension reduction for plotting (usually "pca") \
`subcluster` = indicator; 1 means you want to subcluster, 0 means you don't - adjusts how downstream analysis is conducted on a particular cluster (rather than entire dataset) \
`subclustername` = only required with `subcluster` = 1; this serves as a label for the subcluster \
`clusterres` = clustering resolution you find optimal for your dataset \
`markerdir` = directory with marker-based annotation file (this file serves as reference for marker-based cell annotation \
`reflab` = label for input of reference-based annotation also gets added to output file \
`scanvi_model_path` = directory for scanvi model \
`celltypist_model_path` = directory for celltypist model \
`subdirectory_name` = output directory for reference-based annotations \
`sample_id` =  column name for sample identifier for each subject/patient/rat/mouse \
`annofield1` = first pass column name for annotations \
`annofield2` = second pass column name for annotations

 




 
 
/

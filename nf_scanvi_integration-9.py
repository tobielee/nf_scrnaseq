import os
# os.environ['KMP_DUPLICATE_LIB_OK']='True'

import numpy as np
import matplotlib.pyplot as plt

import warnings
import sys
import argparse
import glob
import pandas as pd
import seaborn as sns
import anndata as ad
import scanpy as sc
import scvi
import anndata2ri
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter

import jax
# Print JAX backend
print(f"Using JAX backend: {jax.default_backend()}")

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

scvi.settings.seed = 42

parser = argparse.ArgumentParser(description='scRNA-seq annotate cells - scanvi and celltypist')
parser.add_argument('--datlabel', help='dataset name', required=True)
parser.add_argument('--indir', help='input directory', required=True)
parser.add_argument('--outdir', help='output directory', required=True)
args = parser.parse_args()

# Assign values to variables
DATASET_LABEL = args.datlabel
indir = args.indir
outdir = args.outdir

# create outdir for plots
sampleplots_path = os.path.join(outdir, DATASET_LABEL)
if not os.path.exists(sampleplots_path):
    os.makedirs(sampleplots_path)
sc.settings.figdir = sampleplots_path  # Set the directory to save figures

# load files from the input directory
rds_files = glob.glob(os.path.join(indir, '*.rds')) # TODO in future try to see if these can just be h5ad files... 
adata_list = []
batch_key = "filebatch"
for file in rds_files:
    with localconverter(anndata2ri.converter):
        r_code = f"""
            suppressPackageStartupMessages(library(Seurat))
            mydata <- readRDS('{file}')
            DefaultAssay(mydata) <- "RNA"
            mydata_sce <- as.SingleCellExperiment(mydata, assay = "RNA")
        """
        robjects.r(r_code)
        query_adata = robjects.globalenv['mydata_sce']
        query_adata.layers['counts'] = query_adata.X.copy()
        query_adata.obs[batch_key] = os.path.basename(file)
        adata_list.append(query_adata)

# merge datasets
merged_adata = ad.concat(adata_list, join= "outer") # TODO maybe I should join outer instead especially for disparate data

def convert_cols_tostring(adata):
    for cat in adata.obs.columns:
        if isinstance(adata.obs[cat].values, pd.Categorical):
            pass
        elif pd.api.types.is_float_dtype(adata.obs[cat]):
            pass
        else:
            print(
                f"Setting obs column {cat} (not categorical neither float) to strings to prevent writing error."
            )
            adata.obs[cat] = adata.obs[cat].astype(str)
    for cat in adata.var.columns:
        if isinstance(adata.var[cat].values, pd.Categorical):
            pass
        elif pd.api.types.is_float_dtype(adata.var[cat]):
            pass
        else:
            print(
                f"Setting var column {cat} (not categorical neither float) to strings to prevent writing error."
            )
            adata.var[cat] = adata.var[cat].astype(str)

convert_cols_tostring(merged_adata)
merged_adata.obs.index = merged_adata.obs.index.astype(str) # some annoyance here
merged_adata.var.index = merged_adata.var.index.astype(str) # some annoyance here
merged_adata.raw = merged_adata  # keep full dimension safe I hope...

sc.pp.highly_variable_genes(merged_adata,flavor="seurat_v3", batch_key=batch_key, n_top_genes=6000, subset=True) # this will autosubset for HVG

# train scVI model
try:
    import torch
    torch_available = True
    cuda_available = torch.cuda.is_available()
except ImportError:
    torch = None
    torch_available = False
    cuda_available = False
    print("Torch is not installed. GPU acceleration for PyTorch will not be available.")

# Determine device for PyTorch if installed
if torch_available:
    device = torch.device("cuda" if cuda_available else "cpu")
    if cuda_available:
        print(f"Using GPU for PyTorch: {torch.cuda.get_device_name(0)}")
    else:
        print("Torch is installed, but no GPU available. Using CPU for PyTorch.")
else:
    device = "cpu"  # Default device if torch is not available
    print("Torch is not installed. Using CPU.")

# Replace `merged_adata` with your actual AnnData object
scvi.model.SCVI.setup_anndata(
    merged_adata,
    layer="counts",
    batch_key="filebatch", # TODO need to have this set for each dataset
    # categorical_covariate_keys=["cell_source", "donor"],
    # continuous_covariate_keys=["percent_mito", "percent_ribo"],
)

model = scvi.model.SCVI(merged_adata, n_layers=2, n_latent=30, gene_likelihood="nb")
if torch_available:
    print(f"Using device for PyTorch: {device}")
    model.to_device(device)  # Send model to the selected device if PyTorch is available
model.train()

scvi_model_out = f"{sampleplots_path}/scvi/"
model.save(scvi_model_out, overwrite=False)
SCVI_LATENT_KEY = "X_scVI"
merged_adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

# set unknowns in new column
label_key = "cell_type1_mod" # TODO set
merged_adata.obs["cell_type1_mod"] = merged_adata.obs["cell_type1"]
merged_adata.obs["cell_type1_mod"] = merged_adata.obs["cell_type1_mod"].apply(
    lambda x: "Unknown" if str(x).startswith("Unknown") else x
)

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=merged_adata,
    labels_key="cell_type1_mod", #TODO need to set this
    unlabeled_category="Unknown", #TODO may need to adjust the column to translate all unknowns
)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
SCANVI_LATENT_KEY = "X_scANVI"
merged_adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(merged_adata)

scanvi_model_out = f"{sampleplots_path}/scanvi/"
scanvi_model.save(scanvi_model_out, overwrite=False)
merged_adata.write_h5ad(f"{DATASET_LABEL}_scanvi_merged.h5ad")


# Plotting umaps
LEIDEN_SCVI = "leiden_scVI"
LEIDEN_SCANVI = "leiden_scANVI"
SCVI_UMAP = "_".join([SCVI_LATENT_KEY, "umap"])
SCANVI_UMAP = "_".join([SCANVI_LATENT_KEY, "umap"])

sc.pp.neighbors(merged_adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(merged_adata, key_added=LEIDEN_SCVI, resolution=0.5)
sc.tl.umap(merged_adata)
merged_adata.obsm[SCVI_UMAP] = merged_adata.obsm["X_umap"].copy()

sc.pp.neighbors(merged_adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(merged_adata, key_added=LEIDEN_SCANVI)
sc.tl.umap(merged_adata)
merged_adata.obsm[SCANVI_UMAP] = merged_adata.obsm["X_umap"].copy()
merged_adata.write_h5ad(f"{DATASET_LABEL}_scanvi_merged.h5ad")


def plot_cell_type_proportions(adata, leiden_key, cell_type_key="cell_type1", title_prefix="Proportion of Cell Types in Leiden Clusters", save_path=None):
    # Create a DataFrame with Leiden clusters and cell types
    df = adata.obs[[leiden_key, cell_type_key]].copy()

    # Calculate the proportions of cell types within each Leiden cluster
    proportions = df.groupby([leiden_key, cell_type_key]).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    # Define a color palette
    unique_cell_types = df[cell_type_key].unique()
    palette = sns.color_palette('tab20', len(unique_cell_types))
    color_dict = {cell_type: color for cell_type, color in zip(unique_cell_types, palette)}
    color_dict['Unknown'] = 'black'  # Set 'Unknown' to black

    # Reorder the columns in proportions DataFrame to match color_dict order
    proportions = proportions[color_dict.keys()]

    # Create the bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    proportions.plot(kind='bar', stacked=True, ax=ax, color=[color_dict[cell_type] for cell_type in proportions.columns])
    ax.set_ylabel('Proportion')
    ax.set_title(f'{title_prefix} ({leiden_key})')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Cell Type')

    # Save the plot if save_path is provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    # # Display the plot
    # plt.show()

sc.settings.figdir = sampleplots_path  # Set the directory to save figures
# TODO MDE plotting
# SCVI_MDE_KEY = "X_scVI_MDE"
# merged_adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])
# sc.pl.embedding(merged_adata, basis=SCVI_MDE_KEY,color=["filebatch", "leiden"],ncols=1, frameon=False)
# sc.pl.embedding(merged_adata, basis=SCVI_MDE_KEY, color=["cell_type1"], ncols=1, frameon=False)
# SCANVI_MDE_KEY = "X_scANVI_MDE"
# merged_adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(merged_adata.obsm[SCANVI_LATENT_KEY])
# sc.pl.embedding(merged_adata, basis=SCANVI_MDE_KEY, color=["cell_type1"], ncols=1, frameon=False)

# Plotting for scVI
sc.pl.embedding(merged_adata, basis='X_scVI_umap', color=[LEIDEN_SCVI], ncols=1, frameon=False, legend_loc="on data", save=f"_{LEIDEN_SCVI}_{DATASET_LABEL}.png")
plot_cell_type_proportions(merged_adata, LEIDEN_SCVI, save_path=f"{sampleplots_path}/{LEIDEN_SCVI}_proportions_{DATASET_LABEL}.png")

# Plotting for scANVI
sc.pl.embedding(merged_adata, basis='X_scANVI_umap', color=[LEIDEN_SCANVI], ncols=1, frameon=False, legend_loc="on data", save=f"_{LEIDEN_SCANVI}_{DATASET_LABEL}.png")
plot_cell_type_proportions(merged_adata, LEIDEN_SCANVI, save_path=f"{sampleplots_path}/{LEIDEN_SCANVI}_proportions_{DATASET_LABEL}.png")



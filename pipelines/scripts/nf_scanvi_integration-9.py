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

parser = argparse.ArgumentParser(description='Integrate annotated single-cell data with ScanVI')
parser.add_argument('--datlabel', help='dataset name', required=True)
parser.add_argument('--indir', help='input directory', required=True)
parser.add_argument('--outdir', help='output directory', required=True)
parser.add_argument('--mode', help='file type to load: rds or h5ad', default='h5ad')
parser.add_argument('--annofield1', help='Name of column var for cell annotations', required=True)
parser.add_argument('--annofield2', help='Name of the derived column used as ScanVI training labels', required=True)

args = parser.parse_args()

# Assign values to variables
DATASET_LABEL = args.datlabel
indir = args.indir
outdir = args.outdir
mode = args.mode
manual_anno_label = args.annofield1
manual_anno_label2 = args.annofield2

# create outdir for plots
if not os.path.exists(outdir):
    os.makedirs(outdir)
sc.settings.figdir = outdir  # Set the directory to save figures

adata_list = []
batch_key = "filebatch"  # Currently using one batch per file; could allow user-defined batches in the future

if mode == "rds":
    rds_files = glob.glob(os.path.join(indir, '*.rds'))
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
    file_keys = rds_files

elif mode == "h5ad":
    h5ad_files = glob.glob(os.path.join(indir, '*.h5ad'))
    for file in h5ad_files:
        adata = ad.read_h5ad(file)
        # make sure counts layer is available
        if adata.raw is not None:
            adata.layers['counts'] = adata.raw.X.copy()
        else:
            adata.layers['counts'] = adata.X.copy()
        adata.obs[batch_key] = os.path.basename(file)
        adata_list.append(adata)
    file_keys = h5ad_files

else:
    raise ValueError(f"Unsupported mode: {mode}")

if not adata_list:
    raise FileNotFoundError(
        f"No .{mode} input files found in {indir}; integration needs at least one annotated input."
    )

# merge datasets (identical for both modes)
merged_adata = ad.concat(
    adata_list,
    join="outer",
    # label=batch_key,
    # keys=[os.path.basename(f) for f in file_keys]
)

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

if manual_anno_label not in merged_adata.obs.columns:
    raise ValueError(
        f"Annotation column '{manual_anno_label}' not found in merged input data."
    )
if not manual_anno_label2:
    raise ValueError("The ScanVI label column name (--annofield2) cannot be empty.")
if manual_anno_label2 == manual_anno_label:
    raise ValueError("--annofield2 must differ from --annofield1 so reviewed labels are preserved.")

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
    batch_key=batch_key,
    # categorical_covariate_keys=["cell_source", "donor"],
    # continuous_covariate_keys=["percent_mito", "percent_ribo"],
)

model = scvi.model.SCVI(merged_adata, n_layers=2, n_latent=30, gene_likelihood="nb")
if torch_available:
    print(f"Using device for PyTorch: {device}")
    model.to_device(device)  # Send model to the selected device if PyTorch is available
model.train()

scvi_model_out = f"{outdir}/scvi/"
model.save(scvi_model_out, overwrite=False)
SCVI_LATENT_KEY = "X_scVI"
merged_adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

# Derive a ScanVI label column without modifying the reviewed source labels.
label_key = manual_anno_label2
unknown_prefix = "Unknown"
label_values = merged_adata.obs[manual_anno_label].fillna(unknown_prefix).astype(str).str.strip()
unknown_mask = label_values.str.lower().str.startswith("unknown") | label_values.eq("")
label_values = label_values.mask(unknown_mask, unknown_prefix)
label_categories = sorted(set(label_values) - {unknown_prefix}) + [unknown_prefix]
merged_adata.obs[label_key] = pd.Categorical(label_values, categories=label_categories)

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=merged_adata,
    labels_key=label_key,
    unlabeled_category=unknown_prefix
)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
SCANVI_LATENT_KEY = "X_scANVI"
merged_adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(merged_adata)

scanvi_model_out = f"{outdir}/scanvi/"
scanvi_model.save(scanvi_model_out, overwrite=False)
merged_adata.write_h5ad(f"{outdir}/{DATASET_LABEL}_scanvi_merged.h5ad")


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
merged_adata.write_h5ad(f"{outdir}/{DATASET_LABEL}_scanvi_merged.h5ad")


def plot_cell_type_proportions(adata, leiden_key, cell_type_key=label_key, title_prefix="Proportion of Cell Types in Leiden Clusters", save_path=None):
    # Create a DataFrame with Leiden clusters and cell types
    df = adata.obs[[leiden_key, cell_type_key]].copy()

    # Calculate the proportions of cell types within each Leiden cluster
    proportions = df.groupby([leiden_key, cell_type_key]).size().unstack(fill_value=0)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    # Define a color palette
    unique_cell_types = df[cell_type_key].unique()
    palette = sns.color_palette('tab20', len(unique_cell_types))
    color_dict = {cell_type: color for cell_type, color in zip(unique_cell_types, palette)}
    if unknown_prefix in unique_cell_types:
        color_dict[unknown_prefix] = 'black'

    # Reorder columns to match the observed cell types and their colors.
    proportions = proportions[list(color_dict)]

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

sc.settings.figdir = outdir  # Set the directory to save figures
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
plot_cell_type_proportions(merged_adata, LEIDEN_SCVI, save_path=f"{outdir}/{LEIDEN_SCVI}_proportions_{DATASET_LABEL}.png")

# Plotting for scANVI
sc.pl.embedding(merged_adata, basis='X_scANVI_umap', color=[LEIDEN_SCANVI], ncols=1, frameon=False, legend_loc="on data", save=f"_{LEIDEN_SCANVI}_{DATASET_LABEL}.png")
plot_cell_type_proportions(merged_adata, LEIDEN_SCANVI, save_path=f"{outdir}/{LEIDEN_SCANVI}_proportions_{DATASET_LABEL}.png")



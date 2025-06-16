import os
import warnings
import sys
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import celltypist
from celltypist import models
import pickle
import anndata2ri
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
import argparse

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

scvi.settings.seed = 42
# scanvi_model.registry_["setup_args"]
# convert mixed datatypes into string (for saving)
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
    

### main function for run 
# TODO I probably need to make it more obvious for passing in input params 
def run_scanvi_and_celltypist(adata_query, scanvi_model_path, celltypist_model_path, sample_id, ref_label):
    # Set counts so I don't need to backpropagate in R to get counts
    adata_query.layers["counts"] = adata_query.X.copy() # this might not be necessary for label transfer
    adata_preserved = adata_query.copy()

    # Set values for mapping model prediction
    scvi.model.SCANVI.prepare_query_anndata(adata_query, scanvi_model_path)
    adata_query.obs["mouse.id"] = adata_query.obs[sample_id]
    adata_query.obs["cell_ontology_class"] = "Unknown"
    
    # Load query data for SCANVI          
    # SCANVI labeling - this must go first since resetting pca hampers modele prediction?
    vae_q = scvi.model.SCANVI.load(scanvi_model_path, adata_query)
    # vae_q = scvi.model.SCANVI.load_query_data(adata_query, scanvi_model_path)
    vae_q.train(max_epochs=100, plan_kwargs=dict(weight_decay=0.0), check_val_every_n_epoch=10)
    adata_query.obsm["X_scANVI"] = vae_q.get_latent_representation()
    adata_query.obs[f'scanvipred_{ref_label}'] = vae_q.predict()
    
    # Celltypist
    adata_celltypist = adata_preserved.copy()
    sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
    sc.pp.log1p(adata_celltypist)
    sc.tl.pca(adata_celltypist, n_comps=50)
    
    predictions = celltypist.annotate(adata_celltypist, model=celltypist_model_path, majority_voting=True, mode='prob match', p_thres=0.5)
    adata_preserved = predictions.to_adata(insert_prob=False, prefix=f'celltypist_{ref_label}.')
    
    convert_cols_tostring(adata_preserved)  
    adata_preserved.obs.index = adata_preserved.obs.index.astype(str)
    print(f"ann obj: {adata_preserved}")
    
    # Save query data with predictions
    if DATA_INTEGRATED:
        query_out = outfile

    else:
        sample_name = adata_preserved.obs[sample_id][0]
        directory, filename = os.path.split(outfile)
        query_out = os.path.join(directory, f"{sample_name}_{filename}")

    print(f"Saving {outfile} to query data to {query_out}") # for debugging
    adata_preserved.write_h5ad(query_out)
    return adata_preserved


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='scRNA-seq annotate cells - scanvi and celltypist')
    parser.add_argument('--datlabel', help='dataset name', required=True)
    parser.add_argument('--infile', help='input file', required=True)
    parser.add_argument('--outfile', help='output file', required=True)
    parser.add_argument('--integrate', help='integration method')
    parser.add_argument('--reflab', help='reference label for column annotation in anndata.obs', required=True)
    parser.add_argument('--scanvi_model_path', help='path to scanvi model', required=True)
    parser.add_argument('--celltypist_model_path', help='path to celltypist model', required=True)
    parser.add_argument('--sample_id', help='sample identifier for each subject/patient/rat/mouse')
    parser.add_argument('--subdirectory_name', help='subdirectory for model outs', default='scanvi-celltypist_outs')

    args = parser.parse_args()

    # Assign values to variables
    DATASET_LABEL = args.datlabel
    infile = args.infile
    outfile = args.outfile
    scanvi_model_path = args.scanvi_model_path
    celltypist_model_path = args.celltypist_model_path
    sample_id = args.sample_id
    ref_label = args.reflab
    DATA_INTEGRATED = args.integrate is not None and args.integrate.upper() != 'NONE'
    print(DATA_INTEGRATED)
    subdirectory_name = args.subdirectory_name

    scanvi_model_path = args.scanvi_model_path
    celltypist_model_path = args.celltypist_model_path

    # Create the subdirectory if it doesn't exist
    # output will be saved into the subdirectory by sample if object is a list/dict otherwise it will be a single h5ad
    if not os.path.exists(subdirectory_name):
        os.makedirs(subdirectory_name)

    ### load query data to annotate
    with localconverter(anndata2ri.converter):
        # counts required for training and annotation
        if DATA_INTEGRATED:
            r_code = f"""
                suppressPackageStartupMessages(library(Seurat))
                mydata <- readRDS('{infile}')
                DefaultAssay(mydata) <- "RNA"
                mydata_sce <- as.SingleCellExperiment(mydata, assay = "RNA")
            """
        else:
            r_code = f"""
                suppressPackageStartupMessages(library(Seurat))
                mydata <- readRDS('{infile}')
                mydata_sce <- lapply(mydata, function(sc) {{as.SingleCellExperiment(sc, assay = "RNA")}})
            """
        robjects.r(r_code)
        query_adata = robjects.globalenv['mydata_sce']
    if DATA_INTEGRATED:
        query_adata = run_scanvi_and_celltypist(query_adata, scanvi_model_path, celltypist_model_path, sample_id, ref_label)
    else:
        for i, ann in enumerate(query_adata):
            adata_query = query_adata[ann]
            query_adata[ann] = run_scanvi_and_celltypist(adata_query, scanvi_model_path, celltypist_model_path, sample_id, ref_label)


    ### save list of anndata objects as dictionary pickle
    # if data is not integrated (assumed as a dict/list) 
    if not DATA_INTEGRATED:
        outfilename = outfile.replace(".h5ad", "_alllabels_dict.pkl")
        output_dict_file = os.path.join(outfilename)
        # Assuming you have a list of AnnData objects named 'ann_data_list'
        ann_data_dict = {}
        # Populate the dictionary with AnnData objects
        for ann in query_adata:
            adata_query = query_adata[ann]
            sample_name = adata_query.obs[sample_id][0]
            ann_data_dict[sample_name] = adata_query
            # query_out = f"{data_dir}rat_breast_model_multi_anno/{sample_name}2.h5ad"
            # adata_query.write_h5ad(query_out)
        # save ann_data_dict as pickle
        with open(output_dict_file, 'wb') as file:
            pickle.dump(ann_data_dict, file)

    # if not DATA_INTEGRATED:
    #     # Load the ann_data_dict back using pickle
    #     with open(output_dict_file, 'rb') as file:
    #         ann_data_dict = pickle.load(file)

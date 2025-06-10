import os

import warnings
import sys
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import jaccard_score, adjusted_rand_score
import seaborn as sns
import matplotlib
matplotlib.use('Agg') # avoid using plot display
import pickle
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import itertools
import argparse

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

sc.set_figure_params(figsize=(4, 4))
plt.rcParams['font.size'] = 10  # Adjust font size
plt.rcParams['font.family'] = 'sans-serif'  # Use a simpler font family

# set common params 
plot_params = {'s': 50, 'frameon': False, 'ncols': 4, 'vmax': 'p99', 'show': False}
# List of methods for annotations to compare


def generate_plots(adata_query, sample_name, ref_label, plot_params = plot_params):
    methods = [
    "Panglao",
    "ToppCell",
    "CellMarker2",
    f"scanvipred_{ref_label}",
    f"celltypist_{ref_label}.majority_voting",
    "seurat_clusters",
    # Add more methods as needed
    ]
    
    plots_path = "plots"
    # Create the subdirectory if it doesn't exist
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)

    sampleplots_path = os.path.join(plots_path, sample_name)
    if not os.path.exists(sampleplots_path):
        os.makedirs(sampleplots_path)
    sc.settings.figdir = sampleplots_path  # Set the directory to save figures

    fig, axes = plt.subplots(3, 2, figsize=(20, 18)) # adjust size as needed
    fig.tight_layout()
    plt.subplots_adjust(wspace=1)
    axes = axes.flatten()

    # check samples batch correction 
    try:
        sc.pl.umap(adata_query, color="sample")
    except RuntimeError as e:
        print(f"RuntimeError during plotting 'sample' UMAP: {e}")   

    # plot annotations for each method
    for method, ax in zip(methods, axes):  
        try:
            if method == "seurat_clusters":
                sc.pl.umap(adata_query, color=method, ax=ax, **plot_params, legend_loc='on data')
            else:
                sc.pl.umap(adata_query, color=method, ax=ax, **plot_params)
        except RuntimeError as e:
            print(f"RuntimeError during plotting '{method}' UMAP: {e}")
     # save figure
    umap_out = os.path.join(sampleplots_path, f"{sample_name}_multianno_umap.png")
    fig.savefig(umap_out, bbox_inches='tight') 
    plt.close(fig)
    
    try:
        sc.pl.umap(adata_query, color='seurat_clusters', s=50, frameon=False, ncols=4, vmax='p99', legend_loc='on data', save=f"_seuratclusters_{sample_name}.png") # need to make this not the highest res - default is highest res
        plt.close()
    except RuntimeError as e:
        print(f"RuntimeError during plotting 'seurat_clusters': {e}")

    # # marker feature plots
    # macrophage_markers2 = ['Cd68',
    #                   'Naaa',
    #                   'Tyrobp',
    #                   'Lyz2',
    #                   'Adgre1',
    #                   'Ccl6',
    #                   'Ptprc',
    #                   'Csf1r',
    #                   'Aif1',
    #                   'Itgb2',
    #                  ]

    # macrophage_markers = [
    #     'Cd68', 'Naaa', 'Tyrobp', 'Lyz2', 'Adgre1', 'Ccl6', 'Ptprc', 'Csf1r', 'Aif1', 'Itgb2',
    #     'Ctss', 'Lgals3', 'Fn1', 'Ifi27l2a', 'Tgfbi', 'Clec7a', 'Ms4a6c', 'Ifi30', 'Mafb', 'Ms4a4c', 'Ccl2', 'S100a4'
    # ]
    # macrophage_markers =  list(set(macrophage_markers).intersection(adata_query.var_names))

    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=macrophage_markers, s=50, frameon=False, ncols=4, vmax='p99', save=f"_macrophagemarks_{sample_name}.png")
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting macrophage markers: {e}")
    #     plt.close()

    # neutrophil_markers_combined = [
    #     'S100a8', 'S100a9', 'Camp', 'Cd177', 'Cxcr2', 'Slpi', 'Slfn4', 'Mmp8',
    #     'Ptprc', 'S100a8', 'S100a9', 'Csf3r', 'Cxcr2', 'Lrg1'
    # ]
    # neutrophil_markers_combined =  list(set(neutrophil_markers_combined).intersection(adata_query.var_names))
    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=neutrophil_markers_combined, s=50, frameon=False, ncols=4, vmax='p99', save=f"_neutrophilmarks_{sample_name}.png")
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting neutrophil markers: {e}")
    #     plt.close()
    
    # ##### additional markers
    # stromal_markers = ['Pecam1', 'Acta2', 'Myl9', 'Mylk']
    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=stromal_markers, title=f'{sample_name} stromal markers', s=50, frameon=False, ncols=4, vmax='p99', save=f"_stromalmarkers_{sample_name}.png")
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting stromal markers: {e}")
    #     plt.close()
    
    
    # rat_nk_sig = ["Klrb1a", "Ncr1", "Klrb1c", "Klrk1", "Ncr3", "Cd244", "Cd226", "Klri2", "Klrc1", "Klrb1b", "Klrb1", "Klre1", "Klri1", "Gzmb", "Prf1"]
    # rat_nk_sig =  list(set(rat_nk_sig).intersection(adata_query.var_names))
    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=rat_nk_sig, title=f'{sample_name} nk markers', s=50, frameon=False, ncols=4, vmax='p99', save=f"_ratnkmarkers_{sample_name}.png")    
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting rat nk markers: {e}")
    #     plt.close()

    # nk_markers = ['Nkg7',
    #             'Klrb1c', 
    #             'Klrk1', 
    #             'Ncr1', 
    #             'Klrg1', 
    #             'Klrb1b', 
    #             'Eomes', 
    #             'Prf1', 
    #             'Gzma', 
    #             'Il18rap', 
    #             'S1pr5', 
    #             'Nkg7', 
    #             'Klrd1', 
    #             'Dock2',
    #             'Xcl1',
    #             'Ncr1']
    # nk_markers_2 = ['Ncr1', 'Nkg7', 'Xcl1', 'Gnly', 'Klrd1', 'Ncam1', 'Gzma', 'Ytgae', 'Havr2c', 'Klrb1c', 'Klrb1', 'Klrc1', 'Cd3g', 'Ccr7', 'Cd3d', 'Cd8a']
    # nk_markers = nk_markers + nk_markers_2
    # nk_markers =  list(set(nk_markers).intersection(adata_query.var_names))
    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=nk_markers, s=50, frameon=False, ncols=4, vmax='p99', save=f"_nkcellmarks_{sample_name}.png")
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting nk cell markers: {e}")
    #     plt.close()

    # pericyte_markers = [
    #     "Pdgfrb",  # PDGFRB
    #     "Aldh1a1",  # ALDH1A1
    #     "Cd44",  # CD44
    #     "Cspg4",  # CSPG4
    #     "Rgs5",  # RGS5
    #     "Cd36"  # CD36
    # ]
    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=pericyte_markers, s=50, frameon=False, ncols=4, vmax='p99', save=f"_pericytemarks_{sample_name}.png")
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting pericyte markers: {e}")
    #     plt.close()
    
    # fibroblast = ['Col1a1',
    #           'Col1a2',
    #           'Col5a1',
    #           'Loxl1',
    #           'Lum',
    #           'Fbln1',
    #           'Fbln2',
    #           'Cd34',
    #           'Pdgfra',
    #           'Fstl1',
    #           'Bgn',
    #           'Col16a1',
    #           'Mmp2',
    # ]

    # with rc_context({'figure.figsize': (3, 3)}):
    #     try:
    #         sc.pl.umap(adata_query, color=fibroblast, s=50, frameon=False, ncols=4, vmax='p99')
    #     except RuntimeError as e:
    #         print(f"RuntimeError during plotting fibroblast markers: {e}")
    #     plt.close()

    # mouse_tumor_cell_markers = [
    #     'Epcam',
    #     'Krt8',
    #     'Krt18',
    #     'Krt17',
    #     'Krt19',
    #     'Tp53',
    #     'Kras',
    #     'Nras',
    #     'Hras',
    #     'Mki67',
    #     'Pcna',
    #     'Cd44',
    #     'Aldh1a1',
    #     'Ccnb1',
    #     'Ccnb2',
    #     'Cdc20',
    #     'Cdc25a',
    #     # Add more mouse tumor cell markers as needed
    # ]
    # mouse_tumor_cell_markers =  list(set(mouse_tumor_cell_markers).intersection(adata_query.var_names))

    # with rc_context({'figure.figsize': (3, 3)}):
    #     sc.pl.umap(adata_query, color= mouse_tumor_cell_markers, s=50, frameon=False, ncols=4, vmax='p99', save = f"cancercellmarkers_featureplot_{sample_name}.png")
    #     plt.close()
    # mouse_immune_cell_markers = [
    #     'Ptprc',
    #     'Cd3d',
    #     'Ly6c',
    #     'Ly6g',
    #     'Cd3e',
    #     'Cd3g',
    #     'Cd4',
    #     'Cd8a',
    #     'Cd8b',
    #     'Trac',
    #     'Trbc1',
    #     'Cd19',
    #     'Ms4a1',
    #     'Cd79a',
    #     'Cd79b',
    #     'Cd14',
    #     'Cd68',
    #     'Itgam',
    #     'Cd56',
    #     'Klrd1',
    #     'Ncr1',
    #     'Itgax',
    #     'Bdca1',
    #     'Bdca3',
    #     'Kit',
    #     'Fcer1a',
    #     'S100a8', 
    #     'S100a9',
    #     # Add more mouse immune cell markers as needed
    # ]
    # mouse_immune_cell_markers =  list(set(mouse_immune_cell_markers).intersection(adata_query.var_names))
    # with rc_context({'figure.figsize': (3, 3)}):
    #     sc.pl.umap(adata_query, color= mouse_immune_cell_markers, s=50, frameon=False, ncols=4, vmax='p99', save = f"immunecellmarkers_featureplot_{sample_name}.png")
    #     plt.close()
    # monocyte_genes = ['Ly6c2', 'Vcan', 'Fn1']
    # with rc_context({'figure.figsize': (3, 3)}):
    #     sc.pl.umap(adata_query, color= monocyte_genes, s=50, frameon=False, ncols=4, vmax='p99')
    #     plt.close()
    # gmp_genes = ['Kit', 'Cd34', 'Csf1r', 'Ly6c1', 'Ly6c2', 'Fcgr3','Fcgr2b']
    # mdp_genes = ['Csf1r', 'Cx3cr1', 'Flt3']
    # lineage_markers = ['Cd3g', 'Cd3d', 'Cd3e','Cd8a', 'Cd19', 'Ly6g6c', 'Ly6g6d', 'Ly6g6e', 'Ly6g6f', 'Ly6g5c', 'Ly6g5b', 'Itgam']
    # with rc_context({'figure.figsize': (3, 3)}):
    #     sc.pl.umap(adata_query, color= mdp_genes, s=50, frameon=False, ncols=4, vmax='p99', save=f"_mdpmarks_{sample_name}.png")
    #     plt.close()
    # with rc_context({'figure.figsize': (3, 3)}):
    #     sc.pl.umap(adata_query, color= gmp_genes, s=50, frameon=False, ncols=4, vmax='p99', save=f"_gmpmarks_{sample_name}.png")
    #     plt.close()
    # with rc_context({'figure.figsize': (3, 3)}):
    #     sc.pl.umap(adata_query, color= lineage_markers, s=50, frameon=False, ncols=4, vmax='p99', save=f"_linmarks_{sample_name}.png")
    #     plt.close()
    # write csv for annotations
    adata_query.obs["seurat_clusters"] = adata_query.obs["seurat_clusters"].astype(str)  # In case they are not strings
    unique_clusters = sorted(adata_query.obs["seurat_clusters"].unique(), key=lambda x: int(x))

    # Create a DataFrame with seurat_cluster and cell_annotation columns
    df = pd.DataFrame({
        'seurat_cluster': unique_clusters,
        'cell_annotation': ''
    })
    csv_filename = f'{sample_name}_annotations.csv'
    df.to_csv(csv_filename, index=False)
        
if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='scRNA-seq annotate cells - scanvi and celltypist')
    parser.add_argument('--datlabel', help='dataset name', required=True)
    parser.add_argument('--infile', help='input file', required=True)
    parser.add_argument('--integrate', help='integration method')
    parser.add_argument('--sample_id', help='sample identifier for each subject/patient/rat/mouse')
    parser.add_argument('--reflab', help='reference label for column annotation in anndata.obs', required=True)

    args = parser.parse_args()

    # Assign values to variables
    DATASET_LABEL = args.datlabel
    infile = args.infile
    sample_id = args.sample_id
    ref_label = args.reflab

    DATA_INTEGRATED = args.integrate is not None and args.integrate.upper() != 'NONE'
    print(DATA_INTEGRATED)

    # load data
    if DATA_INTEGRATED:
        mydata = sc.read_h5ad(infile)
        adata_query = mydata
        generate_plots(adata_query,DATASET_LABEL, ref_label)
    else:
        # Load the pickle file containing the dictionary of AnnData objects
        with open(infile, 'rb') as file:
            mydata = pickle.load(file)
        
        # Ensure that mydata is a dictionary of AnnData objects
        if isinstance(mydata, dict):
            for sample_name, adata_query in mydata.items():
                generate_plots(adata_query, sample_name, ref_label)
        else:
            raise ValueError("Expected a dictionary of AnnData objects in the pickle file")
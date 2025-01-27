import argparse
import scanpy as sc
import pandas as pd

# Function to find the last non-blank value in a column
def get_last_non_blank(series):
    non_blank_series = series.dropna().astype(str)
    non_blank_series = non_blank_series[non_blank_series.str.strip() != '']
    return non_blank_series.iloc[-1] if not non_blank_series.empty else None


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='scRNA-seq annotate cells - scanvi and celltypist')
    parser.add_argument('--annocsv', help='CSV with annotations', required=True)
    parser.add_argument('--adatafile', help='input anndata file (.h5ad)', required=True)
    parser.add_argument('--annofield1', help='Name of column var for cell annotations', required=True)

    args = parser.parse_args()

    # Assign values to variables
    annocsv = args.annocsv
    infile = args.adatafile
    manual_anno_label = args.annofield1

    # read in the manual annotation file
    df = pd.read_csv(annocsv)
    updated_cluster_mapping = dict(zip(df['seurat_cluster'].astype(str), df['cell_annotation']))

    # Create the metadata_mapping for all other columns
    metadata_mapping = {}
    for col in df.columns:
        if col not in ['seurat_cluster', 'cell_annotation']:
            metadata_mapping[col] = get_last_non_blank(df[col])

    print("Updated Cluster Mapping:")
    print(updated_cluster_mapping)
    print("\nMetadata Mapping:")
    print(metadata_mapping)

    # Update the adata object with the new annotations
    adata = sc.read_h5ad(infile)
    adata.obs[manual_anno_label] = [updated_cluster_mapping[cluster] for cluster in adata.obs['seurat_clusters']]
    adata.obs[manual_anno_label].head()

    for key, value in metadata_mapping.items():
        adata.obs[key] = value

    outfile = infile.replace('.h5ad', '_annotated1.h5ad')
    print(f"Saving {outfile} - with annotations") # for debugging

    adata.write_h5ad(outfile)
    
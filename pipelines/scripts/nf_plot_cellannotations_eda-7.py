"""
Generate multi-annotation UMAP panels and a blank annotation CSV
for downstream manual labeling.

Prefers legend_loc="on data" (with truncation) for easier cell annotation.
Falls back gracefully if rendering fails.
"""

import os
import warnings
import argparse
import pickle

import anndata
import scanpy as sc
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------

warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)
warnings.simplefilter("ignore", category=UserWarning)

sc.set_figure_params(figsize=(4, 4))

plt.rcParams["font.size"] = 10
plt.rcParams["font.family"] = "sans-serif"


# Default scatter parameters for UMAP plots
PLOT_PARAMS = {
    "s": 50,
    "frameon": False,
    "show": False,
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _safe_cluster_sort(clusters):
    """
    Sort cluster labels numerically when possible, while allowing
    non-numeric or mixed cluster labels.
    """

    def key_fn(x):
        try:
            return (0, int(x))
        except (ValueError, TypeError):
            return (1, str(x))

    return sorted(clusters, key=key_fn)


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def generate_plots(
    adata_query: anndata.AnnData,
    sample_name: str,
    ref_label: str,
    plot_params: dict = None,
):
    """
    Generate:
      1. A multi-panel UMAP comparing annotation methods
         (prefers "on data" legends with truncation).
      2. A dedicated Seurat-cluster UMAP with a right-margin legend.
      3. A blank CSV for manual cluster annotation.
    """

    if plot_params is None:
        plot_params = PLOT_PARAMS.copy()

    # ------------------------------------------------------------------
    # Determine which annotation columns are available
    # ------------------------------------------------------------------

    preferred_methods = [
        "Panglao",
        "ToppCell",
        "CellMarker2",
        f"scanvipred_{ref_label}",
        f"celltypist_{ref_label}.majority_voting",
        "seurat_clusters",
    ]

    methods = [
        method
        for method in preferred_methods
        if method in adata_query.obs.columns
    ]

    if not methods:
        print(f"[WARN] No annotation columns found for {sample_name}")
        return

    print(f"Generating plots for: {sample_name}")
    print(f"Annotation columns: {methods}")

    # ------------------------------------------------------------------
    # Output directories
    # ------------------------------------------------------------------

    plots_root = "plots"
    sample_dir = os.path.join(plots_root, sample_name)
    os.makedirs(sample_dir, exist_ok=True)

    # Scanpy uses this directory for plots generated with save=
    sc.settings.figdir = sample_dir

    # ------------------------------------------------------------------
    # 1. Multi-panel annotation UMAP
    #    Prefer "on data" + truncation for easier cell annotation.
    # ------------------------------------------------------------------

    n = len(methods)
    ncols = 2
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(14, 5.5 * nrows),
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).flatten()

    # Optional sample/batch UMAP check
    if "sample" in adata_query.obs.columns:
        try:
            sc.pl.umap(adata_query, color="sample", show=False)
            plt.close()
        except Exception as e:
            print(f"[WARN] Sample UMAP failed: {e}")

    for ax, method in zip(axes, methods):
        try:
            # Truncate long category names to reduce raster-overflow risk
            col = adata_query.obs[method].astype(str)
            max_len = col.str.len().max()

            if max_len > 22:
                truncated = col.str.slice(0, 20) + "…"
                adata_query.obs["_tmp_color"] = truncated
                color_key = "_tmp_color"
            else:
                color_key = method

            n_cats = adata_query.obs[color_key].nunique(dropna=True)

            # Prefer labels on the data when the number of categories is manageable
            if n_cats <= 30:
                legend = "on data"
                legend_fontsize = 6
            else:
                legend = "right margin"
                legend_fontsize = 7

            sc.pl.umap(
                adata_query,
                color=color_key,
                ax=ax,
                legend_loc=legend,
                legend_fontsize=legend_fontsize,
                **plot_params,
            )
            ax.set_title(method, fontsize=11)

            # Clean up temporary column
            if color_key == "_tmp_color":
                del adata_query.obs["_tmp_color"]

        except Exception as e:
            print(f"[WARN] UMAP for '{method}' failed ({e}). Trying without legend…")
            try:
                # Clean up any temporary column that may have been left behind
                if "_tmp_color" in adata_query.obs.columns:
                    del adata_query.obs["_tmp_color"]

                sc.pl.umap(
                    adata_query,
                    color=method,
                    ax=ax,
                    legend_loc=None,
                    **plot_params,
                )
                ax.set_title(method, fontsize=11)
            except Exception as e2:
                print(f"[WARN] Even no-legend plot failed: {e2}")
                ax.set_visible(False)

    # Hide unused axes
    for ax in axes[len(methods):]:
        ax.set_visible(False)

    out_png = os.path.join(sample_dir, f"{sample_name}_multianno_umap.png")
    fig.savefig(out_png, bbox_inches="tight", dpi=150)
    plt.close(fig)

    print(f"Saved multi-panel UMAP: {out_png}")

    # ------------------------------------------------------------------
    # 2. Dedicated Seurat-cluster UMAP (right-margin legend)
    # ------------------------------------------------------------------

    if "seurat_clusters" in adata_query.obs.columns:
        try:
            sc.pl.umap(
                adata_query,
                color="seurat_clusters",
                legend_loc="right margin",
                save=f"_seuratclusters_{sample_name}.png",
                **plot_params,
            )
            plt.close()
            print(f"Saved Seurat-cluster UMAP for {sample_name}")

        except Exception as e:
            print(f"[WARN] Dedicated seurat_clusters plot failed: {e}")

    # ------------------------------------------------------------------
    # 3. Blank annotation CSV for manual labeling
    # ------------------------------------------------------------------

    if "seurat_clusters" in adata_query.obs.columns:
        clusters = (
            adata_query.obs["seurat_clusters"]
            .astype(str)
            .unique()
        )
        unique_clusters = _safe_cluster_sort(clusters)

        df = pd.DataFrame(
            {
                "seurat_cluster": unique_clusters,
                "cell_annotation": "",
            }
        )

        csv_path = f"{sample_name}_annotations.csv"
        df.to_csv(csv_path, index=False)
        print(f"Saved blank annotation table: {csv_path}")


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "scRNA-seq multi-annotation UMAP plots "
            "(ScanVI / CellTypist / other methods)"
        )
    )

    parser.add_argument(
        "--datlabel",
        required=True,
        help="Dataset name",
    )

    parser.add_argument(
        "--infile",
        required=True,
        help="Input .h5ad or .pkl file",
    )

    parser.add_argument(
        "--integrate",
        default=None,
        help="Integration method",
    )

    parser.add_argument(
        "--sample_id",
        default=None,
        help="Sample identifier column (currently unused, reserved for future use)",
    )

    parser.add_argument(
        "--reflab",
        required=True,
        help="Reference label used in annotation column names",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Determine whether input is an integrated h5ad or a pickle
    # ------------------------------------------------------------------

    DATA_INTEGRATED = (
        args.integrate is not None
        and str(args.integrate).upper() != "NONE"
    )

    print(f"Integrated data: {DATA_INTEGRATED}")

    # ------------------------------------------------------------------
    # Integrated data: one AnnData object
    # ------------------------------------------------------------------

    if DATA_INTEGRATED:
        adata = sc.read_h5ad(args.infile)
        generate_plots(adata, args.datlabel, args.reflab)

    # ------------------------------------------------------------------
    # Non-integrated data: dictionary of AnnData objects
    # ------------------------------------------------------------------

    else:
        with open(args.infile, "rb") as fh:
            data = pickle.load(fh)

        if not isinstance(data, dict):
            raise ValueError(
                "Expected a dict of AnnData objects in the pickle file"
            )

        for sample_name, adata in data.items():
            generate_plots(adata, sample_name, args.reflab)
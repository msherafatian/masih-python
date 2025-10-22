"""
Scanpy utility functions for data processing.

This module provides helper functions for single-cell RNA-seq analysis
using Scanpy, matching the Seurat workflow from the R version.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Optional, Dict, Any, List
import warnings

from masih.config.settings import AppConfig


def setup_scanpy_settings():
    """Configure Scanpy settings for consistent output."""
    sc.settings.verbosity = 1  # Errors only
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    # Set random seed for reproducibility
    np.random.seed(AppConfig.RANDOM_STATE)


def check_adata_analyses(adata) -> Dict[str, bool]:
    """
    Check which analyses have been performed on an AnnData object.

    Equivalent to check_seurat_analyses in R.

    Args:
        adata: AnnData object

    Returns:
        Dictionary of boolean flags for each analysis
    """
    analyses = {
        'normalized': False,
        'log1p': False,
        'highly_variable_genes': False,
        'scaled': False,
        'pca': False,
        'neighbors': False,
        'umap': False,
        'tsne': False,
        'leiden': False,
        'louvain': False,
        'clusters': False
    }

    # Check for normalization (presence of normalized counts)
    if adata.X is not None:
        # Check if data is normalized (not all integers)
        sample_data = adata.X[:100, :100] if hasattr(
            adata.X, 'toarray') else adata.X[:100, :100]
        if hasattr(sample_data, 'toarray'):
            sample_data = sample_data.toarray()
        analyses['normalized'] = not np.allclose(
            sample_data, sample_data.astype(int), rtol=0, atol=1e-10)

    # Check for log transformation
    if 'log1p' in adata.uns:
        analyses['log1p'] = True

    # Check for highly variable genes - must have the column AND some genes marked as True
    if 'highly_variable' in adata.var.columns:
        n_hvg = adata.var['highly_variable'].sum()
        analyses['highly_variable_genes'] = n_hvg > 0

    # Check for scaled data
    if 'scale.data' in adata.layers or (adata.X is not None and np.abs(adata.X.mean()) < 0.1):
        analyses['scaled'] = True

    # Check for PCA
    if 'X_pca' in adata.obsm:
        analyses['pca'] = True

    # Check for neighbors
    if 'neighbors' in adata.uns:
        analyses['neighbors'] = True

    # Check for UMAP
    if 'X_umap' in adata.obsm:
        analyses['umap'] = True

    # Check for t-SNE
    if 'X_tsne' in adata.obsm:
        analyses['tsne'] = True

    # Check for Leiden clustering
    if 'leiden' in adata.obs.columns:
        analyses['leiden'] = True
        analyses['clusters'] = True

    # Check for Louvain clustering
    if 'louvain' in adata.obs.columns:
        analyses['louvain'] = True
        analyses['clusters'] = True

    return analyses


def calculate_qc_metrics(adata) -> None:
    """
    Calculate QC metrics including mitochondrial percentage.

    Modifies adata in place.

    Args:
        adata: AnnData object
    """
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )


def filter_cells_qc(adata, min_genes: int = 200, max_genes: int = 2500,
                    min_counts: int = 200, max_counts: int = None,
                    max_mt_percent: float = 5) -> tuple:
    """
    Filter cells based on QC metrics.

    Args:
        adata: AnnData object
        min_genes: Minimum number of genes per cell
        max_genes: Maximum number of genes per cell
        min_counts: Minimum number of counts per cell
        max_counts: Maximum number of counts per cell (None for no limit)
        max_mt_percent: Maximum mitochondrial percentage

    Returns:
        Tuple of (filtered adata, n_cells_before, n_cells_after)
    """
    n_cells_before = adata.n_obs

    # Filter based on QC metrics
    filters = (
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['pct_counts_mt'] < max_mt_percent)
    )

    if max_counts is not None:
        filters = filters & (adata.obs['total_counts'] <= max_counts)

    adata = adata[filters, :].copy()

    n_cells_after = adata.n_obs

    return adata, n_cells_before, n_cells_after


def normalize_data(adata, target_sum: float = 1e4) -> None:
    """
    Normalize data (equivalent to NormalizeData in Seurat).

    Args:
        adata: AnnData object
        target_sum: Target sum for normalization
    """
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


def find_variable_genes(adata, n_top_genes: int = 3000) -> None:
    """
    Find highly variable genes (equivalent to FindVariableFeatures).

    Args:
        adata: AnnData object
        n_top_genes: Number of highly variable genes to identify
    """
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor='seurat_v3',
        subset=False
    )


def scale_data(adata, max_value: float = 10) -> None:
    """
    Scale data (equivalent to ScaleData in Seurat).

    Note: This will densify sparse matrices. For large datasets,
    consider using backed mode or only scaling highly variable genes.

    Args:
        adata: AnnData object
        max_value: Maximum value after scaling
    """
    # Scale all genes (simpler and matches Seurat's ScaleData behavior)
    sc.pp.scale(adata, max_value=max_value)


def run_pca(adata, n_comps: int = 50, use_highly_variable: bool = True,
            random_state: int = AppConfig.RANDOM_STATE) -> None:
    """
    Run PCA (equivalent to RunPCA in Seurat).

    Args:
        adata: AnnData object
        n_comps: Number of principal components
        use_highly_variable: Whether to use only highly variable genes
        random_state: Random seed for reproducibility
    """
    # Use updated parameter name to avoid deprecation warning
    if use_highly_variable and 'highly_variable' in adata.var.columns:
        mask_var = 'highly_variable'
    else:
        mask_var = None

    sc.tl.pca(
        adata,
        n_comps=n_comps,
        mask_var=mask_var,
        random_state=random_state
    )


def compute_neighbors(adata, n_neighbors: int = 30, n_pcs: int = 30,
                      random_state: int = AppConfig.RANDOM_STATE) -> None:
    """
    Compute neighborhood graph (equivalent to FindNeighbors).

    Args:
        adata: AnnData object
        n_neighbors: Number of neighbors
        n_pcs: Number of PCs to use
        random_state: Random seed for reproducibility
    """
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        random_state=random_state
    )


def run_umap(adata, random_state: int = AppConfig.RANDOM_STATE) -> None:
    """
    Run UMAP (equivalent to RunUMAP in Seurat).

    Args:
        adata: AnnData object
        random_state: Random seed for reproducibility
    """
    sc.tl.umap(adata, random_state=random_state)


def run_tsne(adata, n_pcs: int = 30,
             random_state: int = AppConfig.RANDOM_STATE) -> None:
    """
    Run t-SNE (equivalent to RunTSNE in Seurat).

    Args:
        adata: AnnData object
        n_pcs: Number of PCs to use
        random_state: Random seed for reproducibility
    """
    sc.tl.tsne(adata, n_pcs=n_pcs, random_state=random_state)


def find_clusters(adata, resolution: float = 0.8,
                  random_state: int = AppConfig.RANDOM_STATE,
                  algorithm: str = 'leiden') -> None:
    """
    Find clusters (equivalent to FindClusters in Seurat).

    Args:
        adata: AnnData object
        resolution: Clustering resolution
        random_state: Random seed for reproducibility
        algorithm: Clustering algorithm ('leiden' or 'louvain')
    """
    if algorithm == 'leiden':
        sc.tl.leiden(
            adata,
            resolution=resolution,
            random_state=random_state
        )
    else:
        sc.tl.louvain(
            adata,
            resolution=resolution,
            random_state=random_state
        )


def score_cell_cycle(adata, s_genes: List[str], g2m_genes: List[str]) -> None:
    """
    Score cell cycle phases (equivalent to CellCycleScoring in Seurat).

    Args:
        adata: AnnData object
        s_genes: List of S phase genes
        g2m_genes: List of G2M phase genes
    """
    # Make sure genes are in the dataset
    s_genes_available = [g for g in s_genes if g in adata.var_names]
    g2m_genes_available = [g for g in g2m_genes if g in adata.var_names]

    if len(s_genes_available) == 0 or len(g2m_genes_available) == 0:
        warnings.warn(
            "Cell cycle genes not found in dataset. Skipping cell cycle scoring.")
        return

    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes_available,
        g2m_genes=g2m_genes_available
    )


# Cell cycle gene lists (from Seurat)
CELL_CYCLE_GENES = {
    's_genes': [
        'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2',
        'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2',
        'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7',
        'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1',
        'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B',
        'BRIP1', 'E2F8'
    ],
    'g2m_genes': [
        'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80',
        'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A',
        'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E',
        'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK',
        'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8',
        'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5',
        'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
    ]
}

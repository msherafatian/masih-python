"""
Marker gene detection utilities for MASIH application.

This module provides functions for differential expression analysis
to identify marker genes for clusters.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Optional, List, Dict, Any
from scipy import stats

from masih.config.settings import AppConfig


def find_all_markers(adata, cluster_key: str = 'leiden',
                     method: str = 'wilcoxon',
                     min_pct: float = 0.1,
                     logfc_threshold: float = 0.25,
                     only_pos: bool = True,
                     max_pval: float = 0.05) -> pd.DataFrame:
    """
    Find marker genes for all clusters (equivalent to FindAllMarkers in Seurat).

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments
        method: Statistical test ('wilcoxon', 't-test', 'logreg')
        min_pct: Minimum fraction of cells expressing the gene
        logfc_threshold: Minimum log fold change threshold
        only_pos: Only return positive markers
        max_pval: Maximum adjusted p-value

    Returns:
        DataFrame with marker genes for all clusters
    """
    import warnings
    warnings.filterwarnings(
        'ignore', message='invalid value encountered in log2')

    # IMPORTANT: Check if we have raw data stored
    # If data was normalized, we need to use the layer with raw counts
    use_raw = False
    if 'raw_counts' in adata.layers:
        # We stored raw counts in a layer during processing
        use_raw = False  # We'll manually handle this
        # Temporarily replace X with raw counts
        original_X = adata.X.copy()
        adata.X = adata.layers['raw_counts'].copy()
    elif adata.raw is not None:
        use_raw = True

    # Run Scanpy's rank_genes_groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        use_raw=use_raw,
        corr_method='benjamini-hochberg',
        pts=True,  # Calculate fraction of cells expressing
        tie_correct=True
    )

    # Restore original X if we modified it
    if 'raw_counts' in adata.layers:
        adata.X = original_X

    # Extract results into a DataFrame
    result = adata.uns['rank_genes_groups']
    clusters = result['names'].dtype.names

    markers_list = []

    for cluster in clusters:
        # Get all results for this cluster
        n_genes = len(result['names'][cluster])

        cluster_markers = pd.DataFrame({
            'gene': result['names'][cluster],
            'avg_log2FC': result['logfoldchanges'][cluster],
            'p_val': result['pvals'][cluster],
            'p_val_adj': result['pvals_adj'][cluster],
            'pct.1': result['pts'][cluster],
            'pct.2': result['pts_rest'][cluster],
            'cluster': cluster
        })

        # Remove NaN values
        cluster_markers = cluster_markers.dropna()

        # Filter based on criteria
        mask = (cluster_markers['p_val_adj'] < max_pval)

        if only_pos:
            mask &= (cluster_markers['avg_log2FC'] > logfc_threshold)
        else:
            mask &= (abs(cluster_markers['avg_log2FC']) > logfc_threshold)

        # Filter by minimum percentage
        mask &= (cluster_markers['pct.1'] > min_pct)

        cluster_markers = cluster_markers[mask]
        markers_list.append(cluster_markers)

    # Combine all clusters
    if len(markers_list) > 0:
        all_markers = pd.concat(markers_list, ignore_index=True)
    else:
        all_markers = pd.DataFrame()

    # Sort by cluster and log fold change
    if len(all_markers) > 0:
        all_markers = all_markers.sort_values(['cluster', 'avg_log2FC'],
                                              ascending=[True, False])

    return all_markers


def find_pairwise_markers(adata, cluster_key: str, ident_1: str, ident_2: str,
                          method: str = 'wilcoxon',
                          min_pct: float = 0.1,
                          logfc_threshold: float = 0.25,
                          only_pos: bool = True) -> pd.DataFrame:
    """
    Find marker genes between two specific clusters (equivalent to FindMarkers).

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments
        ident_1: First cluster identifier
        ident_2: Second cluster identifier
        method: Statistical test
        min_pct: Minimum fraction of cells expressing
        logfc_threshold: Minimum log fold change
        only_pos: Only positive markers

    Returns:
        DataFrame with marker genes
    """
    # Subset to the two clusters
    mask = adata.obs[cluster_key].isin([ident_1, ident_2])
    adata_subset = adata[mask].copy()

    # Run rank_genes_groups on subset
    sc.tl.rank_genes_groups(
        adata_subset,
        groupby=cluster_key,
        groups=[ident_1],
        reference=ident_2,
        method=method,
        use_raw=False,
        corr_method='benjamini-hochberg',
        pts=True
    )

    # Extract results
    result = adata_subset.uns['rank_genes_groups']

    markers = pd.DataFrame({
        'gene': result['names'][ident_1],
        'avg_log2FC': result['logfoldchanges'][ident_1],
        'p_val': result['pvals'][ident_1],
        'p_val_adj': result['pvals_adj'][ident_1],
        'pct.1': result['pts'][ident_1],
        'pct.2': result['pts_rest'][ident_1],
        'cluster': ident_1
    })

    # Filter based on criteria
    mask = (markers['pct.1'] > min_pct)

    if only_pos:
        mask &= (markers['avg_log2FC'] > logfc_threshold)
    else:
        mask &= (abs(markers['avg_log2FC']) > logfc_threshold)

    markers = markers[mask]
    markers = markers.sort_values('avg_log2FC', ascending=False)

    return markers


def get_top_markers(markers: pd.DataFrame, n: int = 10,
                    by_cluster: bool = True) -> pd.DataFrame:
    """
    Get top N marker genes per cluster.

    Args:
        markers: DataFrame of marker genes
        n: Number of top markers per cluster
        by_cluster: Whether to get top N per cluster or overall

    Returns:
        DataFrame with top markers
    """
    if by_cluster and 'cluster' in markers.columns:
        top_markers = markers.groupby('cluster').apply(
            lambda x: x.nlargest(n, 'avg_log2FC')
        ).reset_index(drop=True)
    else:
        top_markers = markers.nlargest(n, 'avg_log2FC')

    return top_markers


def filter_markers(markers: pd.DataFrame,
                   min_logfc: float = 0.25,
                   min_pct: float = 0.1,
                   max_pval: float = 0.05,
                   only_pos: bool = True) -> pd.DataFrame:
    """
    Filter markers by multiple criteria.

    Args:
        markers: DataFrame of markers
        min_logfc: Minimum log fold change
        min_pct: Minimum percent expression
        max_pval: Maximum adjusted p-value
        only_pos: Only positive markers

    Returns:
        Filtered DataFrame
    """
    filtered = markers.copy()

    if 'avg_log2FC' in filtered.columns:
        if only_pos:
            filtered = filtered[filtered['avg_log2FC'] > min_logfc]
        else:
            filtered = filtered[abs(filtered['avg_log2FC']) > min_logfc]

    if 'pct.1' in filtered.columns:
        filtered = filtered[filtered['pct.1'] > min_pct]

    if 'p_val_adj' in filtered.columns:
        filtered = filtered[filtered['p_val_adj'] < max_pval]

    return filtered


def summarize_markers(markers: pd.DataFrame) -> pd.DataFrame:
    """
    Create a summary of marker statistics per cluster.

    Args:
        markers: DataFrame of markers

    Returns:
        Summary DataFrame
    """
    if 'cluster' not in markers.columns:
        return pd.DataFrame()

    summary = markers.groupby('cluster').agg({
        'gene': 'count',
        'avg_log2FC': ['mean', 'median', 'max'],
        'pct.1': 'mean',
        'pct.2': 'mean'
    }).reset_index()

    # Flatten column names
    summary.columns = ['cluster', 'n_markers', 'avg_log2FC',
                       'median_log2FC', 'max_log2FC',
                       'avg_pct1', 'avg_pct2']

    return summary


def get_genes_for_enrichment(markers: pd.DataFrame, cluster: str,
                             n_genes: int = 100,
                             direction: str = 'up') -> List[str]:
    """
    Get marker genes for gene set enrichment analysis.

    Args:
        markers: DataFrame of markers
        cluster: Cluster to get genes for
        n_genes: Number of genes to return
        direction: Direction ('up', 'down', 'both')

    Returns:
        List of gene names
    """
    cluster_markers = markers[markers['cluster'] == cluster].copy()

    if direction == 'up':
        genes = cluster_markers[
            cluster_markers['avg_log2FC'] > 0
        ].nlargest(n_genes, 'avg_log2FC')['gene'].tolist()
    elif direction == 'down':
        genes = cluster_markers[
            cluster_markers['avg_log2FC'] < 0
        ].nsmallest(n_genes, 'avg_log2FC')['gene'].tolist()
    else:  # both
        cluster_markers['abs_log2FC'] = abs(cluster_markers['avg_log2FC'])
        genes = cluster_markers.nlargest(n_genes, 'abs_log2FC')[
            'gene'].tolist()

    return genes

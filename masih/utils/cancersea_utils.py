"""
CancerSEA pathway scoring utilities for MASIH application.

This module provides functions for scoring cancer-related functional states
using the CancerSEA gene sets.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import List, Dict, Optional

from masih.config.settings import AppConfig


def calculate_pathway_score(adata, pathway_name: str, gene_list: List[str],
                            ctrl_size: int = 50) -> str:
    """
    Calculate pathway score for a single CancerSEA pathway.

    Uses Scanpy's score_genes function which is equivalent to
    Seurat's AddModuleScore.

    Args:
        adata: AnnData object
        pathway_name: Name of the pathway
        gene_list: List of genes in the pathway
        ctrl_size: Number of control genes

    Returns:
        Name of the score column added to adata.obs
    """
    # Filter genes to those present in dataset
    available_genes = adata.var_names
    gene_list_filtered = [g for g in gene_list if g in available_genes]

    if len(gene_list_filtered) == 0:
        raise ValueError(
            f"No genes from {pathway_name} pathway found in dataset")

    # Score column name
    score_name = f"{pathway_name}_score"

    # Calculate score using Scanpy
    # This is equivalent to Seurat's AddModuleScore
    sc.tl.score_genes(
        adata,
        gene_list=gene_list_filtered,
        score_name=score_name,
        ctrl_size=ctrl_size,
        use_raw=False
    )

    return score_name


def calculate_multiple_pathways(adata, pathways: List[str],
                                cancersea_pathways: Dict[str, List[str]],
                                ctrl_size: int = 50) -> Dict[str, str]:
    """
    Calculate scores for multiple CancerSEA pathways.

    Args:
        adata: AnnData object
        pathways: List of pathway names to calculate
        cancersea_pathways: Dictionary mapping pathway names to gene lists
        ctrl_size: Number of control genes

    Returns:
        Dictionary mapping pathway names to score column names
    """
    scores_added = {}

    for pathway in pathways:
        if pathway not in cancersea_pathways:
            print(
                f"Warning: Pathway {pathway} not found in CancerSEA gene sets")
            continue

        try:
            gene_list = cancersea_pathways[pathway]
            score_name = calculate_pathway_score(
                adata,
                pathway,
                gene_list,
                ctrl_size
            )
            scores_added[pathway] = score_name
            print(f"Successfully calculated {pathway} score")

        except Exception as e:
            print(f"Error calculating {pathway}: {e}")

    return scores_added


def get_pathway_genes(pathway_name: str,
                      cancersea_pathways: Dict[str, List[str]]) -> List[str]:
    """
    Get gene list for a specific pathway.

    Args:
        pathway_name: Name of the pathway
        cancersea_pathways: Dictionary of pathway gene lists

    Returns:
        List of gene symbols
    """
    return cancersea_pathways.get(pathway_name, [])


def calculate_pathway_averages(adata, score_columns: Dict[str, str],
                               cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Calculate average pathway scores by cluster.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with average scores per cluster
    """
    clusters = sorted(adata.obs[cluster_key].unique())

    avg_scores = pd.DataFrame({
        'Cluster': clusters
    })

    for pathway, score_col in score_columns.items():
        if score_col not in adata.obs.columns:
            continue

        avg_by_cluster = adata.obs.groupby(cluster_key)[score_col].mean()
        avg_scores[pathway] = [avg_by_cluster.get(c, np.nan) for c in clusters]

    return avg_scores


def create_pathway_correlation(adata, score_columns: Dict[str, str]) -> pd.DataFrame:
    """
    Create correlation matrix of pathway scores.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names

    Returns:
        Correlation matrix as DataFrame
    """
    if len(score_columns) < 2:
        raise ValueError("Need at least 2 pathways for correlation")

    # Get score data
    score_data = pd.DataFrame()
    for pathway, score_col in score_columns.items():
        if score_col in adata.obs.columns:
            score_data[pathway] = adata.obs[score_col]

    # Calculate correlation
    return score_data.corr()


def get_pathway_statistics(adata, score_column: str,
                           cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Get statistics for a pathway score across clusters.

    Args:
        adata: AnnData object
        score_column: Name of the score column
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with statistics per cluster
    """
    if score_column not in adata.obs.columns:
        return pd.DataFrame()

    stats = adata.obs.groupby(cluster_key)[score_column].agg([
        'count', 'mean', 'std', 'min', 'max'
    ]).reset_index()

    stats.columns = ['Cluster', 'N_cells', 'Mean', 'Std', 'Min', 'Max']

    return stats

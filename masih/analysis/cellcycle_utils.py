"""
Cell cycle analysis utilities for MASIH application.

This module provides functions for cell cycle phase assignment
and scoring using marker genes.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import List, Optional


def score_cell_cycle(adata, s_genes: Optional[List[str]] = None,
                     g2m_genes: Optional[List[str]] = None) -> None:
    """
    Score cells for cell cycle phase using marker genes.

    This uses Scanpy's score_genes_cell_cycle which is equivalent
    to Seurat's CellCycleScoring.

    Args:
        adata: AnnData object
        s_genes: List of S phase marker genes
        g2m_genes: List of G2M phase marker genes
    """
    # Use default gene lists if not provided
    if s_genes is None:
        s_genes = get_s_phase_genes()

    if g2m_genes is None:
        g2m_genes = get_g2m_phase_genes()

    # Filter genes to those present in dataset
    s_genes_found = [g for g in s_genes if g in adata.var_names]
    g2m_genes_found = [g for g in g2m_genes if g in adata.var_names]

    print(f"Found {len(s_genes_found)}/{len(s_genes)} S phase genes")
    print(f"Found {len(g2m_genes_found)}/{len(g2m_genes)} G2M phase genes")

    # Score cell cycle
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes_found,
        g2m_genes=g2m_genes_found
    )

    # This adds 'S_score', 'G2M_score', and 'phase' to adata.obs


def get_s_phase_genes() -> List[str]:
    """
    Get default S phase marker genes.

    Returns:
        List of S phase gene symbols
    """
    return [
        'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2',
        'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2',
        'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7',
        'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1',
        'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B',
        'BRIP1', 'E2F8'
    ]


def get_g2m_phase_genes() -> List[str]:
    """
    Get default G2M phase marker genes.

    Returns:
        List of G2M phase gene symbols
    """
    return [
        'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80',
        'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A',
        'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E',
        'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK',
        'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8',
        'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5',
        'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
    ]


def calculate_phase_proportions(adata, cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Calculate cell cycle phase proportions per cluster.

    Args:
        adata: AnnData object with 'phase' column
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with phase proportions per cluster
    """
    if 'phase' not in adata.obs.columns:
        raise ValueError(
            "Cell cycle phases not scored. Run score_cell_cycle first.")

    # Calculate counts
    phase_counts = adata.obs.groupby(
        [cluster_key, 'phase']).size().unstack(fill_value=0)

    # Calculate proportions
    phase_props = phase_counts.div(phase_counts.sum(axis=1), axis=0)

    # Reset index for easier plotting
    phase_props = phase_props.reset_index()

    return phase_props


def get_cell_cycle_stats(adata, cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Get detailed cell cycle statistics per cluster.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with statistics
    """
    if 'phase' not in adata.obs.columns:
        return pd.DataFrame()

    clusters = sorted(adata.obs[cluster_key].unique())
    stats_list = []

    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        cluster_data = adata.obs[mask]

        stats = {
            'Cluster': cluster,
            'N_cells': mask.sum(),
            'G1_count': (cluster_data['phase'] == 'G1').sum(),
            'S_count': (cluster_data['phase'] == 'S').sum(),
            'G2M_count': (cluster_data['phase'] == 'G2M').sum(),
            'G1_percent': (cluster_data['phase'] == 'G1').sum() / len(cluster_data) * 100,
            'S_percent': (cluster_data['phase'] == 'S').sum() / len(cluster_data) * 100,
            'G2M_percent': (cluster_data['phase'] == 'G2M').sum() / len(cluster_data) * 100,
            'Mean_S_score': cluster_data['S_score'].mean(),
            'Mean_G2M_score': cluster_data['G2M_score'].mean()
        }

        stats_list.append(stats)

    return pd.DataFrame(stats_list)

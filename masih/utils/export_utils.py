"""
Export utilities for MASIH application.

This module provides functions for generating summary statistics,
methods text, and preparing data for export.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Any
from datetime import datetime


def calculate_cluster_stats(adata, cluster_key: str = 'leiden') -> pd.DataFrame:
    """
    Calculate statistics for each cluster.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with cluster statistics
    """
    clusters = sorted(adata.obs[cluster_key].unique())

    stats_list = []
    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        n_cells = mask.sum()
        pct_cells = (n_cells / len(adata)) * 100

        # Calculate average gene expression
        if hasattr(adata.X, 'toarray'):
            cluster_data = adata.X[mask].toarray()
        else:
            cluster_data = adata.X[mask]

        avg_genes = (cluster_data > 0).sum(axis=1).mean()
        avg_counts = cluster_data.sum(axis=1).mean()

        stats_list.append({
            'Cluster': cluster,
            'N_cells': n_cells,
            'Percent_cells': pct_cells,
            'Avg_genes_per_cell': avg_genes,
            'Avg_counts_per_cell': avg_counts
        })

    return pd.DataFrame(stats_list)


def calculate_cellcycle_proportions(adata, cluster_key: str = 'leiden',
                                    phase_key: str = 'phase') -> pd.DataFrame:
    """
    Calculate cell cycle phase proportions per cluster.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments
        phase_key: Column containing cell cycle phase

    Returns:
        DataFrame with proportions
    """
    if phase_key not in adata.obs.columns:
        return pd.DataFrame()

    clusters = sorted(adata.obs[cluster_key].unique())
    phases = adata.obs[phase_key].unique()

    props_list = []
    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        cluster_phases = adata.obs.loc[mask, phase_key]

        row = {'Cluster': cluster}
        for phase in phases:
            count = (cluster_phases == phase).sum()
            pct = (count / len(cluster_phases)) * 100
            row[f'{phase}_count'] = count
            row[f'{phase}_percent'] = pct

        props_list.append(row)

    return pd.DataFrame(props_list)


def summarize_trajectory(trajectory_data: Dict) -> pd.DataFrame:
    """
    Summarize trajectory analysis results.

    Args:
        trajectory_data: Dictionary containing trajectory results

    Returns:
        DataFrame with trajectory summary
    """
    if not trajectory_data:
        return pd.DataFrame()

    summary = {
        'Parameter': [],
        'Value': []
    }

    if 'method' in trajectory_data:
        summary['Parameter'].append('Method')
        summary['Value'].append(trajectory_data['method'])

    if 'n_components' in trajectory_data:
        summary['Parameter'].append('Number of components')
        summary['Value'].append(trajectory_data['n_components'])

    if 'root_cell' in trajectory_data:
        summary['Parameter'].append('Root cell')
        summary['Value'].append(trajectory_data['root_cell'])

    return pd.DataFrame(summary)


def generate_methods_text(processing_params: Dict,
                          cancersea_pathways: List[str],
                          has_markers: bool,
                          has_trajectory: bool) -> str:
    """
    Generate methods text for publication.

    Args:
        processing_params: Dictionary of processing parameters
        cancersea_pathways: List of calculated pathways
        has_markers: Whether marker analysis was performed
        has_trajectory: Whether trajectory analysis was performed

    Returns:
        Formatted methods text
    """
    text = f"""Single-cell RNA-seq Analysis Methods

Data Processing:
Single-cell RNA-seq data were processed using MASIH (Modular Analysis Suite for Interactive Heterogeneity, v{processing_params.get('version', '1.0.0')}). """

    # QC section
    if processing_params.get('qc_performed', False):
        text += f"""Quality control filtering was performed with the following thresholds: 
minimum {processing_params.get('min_genes', 200)} genes per cell, 
maximum {processing_params.get('max_genes', 2500)} genes per cell, 
and maximum {processing_params.get('max_mt_percent', 5)}% mitochondrial content. """

    # Normalization
    text += f"""Data were normalized using log-normalization with a scale factor of {processing_params.get('target_sum', 10000)}. """

    # Feature selection
    text += f"""Highly variable genes were identified (n={processing_params.get('n_top_genes', 2000)}) for downstream analysis. """

    # Dimensionality reduction
    text += f"""Principal component analysis (PCA) was performed, and the top {processing_params.get('n_pcs', 50)} principal components were used for """

    # Clustering
    text += f"""clustering. Graph-based clustering was performed using the Leiden algorithm with resolution {processing_params.get('resolution', 0.5)} and {processing_params.get('n_neighbors', 15)} neighbors. """

    # UMAP
    text += """Uniform Manifold Approximation and Projection (UMAP) was used for visualization. """

    # Marker genes
    if has_markers:
        text += """Differentially expressed marker genes for each cluster were identified using the Wilcoxon rank-sum test. """

    # CancerSEA
    if cancersea_pathways:
        pathways_str = ", ".join(cancersea_pathways)
        text += f"""Functional state scoring was performed using CancerSEA gene signatures for the following pathways: {pathways_str}. """

    # Trajectory
    if has_trajectory:
        text += """Trajectory inference was performed using PAGA (Partition-based graph abstraction) followed by diffusion pseudotime calculation. """

    text += """\n\nSoftware:
Analysis was conducted using Python 3.x with Scanpy for single-cell analysis, NumPy and Pandas for data manipulation, and Plotly for visualization."""

    return text


def generate_citations() -> str:
    """
    Generate citation text for all tools used.

    Returns:
        Formatted citation text
    """
    citations = """
Software Citations:

Python:
Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Scanpy:
Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology, 19(1), 15.

NumPy:
Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array programming with NumPy. Nature, 585(7825), 357-362.

Pandas:
McKinney, W. (2010). Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference (Vol. 445, pp. 51-56).

Plotly:
Plotly Technologies Inc. (2015). Collaborative data science. Montreal, QC: Plotly Technologies Inc.

Dash:
Plotly Technologies Inc. (2021). Dash: A Python Framework for Building Analytical Web Applications.

UMAP:
McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv preprint arXiv:1802.03426.

Leiden Algorithm:
Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reports, 9(1), 5233.

CancerSEA:
Yuan, H., Yan, M., Zhang, G., et al. (2019). CancerSEA: a cancer single-cell state atlas. Nucleic Acids Research, 47(D1), D900-D908.

PAGA:
Wolf, F. A., Hamey, F. K., Plass, M., et al. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biology, 20(1), 59.
"""
    return citations


def generate_figure_legend(adata, cluster_key: str = 'leiden') -> str:
    """
    Generate template figure legend.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        Formatted figure legend template
    """
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    n_clusters = len(adata.obs[cluster_key].unique())

    legend = f"""Figure Legend Template:

Single-cell RNA-seq analysis of [SAMPLE_NAME]. 
(A) UMAP visualization showing {n_clusters} distinct cell clusters identified from {n_cells:,} cells. 
(B) Expression of marker genes across clusters. 
(C) Functional state analysis using CancerSEA pathways. 
Each point represents a single cell, colored by [cluster assignment/pathway score/gene expression]. 
Data were filtered to retain {n_genes:,} genes and processed as described in Methods.
"""
    return legend


def print_analysis_summary(adata, processing_params: Dict,
                           cancersea_scores: Dict,
                           has_markers: bool,
                           has_trajectory: bool) -> str:
    """
    Generate comprehensive analysis summary.

    Args:
        adata: AnnData object
        processing_params: Processing parameters
        cancersea_scores: Dictionary of pathway scores
        has_markers: Whether markers were calculated
        has_trajectory: Whether trajectory was calculated

    Returns:
        Formatted summary text
    """
    cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

    summary = f"""
═══════════════════════════════════════════════════════════════
                    MASIH Analysis Summary
═══════════════════════════════════════════════════════════════

Dataset Information:
  • Total cells: {adata.n_obs:,}
  • Total genes: {adata.n_vars:,}
  • Number of clusters: {len(adata.obs[cluster_key].unique())}

Processing Parameters:
  • Normalization: Log-normalization (target sum: {processing_params.get('target_sum', 10000)})
  • Variable genes: {processing_params.get('n_top_genes', 2000)}
  • PCA components: {processing_params.get('n_pcs', 50)}
  • Clustering resolution: {processing_params.get('resolution', 0.5)}
  • Neighbors: {processing_params.get('n_neighbors', 15)}

Analyses Completed:
  ✓ Quality control and filtering
  ✓ Normalization and scaling
  ✓ Dimensionality reduction (PCA, UMAP)
  ✓ Clustering
"""

    if has_markers:
        summary += "  ✓ Marker gene identification\n"

    if cancersea_scores:
        summary += f"  ✓ CancerSEA pathway scoring ({len(cancersea_scores)} pathways)\n"

    if has_trajectory:
        summary += "  ✓ Trajectory inference\n"

    summary += f"""
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
═══════════════════════════════════════════════════════════════
"""

    return summary


def prepare_download_data(adata, include_options: List[str],
                          cluster_key: str = 'leiden',
                          cancersea_scores: Optional[Dict] = None,
                          markers_df: Optional[pd.DataFrame] = None,
                          trajectory_data: Optional[Dict] = None) -> Dict[str, pd.DataFrame]:
    """
    Prepare all requested data for download.

    Args:
        adata: AnnData object
        include_options: List of data types to include
        cluster_key: Column containing cluster assignments
        cancersea_scores: Dictionary of pathway scores
        markers_df: Marker genes DataFrame
        trajectory_data: Trajectory results dictionary

    Returns:
        Dictionary mapping sheet names to DataFrames
    """
    sheets = {}

    if 'metadata' in include_options:
        sheets['Cell_Metadata'] = pd.DataFrame(adata.obs)

    if 'cluster_stats' in include_options:
        sheets['Cluster_Statistics'] = calculate_cluster_stats(
            adata, cluster_key)

    if 'cancersea_scores' in include_options and cancersea_scores:
        score_data = pd.DataFrame()
        for pathway, score_col in cancersea_scores.items():
            if score_col in adata.obs.columns:
                score_data[pathway] = adata.obs[score_col]
        if not score_data.empty:
            sheets['CancerSEA_Scores'] = score_data

    if 'pathway_avg' in include_options and cancersea_scores:
        from masih.analysis.cancersea_utils import calculate_pathway_averages
        sheets['Pathway_Averages'] = calculate_pathway_averages(
            adata, cancersea_scores, cluster_key)

    if 'cellcycle_props' in include_options:
        cc_props = calculate_cellcycle_proportions(adata, cluster_key)
        if not cc_props.empty:
            sheets['Cell_Cycle_Proportions'] = cc_props

    if 'markers' in include_options and markers_df is not None:
        sheets['Marker_Genes'] = markers_df

    if 'trajectory' in include_options and trajectory_data:
        sheets['Trajectory_Summary'] = summarize_trajectory(trajectory_data)

    return sheets

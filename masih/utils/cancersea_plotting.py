"""
Visualization utilities for CancerSEA pathway analysis.

This module provides functions for creating heatmaps, violin plots,
and feature plots for pathway scores.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import List, Dict, Optional

from masih.utils.plotting import get_color_palette


def create_pathway_feature_plot(adata, score_column: str, pathway_name: str,
                                reduction: str = 'umap') -> go.Figure:
    """
    Create a feature plot showing pathway score on UMAP/tSNE.

    Args:
        adata: AnnData object
        score_column: Name of the score column in adata.obs
        pathway_name: Display name of the pathway
        reduction: Dimensionality reduction to use ('umap' or 'tsne')

    Returns:
        Plotly figure object
    """
    if score_column not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Score column '{score_column}' not found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    # Get coordinates
    if reduction == 'umap':
        coords = adata.obsm['X_umap']
        x_label, y_label = 'UMAP 1', 'UMAP 2'
    else:
        coords = adata.obsm['X_tsne']
        x_label, y_label = 'tSNE 1', 'tSNE 2'

    # Get scores
    scores = adata.obs[score_column].values

    # Create scatter plot
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=coords[:, 0],
        y=coords[:, 1],
        mode='markers',
        marker=dict(
            size=4,
            color=scores,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Score"),
            line=dict(width=0)
        ),
        text=[f"Score: {s:.3f}" for s in scores],
        hovertemplate='%{text}<extra></extra>',
        showlegend=False
    ))

    fig.update_layout(
        title=f"{pathway_name} Score",
        xaxis_title=x_label,
        yaxis_title=y_label,
        height=500,
        width=600,
        template='plotly_white'
    )

    return fig


def create_pathway_violin_plot(adata, score_column: str, pathway_name: str,
                               cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a violin plot showing pathway scores by cluster.

    Args:
        adata: AnnData object
        score_column: Name of the score column
        pathway_name: Display name of the pathway
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if score_column not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Score column '{score_column}' not found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    # Get data
    scores = adata.obs[score_column].values
    clusters = adata.obs[cluster_key].astype(str).values
    unique_clusters = sorted(adata.obs[cluster_key].unique())

    # Get colors
    colors = get_color_palette(len(unique_clusters), "discrete")
    color_map = dict(zip([str(c) for c in unique_clusters], colors))

    # Create violin plot
    fig = go.Figure()

    for cluster in unique_clusters:
        cluster_str = str(cluster)
        mask = clusters == cluster_str
        cluster_scores = scores[mask]

        fig.add_trace(go.Violin(
            y=cluster_scores,
            x=[cluster_str] * len(cluster_scores),
            name=cluster_str,
            marker=dict(color=color_map[cluster_str]),
            box_visible=True,
            meanline_visible=True,
            showlegend=False
        ))

    fig.update_layout(
        title=f"{pathway_name} Score by Cluster",
        xaxis_title="Cluster",
        yaxis_title="Score",
        height=500,
        width=600,
        template='plotly_white',
        violinmode='group'
    )

    return fig


def create_pathway_heatmap(adata, score_columns: Dict[str, str],
                           cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a heatmap of pathway scores across clusters.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if len(score_columns) < 2:
        fig = go.Figure()
        fig.add_annotation(
            text="Calculate at least 2 pathways to show heatmap",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        return fig

    # Calculate average scores per cluster
    clusters = sorted(adata.obs[cluster_key].unique())
    pathways = list(score_columns.keys())

    score_matrix = np.zeros((len(clusters), len(pathways)))

    for j, (pathway, score_col) in enumerate(score_columns.items()):
        if score_col not in adata.obs.columns:
            continue

        for i, cluster in enumerate(clusters):
            mask = adata.obs[cluster_key] == cluster
            score_matrix[i, j] = adata.obs.loc[mask, score_col].mean()

    # Z-score normalization
    from scipy.stats import zscore
    score_matrix_scaled = zscore(score_matrix, axis=1)

    # Handle NaN values
    score_matrix_scaled = np.nan_to_num(score_matrix_scaled, nan=0.0)

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=score_matrix_scaled,
        x=pathways,
        y=[f"Cluster {c}" for c in clusters],
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(title="Scaled<br>Score"),
        hovertemplate='Pathway: %{x}<br>Cluster: %{y}<br>Score: %{z:.2f}<extra></extra>'
    ))

    fig.update_layout(
        title="CancerSEA Pathway Scores by Cluster",
        xaxis_title="Pathway",
        yaxis_title="Cluster",
        height=max(400, len(clusters) * 40 + 100),
        width=max(600, len(pathways) * 60 + 100),
        xaxis=dict(tickangle=-45)
    )

    return fig


def create_pathway_comparison_plot(adata, score_col1: str, score_col2: str,
                                   pathway1: str, pathway2: str,
                                   cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a scatter plot comparing two pathway scores.

    Args:
        adata: AnnData object
        score_col1: First score column name
        score_col2: Second score column name
        pathway1: First pathway name
        pathway2: Second pathway name
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if score_col1 not in adata.obs.columns or score_col2 not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="One or both score columns not found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    scores1 = adata.obs[score_col1].values
    scores2 = adata.obs[score_col2].values
    clusters = adata.obs[cluster_key].astype(str).values
    unique_clusters = sorted(adata.obs[cluster_key].unique())

    # Get colors
    colors = get_color_palette(len(unique_clusters), "discrete")
    color_map = dict(zip([str(c) for c in unique_clusters], colors))

    # Create scatter plot
    fig = go.Figure()

    for cluster in unique_clusters:
        cluster_str = str(cluster)
        mask = clusters == cluster_str

        fig.add_trace(go.Scatter(
            x=scores1[mask],
            y=scores2[mask],
            mode='markers',
            name=f"Cluster {cluster_str}",
            marker=dict(
                size=5,
                color=color_map[cluster_str],
                opacity=0.6
            ),
            text=[f"Cluster: {cluster_str}" for _ in range(mask.sum())],
            hovertemplate='%{text}<br>%{xaxis.title.text}: %{x:.3f}<br>%{yaxis.title.text}: %{y:.3f}<extra></extra>'
        ))

    fig.update_layout(
        title=f"{pathway1} vs {pathway2}",
        xaxis_title=f"{pathway1} Score",
        yaxis_title=f"{pathway2} Score",
        height=500,
        width=600,
        template='plotly_white'
    )

    return fig

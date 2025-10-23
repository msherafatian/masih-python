"""
Comparative pathway visualization utilities.

This module provides functions for creating comparison plots
between multiple CancerSEA pathways.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Dict, Optional

from masih.utils.plotting import get_color_palette


def create_correlation_matrix_plot(correlation_matrix: pd.DataFrame,
                                   p_values: Optional[pd.DataFrame] = None) -> go.Figure:
    """
    Create interactive correlation matrix heatmap.

    Args:
        correlation_matrix: Correlation matrix DataFrame
        p_values: Optional p-value matrix

    Returns:
        Plotly figure object
    """
    if correlation_matrix.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No correlation data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        return fig

    # Create hover text with correlations and p-values
    hover_text = []
    for i, row in enumerate(correlation_matrix.index):
        hover_row = []
        for j, col in enumerate(correlation_matrix.columns):
            corr = correlation_matrix.iloc[i, j]
            text = f"{row} vs {col}<br>Correlation: {corr:.3f}"

            if p_values is not None and not pd.isna(p_values.iloc[i, j]):
                p_val = p_values.iloc[i, j]
                text += f"<br>P-value: {p_val:.3e}"

            hover_row.append(text)
        hover_text.append(hover_row)

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=correlation_matrix.values,
        x=correlation_matrix.columns,
        y=correlation_matrix.index,
        colorscale='RdBu_r',
        zmid=0,
        zmin=-1,
        zmax=1,
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        colorbar=dict(title="Correlation")
    ))

    # Add correlation coefficients as annotations
    annotations = []
    for i, row in enumerate(correlation_matrix.index):
        for j, col in enumerate(correlation_matrix.columns):
            if i != j:  # Don't annotate diagonal
                corr = correlation_matrix.iloc[i, j]
                annotations.append(
                    dict(
                        x=j,
                        y=i,
                        text=f"{corr:.2f}",
                        showarrow=False,
                        font=dict(size=10, color='black' if abs(
                            corr) < 0.5 else 'white')
                    )
                )

    fig.update_layout(
        title='Pathway Correlation Matrix',
        xaxis=dict(tickangle=-45),
        yaxis=dict(autorange='reversed'),
        height=500,
        width=600,
        annotations=annotations
    )

    return fig


def create_pathway_by_cluster_plot(pathway_data: pd.DataFrame) -> go.Figure:
    """
    Create grouped bar plot showing pathway scores by cluster.

    Args:
        pathway_data: DataFrame with pathway scores per cluster

    Returns:
        Plotly figure object
    """
    if pathway_data.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    print(f"DEBUG: Pathway data columns: {pathway_data.columns.tolist()}")

    # Extract pathway names (columns ending with '_mean')
    pathway_cols = [
        col for col in pathway_data.columns if col.endswith('_mean')]
    pathways = [col.replace('_mean', '') for col in pathway_cols]

    print(f"DEBUG: Found pathways: {pathways}")

    if not pathways:
        fig = go.Figure()
        fig.add_annotation(
            text="No pathway data found",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    fig = go.Figure()

    # Add bar for each pathway
    for pathway in pathways:
        mean_col = f"{pathway}_mean"
        se_col = f"{pathway}_se"

        if mean_col in pathway_data.columns:
            means = pathway_data[mean_col].values

            # Get error bars if available
            if se_col in pathway_data.columns:
                errors = pathway_data[se_col].values
                error_array = errors
            else:
                error_array = None

            fig.add_trace(go.Bar(
                name=pathway,
                x=pathway_data['Cluster'].astype(str),
                y=means,
                error_y=dict(
                    type='data', array=error_array) if error_array is not None else None,
                hovertemplate=f'{pathway}<br>Cluster: %{{x}}<br>Mean: %{{y:.3f}}<extra></extra>'
            ))

    fig.update_layout(
        title='Pathway Scores by Cluster',
        xaxis_title='Cluster',
        yaxis_title='Mean Score',
        barmode='group',  # This creates grouped bars, not stacked
        height=500,
        template='plotly_white',
        xaxis=dict(tickangle=-45),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        )
    )

    return fig


def create_comparison_heatmap(adata, score_columns: Dict[str, str],
                              cluster_key: str = 'leiden') -> go.Figure:
    """
    Create heatmap comparing pathway scores across clusters.

    Args:
        adata: AnnData object
        score_columns: Dictionary mapping pathway names to score column names
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if not score_columns or len(score_columns) < 2:
        fig = go.Figure()
        fig.add_annotation(
            text="Need at least 2 pathways for comparison",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    print(f"DEBUG: Creating heatmap with {len(score_columns)} pathways")

    # Calculate mean scores per cluster
    clusters = sorted(adata.obs[cluster_key].unique())
    pathways = list(score_columns.keys())

    score_matrix = np.zeros((len(clusters), len(pathways)))

    for i, cluster in enumerate(clusters):
        mask = adata.obs[cluster_key] == cluster
        for j, (pathway, score_col) in enumerate(score_columns.items()):
            if score_col in adata.obs.columns:
                score_matrix[i, j] = adata.obs.loc[mask, score_col].mean()
            else:
                print(f"Warning: Score column '{score_col}' not found")

    print(f"DEBUG: Score matrix shape: {score_matrix.shape}")
    print(f"DEBUG: Score matrix:\n{score_matrix}")

    # Z-score normalization across clusters (per pathway)
    from scipy.stats import zscore
    score_matrix_scaled = zscore(score_matrix, axis=0)
    score_matrix_scaled = np.nan_to_num(score_matrix_scaled, nan=0.0)

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=score_matrix_scaled,
        x=pathways,
        y=[f"Cluster {c}" for c in clusters],
        colorscale='Viridis',
        colorbar=dict(title=dict(text="Scaled Score", side='right')),
        hovertemplate='Pathway: %{x}<br>Cluster: %{y}<br>Scaled Score: %{z:.2f}<extra></extra>'
    ))

    # Calculate responsive height and width
    n_clusters = len(clusters)
    n_pathways = len(pathways)

    # Height: 50px per cluster + 150px for margins
    height = max(400, n_clusters * 50 + 150)

    # Width: 100px per pathway + 200px for margins
    width = max(600, n_pathways * 100 + 200)

    fig.update_layout(
        title='Pathway Comparison Heatmap (Z-scored)',
        xaxis_title='Pathway',
        yaxis_title='',
        height=height,
        width=width,
        xaxis=dict(tickangle=-45, side='bottom'),
        margin=dict(l=100, r=100, t=100, b=100)
    )

    return fig


def create_pathway_profile_plot(profile_data: pd.DataFrame) -> go.Figure:
    """
    Create line plot showing pathway profiles across clusters.

    Args:
        profile_data: DataFrame with pathway profiles

    Returns:
        Plotly figure object
    """
    if profile_data.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No profile data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    pathways = profile_data['Pathway'].unique()

    fig = go.Figure()

    for pathway in pathways:
        pathway_data = profile_data[profile_data['Pathway'] == pathway]

        fig.add_trace(go.Scatter(
            x=pathway_data['Cluster'].astype(str),
            y=pathway_data['Mean'],
            mode='lines+markers',
            name=pathway,
            error_y=dict(
                type='data',
                array=pathway_data['SE']
            ) if 'SE' in pathway_data.columns else None
        ))

    fig.update_layout(
        title='Pathway Activity Profiles',
        xaxis_title='Cluster',
        yaxis_title='Mean Score',
        height=500,
        width=800,
        template='plotly_white',
        hovermode='x unified'
    )

    return fig

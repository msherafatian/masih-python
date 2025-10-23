"""
Cell cycle visualization utilities.

This module provides functions for creating cell cycle-related plots.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import Optional

from masih.utils.plotting import get_color_palette


def create_cellcycle_umap(adata, reduction: str = 'umap') -> go.Figure:
    """
    Create UMAP plot colored by cell cycle phase.

    Args:
        adata: AnnData object with 'phase' column
        reduction: Dimensionality reduction to use

    Returns:
        Plotly figure object
    """
    if 'phase' not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="Cell cycle phases not scored",
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

    # Get phases
    phases = adata.obs['phase'].astype(str)
    unique_phases = ['G1', 'S', 'G2M']  # Fixed order

    # Define colors (matching R version)
    phase_colors = {
        'G1': '#E41A1C',   # Red
        'S': '#4DAF4A',    # Green
        'G2M': '#377EB8'   # Blue
    }

    # Create plot
    fig = go.Figure()

    for phase in unique_phases:
        if phase in phases.values:
            mask = phases == phase
            fig.add_trace(go.Scatter(
                x=coords[mask, 0],
                y=coords[mask, 1],
                mode='markers',
                name=phase,
                marker=dict(
                    size=4,
                    color=phase_colors[phase],
                    line=dict(width=0)
                ),
                text=[f"Phase: {phase}"] * mask.sum(),
                hovertemplate='%{text}<extra></extra>'
            ))

    fig.update_layout(
        title='Cell Cycle Phases',
        xaxis_title=x_label,
        yaxis_title=y_label,
        height=500,
        width=600,
        template='plotly_white'
    )

    return fig


def create_cycle_score_plot(adata, score_type: str = 'S_score',
                            reduction: str = 'umap') -> go.Figure:
    """
    Create UMAP plot colored by cell cycle score.

    Args:
        adata: AnnData object
        score_type: Score to plot ('S_score' or 'G2M_score')
        reduction: Dimensionality reduction to use

    Returns:
        Plotly figure object
    """
    if score_type not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text=f"{score_type} not found",
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
    scores = adata.obs[score_type].values

    # Create plot
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

    title = 'S Phase Score' if score_type == 'S_score' else 'G2M Phase Score'

    fig.update_layout(
        title=title,
        xaxis_title=x_label,
        yaxis_title=y_label,
        height=450,
        width=600,
        template='plotly_white'
    )

    return fig


def create_phase_by_cluster_plot(adata, cluster_key: str = 'leiden') -> go.Figure:
    """
    Create stacked bar plot of phase distribution by cluster.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if 'phase' not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="Cell cycle phases not scored",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    from masih.analysis.cellcycle_utils import calculate_phase_proportions

    # Get proportions
    phase_props = calculate_phase_proportions(adata, cluster_key)

    # Define colors
    phase_colors = {
        'G1': '#E41A1C',
        'S': '#4DAF4A',
        'G2M': '#377EB8'
    }

    # Create stacked bar chart
    fig = go.Figure()

    for phase in ['G1', 'S', 'G2M']:
        if phase in phase_props.columns:
            fig.add_trace(go.Bar(
                name=phase,
                x=phase_props[cluster_key].astype(str),
                y=phase_props[phase],
                marker_color=phase_colors[phase]
            ))

    fig.update_layout(
        barmode='stack',
        title='Cell Cycle Phase Distribution by Cluster',
        xaxis_title='Cluster',
        yaxis_title='Proportion',
        height=400,
        template='plotly_white',
        xaxis=dict(tickangle=-45)
    )

    return fig

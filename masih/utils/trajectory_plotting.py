"""
Trajectory visualization utilities.

This module provides functions for creating trajectory-related plots
using PAGA and diffusion pseudotime.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import Optional, List
import scanpy as sc
import matplotlib.pyplot as plt
import io
import base64

from masih.utils.plotting import get_color_palette


def create_paga_trajectory_plot(adata, cluster_key: str = 'leiden') -> go.Figure:
    """
    Create PAGA trajectory plot showing cluster connectivity.

    Args:
        adata: AnnData object with PAGA results
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if 'paga' not in adata.uns:
        fig = go.Figure()
        fig.add_annotation(
            text="PAGA not calculated",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    # Get PAGA positions and connectivity
    pos = adata.uns['paga']['pos']
    connectivities = adata.uns['paga']['connectivities']
    clusters = adata.obs[cluster_key].cat.categories

    # Get colors for clusters
    n_clusters = len(clusters)
    colors = get_color_palette(n_clusters, "discrete")
    color_map = dict(zip(clusters, colors))

    fig = go.Figure()

    # Draw edges (connections)
    threshold = np.percentile(connectivities.data, 50) if len(
        connectivities.data) > 0 else 0

    for i in range(len(clusters)):
        for j in range(i + 1, len(clusters)):
            weight = connectivities[i, j]
            if weight > threshold:
                fig.add_trace(go.Scatter(
                    x=[pos[i, 0], pos[j, 0]],
                    y=[pos[i, 1], pos[j, 1]],
                    mode='lines',
                    line=dict(
                        width=weight * 10,  # Scale width by connectivity
                        color='rgba(150, 150, 150, 0.5)'
                    ),
                    hoverinfo='skip',
                    showlegend=False
                ))

    # Draw nodes (clusters)
    for i, cluster in enumerate(clusters):
        # Count cells in cluster
        n_cells = (adata.obs[cluster_key] == cluster).sum()

        fig.add_trace(go.Scatter(
            x=[pos[i, 0]],
            y=[pos[i, 1]],
            mode='markers+text',
            marker=dict(
                size=20 + np.sqrt(n_cells),
                color=color_map[cluster],
                line=dict(width=2, color='white')
            ),
            text=str(cluster),
            textposition='middle center',
            textfont=dict(size=12, color='white', family='Arial Black'),
            name=str(cluster),
            hovertemplate=f'Cluster: {cluster}<br>Cells: {n_cells}<extra></extra>'
        ))

    fig.update_layout(
        title='Trajectory Analysis (PAGA)',
        xaxis=dict(showgrid=False, zeroline=False,
                   showticklabels=False, title=''),
        yaxis=dict(showgrid=False, zeroline=False,
                   showticklabels=False, title=''),
        height=600,
        width=800,
        showlegend=True,
        hovermode='closest',
        plot_bgcolor='white'
    )

    return fig


def create_pseudotime_plot(adata, reduction: str = 'umap') -> go.Figure:
    """
    Create UMAP/PCA plot colored by pseudotime.

    Args:
        adata: AnnData object with pseudotime
        reduction: Dimensionality reduction to use

    Returns:
        Plotly figure object
    """
    if 'dpt_pseudotime' not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="Pseudotime not calculated",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    # Get coordinates
    if reduction == 'umap' and 'X_umap' in adata.obsm:
        coords = adata.obsm['X_umap']
        x_label, y_label = 'UMAP 1', 'UMAP 2'
    elif 'X_pca' in adata.obsm:
        coords = adata.obsm['X_pca'][:, :2]
        x_label, y_label = 'PC 1', 'PC 2'
    else:
        fig = go.Figure()
        fig.add_annotation(
            text="No coordinates available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Get pseudotime
    pseudotime = adata.obs['dpt_pseudotime'].values

    # Create plot
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=coords[:, 0],
        y=coords[:, 1],
        mode='markers',
        marker=dict(
            size=4,
            color=pseudotime,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title=dict(text="Pseudotime",
                                     side='right')),
            line=dict(width=0)
        ),
        text=[f"Pseudotime: {pt:.3f}" if not np.isnan(pt) else "NA"
              for pt in pseudotime],
        hovertemplate='%{text}<extra></extra>',
        showlegend=False
    ))

    fig.update_layout(
        title='Pseudotime',
        xaxis_title=x_label,
        yaxis_title=y_label,
        height=600,
        width=800,
        template='plotly_white'
    )

    return fig


def create_stratified_trajectory_plot(adata, cluster_key: str = 'leiden',
                                      neotype_col: Optional[str] = None) -> go.Figure:
    """
    Create stratified trajectory plot (clusters with trajectory overlay).

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments
        neotype_col: Optional column for additional grouping

    Returns:
        Plotly figure object
    """
    if 'paga' not in adata.uns:
        fig = go.Figure()
        fig.add_annotation(
            text="PAGA not calculated",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Get coordinates
    if 'X_pca' in adata.obsm:
        coords = adata.obsm['X_pca'][:, :2]
        x_label, y_label = 'PC 1', 'PC 2'
    elif 'X_umap' in adata.obsm:
        coords = adata.obsm['X_umap']
        x_label, y_label = 'UMAP 1', 'UMAP 2'
    else:
        fig = go.Figure()
        fig.add_annotation(
            text="No coordinates available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Get clusters
    clusters = adata.obs[cluster_key].astype(str)
    unique_clusters = sorted(adata.obs[cluster_key].unique())

    # Get colors
    colors = get_color_palette(len(unique_clusters), "discrete")
    color_map = dict(zip([str(c) for c in unique_clusters], colors))

    fig = go.Figure()

    # Plot cells colored by cluster
    for cluster in unique_clusters:
        mask = clusters == str(cluster)

        # Determine marker shape based on neotype if provided
        if neotype_col and neotype_col in adata.obs.columns:
            neotypes = adata.obs.loc[mask, neotype_col].astype(str)
            unique_neotypes = neotypes.unique()

            for neotype in unique_neotypes:
                neotype_mask = mask & (
                    adata.obs[neotype_col].astype(str) == neotype)

                fig.add_trace(go.Scatter(
                    x=coords[neotype_mask, 0],
                    y=coords[neotype_mask, 1],
                    mode='markers',
                    name=f"{cluster} ({neotype})",
                    marker=dict(
                        size=6,
                        color=color_map[str(cluster)],
                        symbol='circle' if neotype == unique_neotypes[0] else 'square',
                        line=dict(width=0.5, color='white')
                    ),
                    text=[
                        f"Cluster: {cluster}<br>Type: {neotype}"] * neotype_mask.sum(),
                    hovertemplate='%{text}<extra></extra>'
                ))
        else:
            fig.add_trace(go.Scatter(
                x=coords[mask, 0],
                y=coords[mask, 1],
                mode='markers',
                name=str(cluster),
                marker=dict(
                    size=6,
                    color=color_map[str(cluster)],
                    line=dict(width=0.5, color='white')
                ),
                text=[f"Cluster: {cluster}"] * mask.sum(),
                hovertemplate='%{text}<extra></extra>'
            ))

    # Overlay PAGA trajectory
    pos = adata.uns['paga']['pos']
    connectivities = adata.uns['paga']['connectivities']
    threshold = np.percentile(connectivities.data, 50) if len(
        connectivities.data) > 0 else 0

    for i in range(len(unique_clusters)):
        for j in range(i + 1, len(unique_clusters)):
            weight = connectivities[i, j]
            if weight > threshold:
                fig.add_trace(go.Scatter(
                    x=[pos[i, 0], pos[j, 0]],
                    y=[pos[i, 1], pos[j, 1]],
                    mode='lines',
                    line=dict(width=2, color='black'),
                    hoverinfo='skip',
                    showlegend=False
                ))

    fig.update_layout(
        title='Stratified Trajectory',
        xaxis_title=x_label,
        yaxis_title=y_label,
        height=600,
        width=800,
        template='plotly_white'
    )

    return fig


def create_trajectory_matplotlib_plot(adata, cluster_key: str = 'leiden') -> str:
    """
    Create trajectory plot using matplotlib (for export compatibility).
    Returns base64-encoded image.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        Base64-encoded PNG image
    """
    if 'paga' not in adata.uns:
        return None

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Use scanpy's built-in PAGA plot
    sc.pl.paga(
        adata,
        color=cluster_key,
        ax=ax,
        show=False,
        frameon=False,
        node_size_scale=2,
        edge_width_scale=0.5
    )

    ax.set_title('Trajectory Analysis', fontsize=14, fontweight='bold')

    # Convert to base64
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)

    return img_base64


def create_pseudotime_matplotlib_plot(adata, reduction: str = 'umap') -> str:
    """
    Create pseudotime plot using matplotlib (for export).

    Args:
        adata: AnnData object
        reduction: Dimensionality reduction

    Returns:
        Base64-encoded PNG image
    """
    if 'dpt_pseudotime' not in adata.obs.columns:
        return None

    fig, ax = plt.subplots(figsize=(10, 8))

    # Get coordinates
    if reduction == 'umap' and 'X_umap' in adata.obsm:
        coords = adata.obsm['X_umap']
    else:
        coords = adata.obsm['X_pca'][:, :2]

    pseudotime = adata.obs['dpt_pseudotime'].values

    # Create scatter plot
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=pseudotime,
        cmap='viridis',
        s=10,
        edgecolors='none'
    )

    plt.colorbar(scatter, ax=ax, label='Pseudotime')
    ax.set_xlabel(f'{reduction.upper()} 1' if reduction == 'umap' else 'PC 1')
    ax.set_ylabel(f'{reduction.upper()} 2' if reduction == 'umap' else 'PC 2')
    ax.set_title('Pseudotime', fontsize=14, fontweight='bold')

    # Convert to base64
    buf = io.BytesIO()
    plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)

    return img_base64

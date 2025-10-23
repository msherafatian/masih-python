"""
Plot export utilities for MASIH application.

This module provides functions for creating exportable plots
with customization options.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from typing import Optional, Dict, Any
import io
import base64

# Import only what exists
from masih.utils.plotting import (
    create_feature_plot,
    get_color_palette
)
from masih.utils.marker_plotting import create_marker_heatmap
from masih.utils.cancersea_plotting import (
    create_pathway_feature_plot,
    create_pathway_heatmap
)

"""
Plot export utilities for MASIH application.

This module provides functions for creating exportable plots
with customization options.
"""


# Import only what exists


# Add the create_umap_plot function here since it doesn't exist yet
def create_umap_plot(adata, color_by: str = 'leiden',
                     reduction: str = 'umap') -> go.Figure:
    """
    Create a UMAP plot colored by a variable.

    Args:
        adata: AnnData object
        color_by: Column to color by
        reduction: Reduction to use ('umap' or 'tsne')

    Returns:
        Plotly figure object
    """
    # Get coordinates
    if reduction == 'umap':
        coords = adata.obsm['X_umap']
        x_label, y_label = 'UMAP 1', 'UMAP 2'
    else:
        coords = adata.obsm['X_tsne']
        x_label, y_label = 'tSNE 1', 'tSNE 2'

    # Get colors
    if color_by in adata.obs.columns:
        color_data = adata.obs[color_by].astype(str)

        # Check if categorical or continuous
        if adata.obs[color_by].dtype.name == 'category' or len(color_data.unique()) < 20:
            # Categorical
            unique_vals = sorted(color_data.unique())
            colors = get_color_palette(len(unique_vals), "discrete")
            color_map = dict(zip(unique_vals, colors))

            fig = go.Figure()
            for val in unique_vals:
                mask = color_data == val
                fig.add_trace(go.Scatter(
                    x=coords[mask, 0],
                    y=coords[mask, 1],
                    mode='markers',
                    name=str(val),
                    marker=dict(size=4, color=color_map[val]),
                    showlegend=True
                ))
        else:
            # Continuous
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=coords[:, 0],
                y=coords[:, 1],
                mode='markers',
                marker=dict(
                    size=4,
                    color=adata.obs[color_by],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title=color_by)
                ),
                showlegend=False
            ))
    else:
        # Default: all same color
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=coords[:, 0],
            y=coords[:, 1],
            mode='markers',
            marker=dict(size=4, color='blue'),
            showlegend=False
        ))

    fig.update_layout(
        xaxis_title=x_label,
        yaxis_title=y_label,
        template='plotly_white',
        width=600,
        height=500
    )

    return fig


def create_export_plot(adata, plot_type: str,
                       cluster_key: str = 'leiden',
                       gene: Optional[str] = None,
                       pathway: Optional[str] = None,
                       cancersea_scores: Optional[Dict] = None,
                       markers_genes: Optional[list] = None,
                       custom_title: Optional[str] = None,
                       white_bg: bool = True) -> go.Figure:
    """
    Create a plot configured for export.

    Args:
        adata: AnnData object
        plot_type: Type of plot to create
        cluster_key: Column containing cluster assignments
        gene: Gene name (for gene expression plots)
        pathway: Pathway name (for pathway plots)
        cancersea_scores: Dictionary of pathway scores
        markers_genes: List of marker genes for heatmap
        custom_title: Custom title for the plot
        white_bg: Use white background

    Returns:
        Plotly figure object
    """
    fig = None

    # Cluster plot
    if plot_type == 'cluster':
        fig = create_umap_plot(adata, color_by=cluster_key)
        if custom_title:
            fig.update_layout(title=custom_title)

    # Cell cycle plot
    elif plot_type == 'cellcycle' and 'phase' in adata.obs.columns:
        fig = create_umap_plot(adata, color_by='phase')
        if custom_title:
            fig.update_layout(title=custom_title)

    # Cell cycle distribution
    elif plot_type == 'cellcycle_distribution' and 'phase' in adata.obs.columns:
        fig = create_cellcycle_distribution_plot(adata, cluster_key)
        if custom_title:
            fig.update_layout(title=custom_title)

    # Gene expression
    elif plot_type == 'gene' and gene:
        fig = create_feature_plot(
            adata, gene, title=custom_title or f"{gene} Expression")

    # CancerSEA feature plot
    elif plot_type == 'cancersea_feature' and pathway and cancersea_scores:
        if pathway in cancersea_scores:
            score_col = cancersea_scores[pathway]
            fig = create_pathway_feature_plot(adata, score_col, pathway)
            if custom_title:
                fig.update_layout(title=custom_title)

    # CancerSEA heatmap
    elif plot_type == 'cancersea_heatmap' and cancersea_scores:
        fig = create_pathway_heatmap(adata, cancersea_scores, cluster_key)
        if custom_title:
            fig.update_layout(title=custom_title)

    # Pathway correlation
    elif plot_type == 'correlation_matrix' and cancersea_scores:
        from masih.utils.cancersea_utils import create_pathway_correlation
        corr_matrix = create_pathway_correlation(adata, cancersea_scores)
        fig = create_correlation_heatmap(corr_matrix)
        if custom_title:
            fig.update_layout(title=custom_title)

    # Marker heatmap
    elif plot_type == 'marker_heatmap' and markers_genes:
        fig = create_marker_heatmap(adata, markers_genes, cluster_key)
        if custom_title:
            fig.update_layout(title=custom_title)

    # Cluster tree
    elif plot_type == 'tree':
        fig = create_cluster_tree_plot(adata, cluster_key)
        if custom_title:
            fig.update_layout(title=custom_title)

    # Apply white background if requested
    if fig and white_bg:
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(color='black')
        )

    return fig


def create_cellcycle_distribution_plot(adata, cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a bar plot showing cell cycle distribution per cluster.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if 'phase' not in adata.obs.columns:
        fig = go.Figure()
        fig.add_annotation(
            text="Cell cycle phase not available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    # Calculate proportions
    clusters = sorted(adata.obs[cluster_key].unique())
    phases = sorted(adata.obs['phase'].unique())

    data = []
    for phase in phases:
        proportions = []
        for cluster in clusters:
            mask = (adata.obs[cluster_key] == cluster) & (
                adata.obs['phase'] == phase)
            cluster_total = (adata.obs[cluster_key] == cluster).sum()
            proportion = (mask.sum() / cluster_total) * \
                100 if cluster_total > 0 else 0
            proportions.append(proportion)

        data.append(go.Bar(
            name=phase,
            x=[str(c) for c in clusters],
            y=proportions
        ))

    fig = go.Figure(data=data)
    fig.update_layout(
        barmode='stack',
        title='Cell Cycle Distribution by Cluster',
        xaxis_title='Cluster',
        yaxis_title='Percentage',
        template='plotly_white'
    )

    return fig


def create_correlation_heatmap(corr_matrix: pd.DataFrame) -> go.Figure:
    """
    Create a correlation heatmap.

    Args:
        corr_matrix: Correlation matrix DataFrame

    Returns:
        Plotly figure object
    """
    fig = go.Figure(data=go.Heatmap(
        z=corr_matrix.values,
        x=corr_matrix.columns,
        y=corr_matrix.index,
        colorscale='RdBu_r',
        zmid=0,
        zmin=-1,
        zmax=1,
        colorbar=dict(title="Correlation"),
        hovertemplate='%{x} vs %{y}<br>Correlation: %{z:.3f}<extra></extra>'
    ))

    fig.update_layout(
        title='Pathway Correlation Matrix',
        xaxis_title='',
        yaxis_title='',
        width=600,
        height=600,
        xaxis=dict(tickangle=-45)
    )

    return fig


def create_cluster_tree_plot(adata, cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a cluster dendrogram/tree plot.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    # Calculate cluster similarity based on average gene expression
    clusters = sorted(adata.obs[cluster_key].unique())
    n_clusters = len(clusters)

    # Calculate average expression per cluster
    cluster_means = []
    for cluster in clusters:
        mask = adata.obs[cluster_key] == cluster
        if hasattr(adata.X, 'toarray'):
            cluster_data = adata.X[mask].toarray()
        else:
            cluster_data = adata.X[mask]
        cluster_means.append(cluster_data.mean(axis=0))

    cluster_means = np.array(cluster_means)

    # Calculate correlation between clusters
    from scipy.stats import spearmanr
    from scipy.cluster.hierarchy import dendrogram, linkage

    # Perform hierarchical clustering
    linkage_matrix = linkage(
        cluster_means, method='average', metric='correlation')

    # Create dendrogram
    dendro = dendrogram(linkage_matrix, labels=[
                        str(c) for c in clusters], no_plot=True)

    # Extract dendrogram data for plotting
    icoord = np.array(dendro['icoord'])
    dcoord = np.array(dendro['dcoord'])

    # Create plot
    fig = go.Figure()

    for i in range(len(icoord)):
        fig.add_trace(go.Scatter(
            x=icoord[i],
            y=dcoord[i],
            mode='lines',
            line=dict(color='black', width=2),
            hoverinfo='none',
            showlegend=False
        ))

    # Add cluster labels
    labels_pos = {label: pos for pos, label in zip(
        dendro['ivl'], range(len(dendro['ivl'])))}

    fig.update_layout(
        title='Cluster Hierarchy Tree',
        xaxis=dict(
            tickvals=list(range(5, len(clusters) * 10 + 5, 10)),
            ticktext=dendro['ivl'],
            title='Cluster'
        ),
        yaxis=dict(title='Distance'),
        showlegend=False,
        template='plotly_white',
        height=500,
        width=800
    )

    return fig


def fig_to_base64(fig: go.Figure, format: str = 'png',
                  width: int = 800, height: int = 600,
                  scale: float = 2) -> str:
    """
    Convert Plotly figure to base64-encoded image.

    Args:
        fig: Plotly figure
        format: Image format ('png', 'jpeg', 'svg', 'pdf')
        width: Width in pixels
        height: Height in pixels
        scale: Scale factor for resolution

    Returns:
        Base64-encoded image string
    """
    img_bytes = fig.to_image(
        format=format,
        width=width,
        height=height,
        scale=scale
    )

    img_base64 = base64.b64encode(img_bytes).decode()
    return img_base64


def get_available_plot_types(adata, cancersea_scores: Optional[Dict] = None,
                             has_markers: bool = False,
                             has_trajectory: bool = False) -> list:
    """
    Get list of available plot types based on current analysis state.

    Args:
        adata: AnnData object
        cancersea_scores: Dictionary of pathway scores
        has_markers: Whether markers have been calculated
        has_trajectory: Whether trajectory has been calculated

    Returns:
        List of available plot type identifiers
    """
    plot_types = ['cluster']

    if 'phase' in adata.obs.columns:
        plot_types.extend(['cellcycle', 'cellcycle_distribution'])

    if cancersea_scores and len(cancersea_scores) > 0:
        plot_types.extend(['cancersea_feature'])

        if len(cancersea_scores) >= 2:
            plot_types.extend(['cancersea_heatmap', 'correlation_matrix'])

    if has_markers:
        plot_types.append('marker_heatmap')

    if has_trajectory:
        plot_types.extend(['trajectory_main', 'pseudotime'])

    plot_types.append('tree')

    return plot_types

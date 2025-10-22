"""
Plotting utility functions for MASIH application.

This module provides consistent plotting styles and helper functions
for creating publication-quality visualizations.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, List, Dict, Any

from masih.config.settings import AppConfig


def get_color_palette(n: int, palette: str = "discrete") -> List[str]:
    """
    Get a consistent color palette.

    Args:
        n: Number of colors needed
        palette: Type of palette ("discrete", "continuous", "diverging")

    Returns:
        List of color strings
    """
    if palette == "discrete":
        if n <= 10:
            # Use Plotly's default discrete colors
            colors = px.colors.qualitative.Plotly
        elif n <= 24:
            colors = px.colors.qualitative.Light24
        else:
            # Generate colors using a colorscale
            colors = px.colors.sample_colorscale(
                "viridis", [i/(n-1) for i in range(n)])
        return colors[:n]

    elif palette == "continuous":
        return px.colors.sample_colorscale("viridis", [i/(n-1) for i in range(n)])

    elif palette == "diverging":
        return px.colors.sample_colorscale("RdBu_r", [i/(n-1) for i in range(n)])

    return px.colors.qualitative.Plotly[:n]


def create_dim_plot(adata, reduction: str = "umap", color_by: str = "leiden",
                    title: str = None, show_legend: bool = True,
                    point_size: int = 3) -> go.Figure:
    """
    Create a dimensionality reduction plot (equivalent to DimPlot in Seurat).

    Args:
        adata: AnnData object
        reduction: Reduction to plot ("umap", "tsne", "pca")
        color_by: Column to color by
        title: Plot title
        show_legend: Whether to show legend
        point_size: Size of points

    Returns:
        Plotly figure object
    """
    # Get coordinates
    reduction_key = f"X_{reduction}"
    if reduction_key not in adata.obsm:
        raise ValueError(f"Reduction '{reduction}' not found in data")

    coords = adata.obsm[reduction_key]

    # Get colors
    if color_by not in adata.obs.columns:
        raise ValueError(f"Column '{color_by}' not found in metadata")

    color_data = adata.obs[color_by].astype(str)

    # Create dataframe for plotting
    plot_df = pd.DataFrame({
        f'{reduction.upper()}_1': coords[:, 0],
        f'{reduction.upper()}_2': coords[:, 1],
        'color': color_data,
        'cell_id': adata.obs_names
    })

    # Get unique categories and assign colors
    categories = sorted(plot_df['color'].unique())
    n_categories = len(categories)
    colors = get_color_palette(n_categories, "discrete")
    color_map = dict(zip(categories, colors))

    # Create figure
    fig = go.Figure()

    for i, category in enumerate(categories):
        mask = plot_df['color'] == category
        fig.add_trace(go.Scatter(
            x=plot_df.loc[mask, f'{reduction.upper()}_1'],
            y=plot_df.loc[mask, f'{reduction.upper()}_2'],
            mode='markers',
            name=str(category),
            marker=dict(
                size=point_size,
                color=color_map[category],
                line=dict(width=0)
            ),
            text=plot_df.loc[mask, 'cell_id'],
            hovertemplate='<b>%{text}</b><br>' +
            f'{reduction.upper()}_1: %{{x:.2f}}<br>' +
            f'{reduction.upper()}_2: %{{y:.2f}}<br>' +
            f'{color_by}: {category}<extra></extra>'
        ))

    # Update layout
    fig.update_layout(
        title=title or f"Cluster Distribution ({reduction.upper()})",
        xaxis_title=f"{reduction.upper()}_1",
        yaxis_title=f"{reduction.upper()}_2",
        showlegend=show_legend,
        hovermode='closest',
        plot_bgcolor='white',
        width=700,
        height=600,
        margin=dict(l=40, r=40, t=60, b=40)
    )

    # Update axes
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')

    return fig


def create_feature_plot(adata, feature: str, reduction: str = "umap",
                        title: str = None, colorscale: str = "Viridis") -> go.Figure:
    """
    Create a feature plot showing gene expression (equivalent to FeaturePlot).

    Args:
        adata: AnnData object
        feature: Gene or feature to plot
        reduction: Reduction to plot
        title: Plot title
        colorscale: Plotly colorscale name

    Returns:
        Plotly figure object
    """
    # Get coordinates
    reduction_key = f"X_{reduction}"
    if reduction_key not in adata.obsm:
        raise ValueError(f"Reduction '{reduction}' not found in data")

    coords = adata.obsm[reduction_key]

    # Get feature values
    if feature in adata.var_names:
        # Gene expression
        gene_idx = adata.var_names.get_loc(feature)
        if hasattr(adata.X, 'toarray'):
            values = adata.X[:, gene_idx].toarray().flatten()
        else:
            values = adata.X[:, gene_idx].flatten()
    elif feature in adata.obs.columns:
        # Metadata column
        values = adata.obs[feature].values
    else:
        raise ValueError(f"Feature '{feature}' not found")

    # Create figure
    fig = go.Figure(data=go.Scatter(
        x=coords[:, 0],
        y=coords[:, 1],
        mode='markers',
        marker=dict(
            size=3,
            color=values,
            colorscale=colorscale,
            showscale=True,
            colorbar=dict(title=feature),
            line=dict(width=0)
        ),
        text=adata.obs_names,
        hovertemplate='<b>%{text}</b><br>' +
        f'{reduction.upper()}_1: %{{x:.2f}}<br>' +
        f'{reduction.upper()}_2: %{{y:.2f}}<br>' +
        f'{feature}: %{{marker.color:.2f}}<extra></extra>'
    ))

    # Update layout
    fig.update_layout(
        title=title or f"{feature} Expression",
        xaxis_title=f"{reduction.upper()}_1",
        yaxis_title=f"{reduction.upper()}_2",
        hovermode='closest',
        plot_bgcolor='white',
        width=700,
        height=600
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')

    return fig


def create_cluster_tree_plot(adata, cluster_key: str = "leiden") -> go.Figure:
    """
    Create a hierarchical cluster tree plot.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist

    # Get cluster assignments
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found")

    clusters = adata.obs[cluster_key].astype(str)
    unique_clusters = sorted(clusters.unique())

    # Calculate cluster centroids in PCA space
    if 'X_pca' not in adata.obsm:
        raise ValueError("PCA not found. Run PCA first.")

    pca_coords = adata.obsm['X_pca'][:, :30]  # Use first 30 PCs

    # Calculate mean PCA coordinates for each cluster
    centroids = []
    for cluster in unique_clusters:
        mask = clusters == cluster
        centroid = pca_coords[mask].mean(axis=0)
        centroids.append(centroid)

    centroids = np.array(centroids)

    # Perform hierarchical clustering
    linkage_matrix = linkage(centroids, method='ward')

    # Create dendrogram
    fig = go.Figure()

    # Use scipy's dendrogram function to get coordinates
    dend = dendrogram(linkage_matrix, labels=unique_clusters, no_plot=True)

    # Extract dendrogram data
    icoord = np.array(dend['icoord'])
    dcoord = np.array(dend['dcoord'])

    # Plot dendrogram lines
    for i in range(len(icoord)):
        fig.add_trace(go.Scatter(
            x=icoord[i],
            y=dcoord[i],
            mode='lines',
            line=dict(color='black', width=2),
            hoverinfo='skip',
            showlegend=False
        ))

    # Add cluster labels
    labels = dend['ivl']
    label_coords = np.array(dend['icoord'])[:, [0, -1]].mean(axis=1)

    # Update layout
    fig.update_layout(
        title="Cluster Dendrogram",
        xaxis=dict(
            title="Cluster",
            ticktext=labels,
            tickvals=[10 + i * 10 for i in range(len(labels))],
            showgrid=False
        ),
        yaxis=dict(
            title="Height",
            showgrid=True,
            gridcolor='lightgray'
        ),
        plot_bgcolor='white',
        height=600,
        showlegend=False
    )

    return fig


def calculate_cluster_stats(adata, cluster_key: str = "leiden") -> pd.DataFrame:
    """
    Calculate cluster statistics.

    Args:
        adata: AnnData object
        cluster_key: Column containing cluster assignments

    Returns:
        DataFrame with cluster statistics
    """
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found")

    clusters = adata.obs[cluster_key]
    counts = clusters.value_counts().sort_index()
    percentages = (counts / len(clusters) * 100).round(2)

    # Build base stats
    stats_df = pd.DataFrame({
        'Cluster': counts.index.astype(str),
        'Cell_Count': counts.values,
        'Percentage': percentages.values
    })

    # Add optional columns if they exist
    if 'n_genes_by_counts' in adata.obs.columns:
        mean_genes = adata.obs.groupby(cluster_key, observed=True)[
            'n_genes_by_counts'].mean().round(0)
        stats_df['Mean_Genes'] = mean_genes.reindex(counts.index).values

    if 'total_counts' in adata.obs.columns:
        mean_counts = adata.obs.groupby(cluster_key, observed=True)[
            'total_counts'].mean().round(0)
        stats_df['Mean_Counts'] = mean_counts.reindex(counts.index).values

    # Add mean mitochondrial percentage if available
    if 'pct_counts_mt' in adata.obs.columns:
        mean_mt = adata.obs.groupby(cluster_key, observed=True)[
            'pct_counts_mt'].mean().round(2)
        stats_df['Mean_MT_Pct'] = mean_mt.reindex(counts.index).values

    return stats_df

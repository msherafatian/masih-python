"""
Marker gene visualization utilities.

This module provides functions for creating heatmaps, dot plots,
and other visualizations specific to marker genes.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Optional, List

from masih.utils.plotting import get_color_palette


def create_marker_heatmap(adata, genes: List[str], cluster_key: str = 'leiden',
                          max_genes: int = 50) -> go.Figure:
    """
    Create a heatmap of marker gene expression across clusters.

    Args:
        adata: AnnData object
        genes: List of genes to plot
        cluster_key: Column containing cluster assignments
        max_genes: Maximum number of genes to display

    Returns:
        Plotly figure object
    """
    # Limit number of genes
    if len(genes) > max_genes:
        genes = genes[:max_genes]

    # Get gene indices
    gene_indices = [adata.var_names.get_loc(
        g) for g in genes if g in adata.var_names]
    genes_found = [genes[i]
                   for i, g in enumerate(genes) if g in adata.var_names]

    if len(genes_found) == 0:
        # Return empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text="No genes found in dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        return fig

    # Get expression data
    if hasattr(adata.X, 'toarray'):
        expr_data = adata.X[:, gene_indices].toarray()
    else:
        expr_data = adata.X[:, gene_indices]

    # Get cluster assignments
    clusters = adata.obs[cluster_key].values
    unique_clusters = sorted(adata.obs[cluster_key].unique())

    # Calculate mean expression per cluster
    cluster_means = []
    for cluster in unique_clusters:
        mask = clusters == cluster
        cluster_mean = expr_data[mask].mean(axis=0)
        cluster_means.append(cluster_mean)

    cluster_means = np.array(cluster_means)

    # Z-score normalization across clusters
    from scipy.stats import zscore
    cluster_means_scaled = zscore(cluster_means, axis=0)

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=cluster_means_scaled.T,
        x=[str(c) for c in unique_clusters],
        y=genes_found,
        colorscale='RdBu_r',
        zmid=0,
        colorbar=dict(title="Scaled<br>Expression"),
        hovertemplate='Cluster: %{x}<br>Gene: %{y}<br>Expr: %{z:.2f}<extra></extra>'
    ))

    # Calculate appropriate height based on number of genes
    height = max(600, len(genes_found) * 20 + 100)

    fig.update_layout(
        title="Top Marker Genes Heatmap",
        xaxis_title="Cluster",
        yaxis_title="Gene",
        height=height,
        width=800,
        yaxis=dict(tickfont=dict(size=max(6, min(10, 300/len(genes_found))))),
        xaxis=dict(side='bottom')
    )

    return fig


def create_dot_plot(adata, genes: List[str], cluster_key: str = 'leiden',
                    max_genes: int = 30) -> go.Figure:
    """
    Create a dot plot showing marker gene expression and percentage.

    Args:
        adata: AnnData object
        genes: List of genes to plot
        cluster_key: Column containing cluster assignments
        max_genes: Maximum number of genes to display

    Returns:
        Plotly figure object
    """
    # Limit number of genes
    if len(genes) > max_genes:
        genes = genes[:max_genes]

    # Get gene indices
    gene_indices = [adata.var_names.get_loc(
        g) for g in genes if g in adata.var_names]
    genes_found = [genes[i]
                   for i, g in enumerate(genes) if g in adata.var_names]

    if len(genes_found) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No genes found in dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="gray")
        )
        return fig

    # Get expression data
    if hasattr(adata.X, 'toarray'):
        expr_data = adata.X[:, gene_indices].toarray()
    else:
        expr_data = adata.X[:, gene_indices]

    # Get clusters
    clusters = adata.obs[cluster_key].values
    unique_clusters = sorted(adata.obs[cluster_key].unique())

    # Calculate mean expression and percentage for each cluster
    plot_data = []

    for i, gene in enumerate(genes_found):
        for cluster in unique_clusters:
            mask = clusters == cluster
            gene_expr = expr_data[mask, i]

            mean_expr = gene_expr.mean()
            pct_expr = (gene_expr > 0).sum() / len(gene_expr) * 100

            plot_data.append({
                'gene': gene,
                'cluster': str(cluster),
                'mean_expr': mean_expr,
                'pct_expr': pct_expr
            })

    plot_df = pd.DataFrame(plot_data)

    # Create figure
    fig = go.Figure()

    # Add scatter trace for each gene
    for gene in genes_found:
        gene_data = plot_df[plot_df['gene'] == gene]

        fig.add_trace(go.Scatter(
            x=gene_data['cluster'],
            y=[gene] * len(gene_data),
            mode='markers',
            marker=dict(
                size=gene_data['pct_expr'],
                sizemode='diameter',
                sizeref=2,
                color=gene_data['mean_expr'],
                colorscale='Viridis',
                showscale=True if gene == genes_found[0] else False,
                colorbar=dict(title="Mean<br>Expression", x=1.05),
                line=dict(width=1, color='white')
            ),
            text=gene_data.apply(lambda x:
                                 f"Cluster: {x['cluster']}<br>" +
                                 f"Gene: {x['gene']}<br>" +
                                 f"Mean Expr: {x['mean_expr']:.2f}<br>" +
                                 f"% Expressed: {x['pct_expr']:.1f}%", axis=1),
            hovertemplate='%{text}<extra></extra>',
            showlegend=False
        ))

    # Calculate height
    height = max(600, len(genes_found) * 30 + 100)

    fig.update_layout(
        title="Marker Genes Dot Plot",
        xaxis_title="Cluster",
        yaxis_title="Gene",
        height=height,
        width=800,
        yaxis=dict(
            categoryorder='array',
            categoryarray=genes_found[::-1],  # Reverse order
            tickfont=dict(size=10)
        ),
        hovermode='closest'
    )

    return fig


def create_violin_plot(adata, gene: str, cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a violin plot for a specific gene across clusters.

    Args:
        adata: AnnData object
        gene: Gene name
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if gene not in adata.var_names:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Gene '{gene}' not found in dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    # Get gene expression
    gene_idx = adata.var_names.get_loc(gene)
    if hasattr(adata.X, 'toarray'):
        expr = adata.X[:, gene_idx].toarray().flatten()
    else:
        expr = adata.X[:, gene_idx].flatten()

    # Get clusters
    clusters = adata.obs[cluster_key].astype(str)
    unique_clusters = sorted(clusters.unique())

    # Get colors
    n_clusters = len(unique_clusters)
    colors = get_color_palette(n_clusters, "discrete")
    color_map = dict(zip(unique_clusters, colors))

    # Create figure
    fig = go.Figure()

    for cluster in unique_clusters:
        mask = clusters == cluster
        cluster_expr = expr[mask]

        fig.add_trace(go.Violin(
            y=cluster_expr,
            x=[cluster] * len(cluster_expr),
            name=cluster,
            marker=dict(color=color_map[cluster]),
            box_visible=True,
            meanline_visible=True,
            showlegend=False
        ))

    fig.update_layout(
        title=f"{gene} Expression by Cluster",
        xaxis_title="Cluster",
        yaxis_title="Expression",
        height=500,
        width=800,
        violinmode='group'
    )

    return fig


def create_ridge_plot(adata, gene: str, cluster_key: str = 'leiden') -> go.Figure:
    """
    Create a ridge plot (stacked density plots) for a gene.

    Args:
        adata: AnnData object
        gene: Gene name
        cluster_key: Column containing cluster assignments

    Returns:
        Plotly figure object
    """
    if gene not in adata.var_names:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Gene '{gene}' not found in dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=14, color="red")
        )
        return fig

    # Get gene expression
    gene_idx = adata.var_names.get_loc(gene)
    if hasattr(adata.X, 'toarray'):
        expr = adata.X[:, gene_idx].toarray().flatten()
    else:
        expr = adata.X[:, gene_idx].flatten()

    # Get clusters
    clusters = adata.obs[cluster_key].astype(str)
    unique_clusters = sorted(clusters.unique(), reverse=True)  # Top to bottom

    # Get colors
    colors = get_color_palette(len(unique_clusters), "discrete")
    color_map = dict(zip(unique_clusters, colors))

    fig = go.Figure()

    for i, cluster in enumerate(unique_clusters):
        mask = clusters == cluster
        cluster_expr = expr[mask]

        # Create density using histogram
        fig.add_trace(go.Violin(
            x=cluster_expr,
            y=[i] * len(cluster_expr),
            name=cluster,
            orientation='h',
            marker=dict(color=color_map[cluster]),
            showlegend=False,
            width=0.7
        ))

    fig.update_layout(
        title=f"{gene} Expression Distribution",
        xaxis_title="Expression",
        yaxis=dict(
            ticktext=unique_clusters,
            tickvals=list(range(len(unique_clusters))),
            title="Cluster"
        ),
        height=max(400, len(unique_clusters) * 50),
        width=800,
        violinmode='overlay'
    )

    return fig

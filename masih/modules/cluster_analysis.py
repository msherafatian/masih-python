"""
Cluster Analysis Module for MASIH application.

This module provides interactive visualization and analysis of clustering results,
including dimensionality reduction plots, cluster trees, and statistics.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context, dash_table
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import pandas as pd

from masih.config.settings import AppConfig
from masih.utils.data_store import get_adata
from masih.utils.plotting import (
    create_dim_plot,
    create_feature_plot,
    create_cluster_tree_plot,
    calculate_cluster_stats
)


def create_cluster_analysis_layout():
    """
    Create the layout for the Cluster Analysis module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H3("Cluster Analysis", className="mb-3"),
                html.P(
                    "Visualize and analyze clustering results with interactive plots.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Main visualization row
        dbc.Row([
            # Cluster visualization
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([
                        html.H5("Cluster Visualization", className="mb-0")
                    ]),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Reduction:"),
                                dcc.Dropdown(
                                    id='cluster-reduction-select',
                                    options=[],
                                    value='umap',
                                    clearable=False
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Color by:"),
                                dcc.Dropdown(
                                    id='cluster-color-by-select',
                                    options=[],
                                    value='leiden',
                                    clearable=False
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Point size:"),
                                dcc.Slider(
                                    id='cluster-point-size',
                                    min=1,
                                    max=10,
                                    step=1,
                                    value=3,
                                    marks={i: str(i) for i in [1, 3, 5, 7, 10]}
                                )
                            ], md=4)
                        ], className="mb-3"),

                        dcc.Loading(
                            id="loading-cluster-plot",
                            type="default",
                            children=dcc.Graph(
                                id='cluster-plot',
                                config={'displayModeBar': True},
                                style={'height': '600px'}
                            )
                        )
                    ])
                ])
            ], md=8),

            # Cluster tree
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([
                        html.H5("Cluster Tree", className="mb-0")
                    ]),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-cluster-tree",
                            type="default",
                            children=dcc.Graph(
                                id='cluster-tree-plot',
                                config={'displayModeBar': False},
                                style={'height': '600px'}
                            )
                        )
                    ])
                ])
            ], md=4)
        ], className="mb-4"),

        # Cluster statistics
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader([
                        html.H5("Cluster Statistics", className="mb-0")
                    ]),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-cluster-stats",
                            type="default",
                            children=dash_table.DataTable(
                                id='cluster-stats-table',
                                columns=[],
                                data=[],
                                style_table={'overflowX': 'auto'},
                                style_cell={
                                    'textAlign': 'left',
                                    'padding': '10px'
                                },
                                style_header={
                                    'backgroundColor': 'rgb(230, 230, 230)',
                                    'fontWeight': 'bold'
                                },
                                style_data_conditional=[
                                    {
                                        'if': {'row_index': 'odd'},
                                        'backgroundColor': 'rgb(248, 248, 248)'
                                    }
                                ],
                                page_size=15,
                                sort_action='native',
                                filter_action='native'
                            )
                        )
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


def register_cluster_analysis_callbacks(app):
    """
    Register all callbacks for the cluster analysis module.

    Args:
        app: Dash application instance
    """

    # Callback to populate dropdowns when data is loaded
    @app.callback(
        [Output('cluster-reduction-select', 'options'),
         Output('cluster-reduction-select', 'value'),
         Output('cluster-color-by-select', 'options'),
         Output('cluster-color-by-select', 'value')],
        Input('adata-store', 'data')
    )
    def update_dropdown_options(adata_metadata):
        """Update dropdown options based on available data."""
        if adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get available reductions
            reductions = []
            if 'X_umap' in adata.obsm:
                reductions.append({'label': 'UMAP', 'value': 'umap'})
            if 'X_tsne' in adata.obsm:
                reductions.append({'label': 't-SNE', 'value': 'tsne'})
            if 'X_pca' in adata.obsm:
                reductions.append({'label': 'PCA', 'value': 'pca'})

            reduction_value = reductions[0]['value'] if reductions else 'umap'

            # Get metadata columns for coloring
            color_options = []

            # Add cluster columns
            if 'leiden' in adata.obs.columns:
                color_options.append(
                    {'label': 'Leiden Clusters', 'value': 'leiden'})
            if 'louvain' in adata.obs.columns:
                color_options.append(
                    {'label': 'Louvain Clusters', 'value': 'louvain'})

            # Add cell cycle phase
            if 'phase' in adata.obs.columns:
                color_options.append(
                    {'label': 'Cell Cycle Phase', 'value': 'phase'})

            # Add other categorical columns
            for col in adata.obs.columns:
                if col not in ['leiden', 'louvain', 'phase', 'n_genes', 'n_genes_by_counts',
                               'total_counts', 'total_counts_mt', 'pct_counts_mt']:
                    # Check if categorical or few unique values
                    if adata.obs[col].dtype == 'object' or adata.obs[col].nunique() < 20:
                        color_options.append({'label': col, 'value': col})

            color_value = 'leiden' if 'leiden' in adata.obs.columns else color_options[
                0]['value']

            return reductions, reduction_value, color_options, color_value

        except Exception as e:
            print(f"Error updating dropdowns: {e}")
            raise PreventUpdate

    # Callback to update cluster plot

    @app.callback(
        Output('cluster-plot', 'figure'),
        [Input('cluster-reduction-select', 'value'),
         Input('cluster-color-by-select', 'value'),
         Input('cluster-point-size', 'value')],
        State('adata-store', 'data')
    )
    def update_cluster_plot(reduction, color_by, point_size, adata_metadata):
        """Update the main cluster visualization plot."""
        if adata_metadata is None or reduction is None or color_by is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Create the plot
            fig = create_dim_plot(
                adata,
                reduction=reduction,
                color_by=color_by,
                point_size=point_size,
                title=f"Cluster Distribution ({reduction.upper()})"
            )

            return fig

        except Exception as e:
            print(f"Error creating cluster plot: {e}")
            # Return empty figure with error message
            fig = go.Figure()
            fig.add_annotation(
                text=f"Error creating plot: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="red")
            )
            return fig

    # Callback to update cluster tree

    @app.callback(
        Output('cluster-tree-plot', 'figure'),
        [Input('adata-store', 'data'),
         Input('cluster-color-by-select', 'value')]
    )
    def update_cluster_tree(adata_metadata, cluster_key):
        """Update the cluster tree dendrogram."""
        if adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Use the selected cluster key, default to leiden
            if cluster_key is None or cluster_key not in adata.obs.columns:
                cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Create the tree plot
            fig = create_cluster_tree_plot(adata, cluster_key=cluster_key)

            return fig

        except Exception as e:
            print(f"Error creating cluster tree: {e}")
            # Return empty figure with error message
            fig = go.Figure()
            fig.add_annotation(
                text=f"Unable to create cluster tree.\n{str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=12, color="gray")
            )
            fig.update_layout(height=600)
            return fig

# Callback to update cluster statistics table
    @app.callback(
        [Output('cluster-stats-table', 'columns'),
         Output('cluster-stats-table', 'data')],
        [Input('adata-store', 'data'),
         Input('cluster-color-by-select', 'value')]
    )
    def update_cluster_stats(adata_metadata, cluster_key):
        """Update the cluster statistics table."""
        if adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Use the selected cluster key, default to leiden
            if cluster_key is None or cluster_key not in adata.obs.columns:
                cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Calculate statistics
            stats_df = calculate_cluster_stats(adata, cluster_key=cluster_key)

            # Dynamically create columns based on what's in the dataframe
            columns = [{'name': col.replace('_', ' '), 'id': col}
                       for col in stats_df.columns]

            data = stats_df.to_dict('records')

            return columns, data

        except Exception as e:
            print(f"Error calculating cluster stats: {e}")
            import traceback
            traceback.print_exc()
            return [], []

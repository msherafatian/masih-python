"""
Marker Genes Module for MASIH application.

This module provides differential expression analysis to identify
marker genes that characterize each cluster.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context, dash_table
from dash.exceptions import PreventUpdate
import pandas as pd
import io
import base64
import plotly.graph_objects as go

from masih.config.settings import AppConfig
from masih.utils.data_store import get_adata, data_store
from masih.analysis.marker_utils import (
    find_all_markers,
    find_pairwise_markers,
    get_top_markers,
    filter_markers,
    summarize_markers
)
from masih.utils.marker_plotting import (
    create_marker_heatmap,
    create_dot_plot,
    create_violin_plot,
    create_ridge_plot
)
from masih.utils.plotting import create_feature_plot


def create_markers_layout():
    """
    Create the layout for the Marker Genes module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval to trigger reload when tab becomes visible
        dcc.Interval(
            id='markers-reload-interval',
            interval=500,  # Check every 500ms
            max_intervals=1,  # Fire only once
            disabled=False
        ),
        dbc.Row([
            dbc.Col([
                html.H3("Marker Gene Analysis", className="mb-3"),
                html.P(
                    "Identify genes that characterize each cluster through differential expression analysis.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Settings Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Marker Gene Analysis Settings", className="mb-0")),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Statistical test:"),
                                dcc.Dropdown(
                                    id='marker-test-method',
                                    options=[
                                        {'label': 'Wilcoxon Rank Sum',
                                            'value': 'wilcoxon'},
                                        {'label': "Student's t-test",
                                            'value': 't-test'},
                                        {'label': 'Logistic Regression',
                                            'value': 'logreg'}
                                    ],
                                    value='wilcoxon',
                                    clearable=False
                                )
                            ], md=3),

                            dbc.Col([
                                dbc.Label("Comparison type:"),
                                dbc.RadioItems(
                                    id='marker-comparison-type',
                                    options=[
                                        {'label': 'One vs All',
                                            'value': 'one_vs_all'},
                                        {'label': 'Pairwise', 'value': 'pairwise'}
                                    ],
                                    value='one_vs_all',
                                    inline=False
                                )
                            ], md=3),

                            dbc.Col([
                                dbc.Label("Min log fold change:"),
                                dbc.Input(
                                    id='marker-min-logfc',
                                    type='number',
                                    value=AppConfig.DEFAULT_MARKER_LOGFC,
                                    min=0,
                                    max=2,
                                    step=0.1
                                )
                            ], md=3),

                            dbc.Col([
                                dbc.Label("Min percent expressed:"),
                                dbc.Input(
                                    id='marker-min-pct',
                                    type='number',
                                    value=AppConfig.DEFAULT_MARKER_MINPCT,
                                    min=0,
                                    max=1,
                                    step=0.05
                                )
                            ], md=3)
                        ], className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Checklist(
                                    id='marker-only-pos',
                                    options=[
                                        {'label': ' Only positive markers', 'value': 'only_pos'}],
                                    value=['only_pos'],
                                    inline=True
                                )
                            ], md=3),

                            dbc.Col([
                                dbc.Label("Top N genes per cluster:"),
                                dbc.Input(
                                    id='marker-top-n',
                                    type='number',
                                    value=AppConfig.DEFAULT_MARKER_TOPN,
                                    min=1,
                                    max=50,
                                    step=1
                                )
                            ], md=3),

                            dbc.Col([
                                dbc.Label("Max genes to display:"),
                                dbc.Input(
                                    id='marker-max-display',
                                    type='number',
                                    value=AppConfig.DEFAULT_MARKER_DISPLAY_GENES,
                                    min=10,
                                    max=100,
                                    step=5
                                )
                            ], md=3),

                            dbc.Col([
                                html.Div(id='marker-pairwise-selectors', style={'display': 'none'}, children=[
                                    dbc.Label("Cluster 1:"),
                                    dcc.Dropdown(
                                        id='marker-cluster1', options=[], clearable=False),
                                    dbc.Label("Cluster 2:", className="mt-2"),
                                    dcc.Dropdown(
                                        id='marker-cluster2', options=[], clearable=False)
                                ])
                            ], md=3)
                        ]),

                        html.Hr(),

                        dbc.Button(
                            "Find Marker Genes",
                            id='calculate-markers-button',
                            color='success',
                            size='lg',
                            className="mt-2"
                        ),

                        html.Div(id='marker-calculation-status',
                                 className="mt-3"),

                        # Quick actions (shown after calculation)
                        html.Div(id='marker-quick-actions', style={'display': 'none'}, children=[
                            html.Hr(),
                            html.H5("Quick Actions:", className="mt-3"),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        "Download All Markers (Excel)",
                                        id='download-all-markers-button',
                                        color='info',
                                        className="w-100"
                                    ),
                                    dcc.Download(id='download-all-markers')
                                ], md=4),

                                dbc.Col([
                                    dbc.Button(
                                        "Download Top Markers (CSV)",
                                        id='download-top-markers-button',
                                        color='info',
                                        className="w-100"
                                    ),
                                    dcc.Download(id='download-top-markers')
                                ], md=4),

                                dbc.Col([
                                    dbc.Button(
                                        "Clear Results",
                                        id='clear-markers-button',
                                        color='warning',
                                        className="w-100"
                                    )
                                ], md=4)
                            ])
                        ])
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Marker Table
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Marker Gene Table", className="mb-0")),
                    dbc.CardBody([
                        html.Div(id='marker-table-container', children=[
                            html.H5(
                                "No markers calculated yet. Click 'Find Marker Genes' to start analysis.",
                                style={'color': '#999',
                                       'textAlign': 'center', 'padding': '50px'}
                            )
                        ])
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Visualizations
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Top Markers Heatmap", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-marker-heatmap",
                            children=dcc.Graph(
                                id='marker-heatmap', config={'displayModeBar': True})
                        )
                    ]),
                    dbc.CardFooter(id='marker-heatmap-info', children="")
                ], className="mb-4")
            ], width=12)
        ]),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Marker Dot Plot", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-marker-dotplot",
                            children=dcc.Graph(
                                id='marker-dotplot', config={'displayModeBar': True})
                        )
                    ])
                ], className="mb-4")
            ], width=12)
        ]),

        # Selected Marker Visualization
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Selected Marker Visualization", className="mb-0")),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Select marker gene:"),
                                dcc.Dropdown(
                                    id='visualize-marker-gene',
                                    options=[],
                                    placeholder="Select a gene..."
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Plot type:"),
                                dbc.RadioItems(
                                    id='marker-plot-type',
                                    options=[
                                        {'label': 'Feature Plot',
                                            'value': 'feature'},
                                        {'label': 'Violin Plot',
                                            'value': 'violin'},
                                        {'label': 'Ridge Plot', 'value': 'ridge'}
                                    ],
                                    value='feature',
                                    inline=True
                                )
                            ], md=4),

                            dbc.Col([
                                html.Br(),
                                dbc.Button(
                                    "Update Plot",
                                    id='update-marker-plot-button',
                                    color='primary',
                                    className="mt-2"
                                )
                            ], md=4)
                        ], className="mb-3"),

                        dcc.Loading(
                            id="loading-selected-marker",
                            children=dcc.Graph(
                                id='selected-marker-plot',
                                config={'displayModeBar': True},
                                style={'height': '500px'}
                            )
                        )
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


def _rehydrate_markers_view(adata, markers_data):
    """
    Shared function to rebuild markers view from stored data.

    This ensures both callbacks (initial calculation and tab-switch reload)
    produce identical output.

    Args:
        adata: AnnData object
        markers_data: Dictionary containing markers, top_markers, and parameters

    Returns:
        Tuple of (table, quick_actions_style, heatmap_fig, dotplot_fig, gene_options, heatmap_info)
    """
    markers_df = pd.DataFrame(markers_data['markers'])
    params = markers_data['parameters']
    top_markers_df = pd.DataFrame(markers_data['top_markers'])

    if markers_df.empty:
        raise PreventUpdate

    cluster_key = params['cluster_key']
    max_display = params['max_display']

    # Create table
    table = dash_table.DataTable(
        id='marker-genes-table',
        columns=[
            {'name': 'Cluster', 'id': 'cluster'},
            {'name': 'Gene', 'id': 'gene'},
            {'name': 'Log2FC', 'id': 'avg_log2FC',
                'type': 'numeric', 'format': {'specifier': '.3f'}},
            {'name': 'P-value', 'id': 'p_val', 'type': 'numeric',
                'format': {'specifier': '.2e'}},
            {'name': 'Adj P-value', 'id': 'p_val_adj',
                'type': 'numeric', 'format': {'specifier': '.2e'}},
            {'name': 'Pct.1', 'id': 'pct.1', 'type': 'numeric',
                'format': {'specifier': '.3f'}},
            {'name': 'Pct.2', 'id': 'pct.2', 'type': 'numeric',
                'format': {'specifier': '.3f'}}
        ],
        data=markers_df.round(3).to_dict('records'),
        style_table={'overflowX': 'auto'},
        style_cell={'textAlign': 'left',
                    'padding': '10px', 'minWidth': '100px'},
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
        style_data_conditional=[
            {'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'},
            {'if': {'filter_query': '{avg_log2FC} > 1', 'column_id': 'avg_log2FC'},
             'backgroundColor': '#d4edda', 'color': '#155724'}
        ],
        page_size=25,
        sort_action='native',
        filter_action='native',
        export_format='csv'
    )

    # Create visualizations
    top_genes = top_markers_df['gene'].unique()[:max_display]
    heatmap_fig = create_marker_heatmap(
        adata, list(top_genes), cluster_key, max_display)
    heatmap_info = f"Showing top {len(top_genes)} marker genes"

    dotplot_genes = top_markers_df.groupby('cluster').head(3)[
        'gene'].unique()[:30]
    dotplot_fig = create_dot_plot(adata, list(
        dotplot_genes), cluster_key, max_genes=30)

    gene_options = [{'label': g, 'value': g}
                    for g in markers_df['gene'].unique()]

    return table, {'display': 'block'}, heatmap_fig, dotplot_fig, gene_options, heatmap_info


def register_markers_callbacks(app):
    """
    Register all callbacks for the marker genes module.

    Args:
        app: Dash application instance
    """

    from datetime import datetime
    import traceback

    # Callback to show/hide pairwise selectors

    @app.callback(
        Output('marker-pairwise-selectors', 'style'),
        Input('marker-comparison-type', 'value')
    )
    def toggle_pairwise_selectors(comparison_type):
        """Show pairwise cluster selectors when pairwise is selected."""
        if comparison_type == 'pairwise':
            return {'display': 'block'}
        return {'display': 'none'}

    # Callback to populate cluster dropdowns

    @app.callback(
        [Output('marker-cluster1', 'options'),
         Output('marker-cluster1', 'value'),
         Output('marker-cluster2', 'options'),
         Output('marker-cluster2', 'value')],
        Input('adata-store', 'data')
    )
    def update_cluster_options(adata_metadata):
        """Update cluster dropdown options."""
        if adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'
            clusters = sorted(adata.obs[cluster_key].unique())

            options = [{'label': str(c), 'value': str(c)} for c in clusters]

            # Default values
            value1 = str(clusters[0]) if len(clusters) > 0 else None
            value2 = str(clusters[1]) if len(clusters) > 1 else value1

            return options, value1, options, value2

        except Exception as e:
            print(f"Error updating cluster options: {e}")
            raise PreventUpdate

   # Main callback to calculate markers
    @app.callback(
        [Output('marker-calculation-status', 'children'),
         Output('marker-quick-actions', 'style'),
         Output('marker-table-container', 'children'),
         Output('marker-heatmap', 'figure'),
         Output('marker-dotplot', 'figure'),
         Output('visualize-marker-gene', 'options'),
         Output('marker-heatmap-info', 'children'),
         Output('markers-store', 'data'),  # NEW - store the results
         # NEW - update flags
         Output('processed-flags', 'data', allow_duplicate=True)],
        Input('calculate-markers-button', 'n_clicks'),
        [State('adata-store', 'data'),
         State('marker-test-method', 'value'),
         State('marker-comparison-type', 'value'),
         State('marker-min-logfc', 'value'),
         State('marker-min-pct', 'value'),
         State('marker-only-pos', 'value'),
         State('marker-top-n', 'value'),
         State('marker-max-display', 'value'),
         State('marker-cluster1', 'value'),
         State('marker-cluster2', 'value'),
         State('processed-flags', 'data')],
        prevent_initial_call=True
    )
    def calculate_markers(n_clicks, adata_metadata, method, comparison_type,
                          min_logfc, min_pct, only_pos, top_n, max_display,
                          cluster1, cluster2, flags):
        """Calculate marker genes."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            only_pos_bool = 'only_pos' in only_pos if only_pos else False

            # Calculate markers
            if comparison_type == 'one_vs_all':
                markers = find_all_markers(
                    adata,
                    cluster_key=cluster_key,
                    method=method,
                    min_pct=min_pct,
                    logfc_threshold=min_logfc,
                    only_pos=only_pos_bool
                )
            else:  # pairwise
                markers = find_pairwise_markers(
                    adata,
                    cluster_key=cluster_key,
                    ident_1=cluster1,
                    ident_2=cluster2,
                    method=method,
                    min_pct=min_pct,
                    logfc_threshold=min_logfc,
                    only_pos=only_pos_bool
                )

            if len(markers) == 0:
                status = dbc.Alert(
                    "No markers found with current criteria. Try relaxing the thresholds.",
                    color="warning"
                )
                empty_markers = {'markers': [],
                                 'top_markers': [], 'parameters': {}}
                return status, {'display': 'none'}, html.H5("No markers found",
                                                            style={'color': '#999', 'textAlign': 'center', 'padding': '50px'}), {}, {}, [], "", empty_markers, flags

            # Get top markers
            top_markers = get_top_markers(markers, n=top_n, by_cluster=True)

            # Store in markers-store for persistence
            markers_data = {
                'markers': markers.to_dict('records'),
                'top_markers': top_markers.to_dict('records'),
                'parameters': {
                    'method': method,
                    'comparison_type': comparison_type,
                    'min_logfc': min_logfc,
                    'min_pct': min_pct,
                    'only_pos': only_pos_bool,
                    'top_n': top_n,
                    'max_display': max_display,
                    'cluster_key': cluster_key
                }
            }

            # Also store for download
            data_store.store_additional_data(
                session_id, 'markers', markers.to_dict('records'))
            data_store.store_additional_data(
                session_id, 'top_markers', top_markers.to_dict('records'))

            # Create table
            table = dash_table.DataTable(
                id='marker-genes-table',
                columns=[
                    {'name': 'Cluster', 'id': 'cluster'},
                    {'name': 'Gene', 'id': 'gene'},
                    {'name': 'Log2FC', 'id': 'avg_log2FC',
                        'type': 'numeric', 'format': {'specifier': '.3f'}},
                    {'name': 'P-value', 'id': 'p_val', 'type': 'numeric',
                        'format': {'specifier': '.2e'}},
                    {'name': 'Adj P-value', 'id': 'p_val_adj',
                        'type': 'numeric', 'format': {'specifier': '.2e'}},
                    {'name': 'Pct.1', 'id': 'pct.1', 'type': 'numeric',
                        'format': {'specifier': '.3f'}},
                    {'name': 'Pct.2', 'id': 'pct.2', 'type': 'numeric',
                        'format': {'specifier': '.3f'}}
                ],
                data=markers.round(3).to_dict('records'),
                style_table={'overflowX': 'auto'},
                style_cell={
                    'textAlign': 'left',
                    'padding': '10px',
                    'minWidth': '100px'
                },
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold'
                },
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(248, 248, 248)'
                    },
                    {
                        'if': {
                            'filter_query': '{avg_log2FC} > 1',
                            'column_id': 'avg_log2FC'
                        },
                        'backgroundColor': '#d4edda',
                        'color': '#155724'
                    }
                ],
                page_size=25,
                sort_action='native',
                filter_action='native',
                export_format='csv'
            )

            # Create heatmap
            top_genes = top_markers['gene'].unique()[:max_display]
            heatmap_fig = create_marker_heatmap(
                adata, list(top_genes), cluster_key, max_display)
            heatmap_info = f"Showing top {len(top_genes)} marker genes"

            # Create dot plot
            dotplot_genes = top_markers.groupby('cluster').head(3)[
                'gene'].unique()[:30]
            dotplot_fig = create_dot_plot(adata, list(
                dotplot_genes), cluster_key, max_genes=30)

            # Gene options for visualization
            gene_options = [{'label': g, 'value': g}
                            for g in markers['gene'].unique()]

            # Status message
            status = dbc.Alert([
                html.I(className="fas fa-check-circle me-2"),
                f"Found {len(markers)} marker genes across {markers['cluster'].nunique()} cluster(s)!"
            ], color="success")

            # Update flags
            if flags is None:
                flags = {}
            flags['markers_calculated'] = True

            return (status, {'display': 'block'}, table, heatmap_fig,
                    dotplot_fig, gene_options, heatmap_info, markers_data, flags)

        except Exception as e:
            print(f"Error calculating markers: {e}")
            traceback.print_exc()

            status = dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error calculating markers: {str(e)}"
            ], color="danger")

            return status, {'display': 'none'}, html.H5("Error",
                                                        style={'color': 'red', 'textAlign': 'center'}), {}, {}, [], "", None, flags

 # Callback 1: Load when markers are first calculated
    @app.callback(
        [Output('marker-table-container', 'children', allow_duplicate=True),
         Output('marker-quick-actions', 'style', allow_duplicate=True),
         Output('marker-heatmap', 'figure', allow_duplicate=True),
         Output('marker-dotplot', 'figure', allow_duplicate=True),
         Output('visualize-marker-gene', 'options', allow_duplicate=True),
         Output('marker-heatmap-info', 'children', allow_duplicate=True)],
        Input('markers-store', 'data'),
        State('adata-store', 'data'),
        prevent_initial_call=True
    )
    def load_stored_markers(markers_data, adata_metadata):
        """Load markers immediately after calculation (markers-store changes)."""
        if not markers_data or not adata_metadata:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            return _rehydrate_markers_view(adata, markers_data)

        except Exception as e:
            print(f"Error loading stored markers: {e}")
            import traceback
            traceback.print_exc()
            raise PreventUpdate

    # Callback 2: Reload when switching back to markers tab

    @app.callback(
        [Output('marker-table-container', 'children', allow_duplicate=True),
         Output('marker-quick-actions', 'style', allow_duplicate=True),
         Output('marker-heatmap', 'figure', allow_duplicate=True),
         Output('marker-dotplot', 'figure', allow_duplicate=True),
         Output('visualize-marker-gene', 'options', allow_duplicate=True),
         Output('marker-heatmap-info', 'children', allow_duplicate=True)],
        Input('main-tabs', 'active_tab'),
        [State('markers-store', 'data'),
         State('adata-store', 'data')],
        prevent_initial_call=True
    )
    def reload_markers_on_tab_switch(active_tab, markers_data, adata_metadata):
        """Reload markers when returning to markers tab."""
        if active_tab != 'markers':
            raise PreventUpdate

        if not markers_data or not adata_metadata:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            print("DEBUG: Reloading markers on tab switch")
            return _rehydrate_markers_view(adata, markers_data)

        except Exception as e:
            print(f"Error reloading markers on tab switch: {e}")
            import traceback
            traceback.print_exc()
            raise PreventUpdate

    # Callback to update selected marker plot

    @app.callback(
        Output('selected-marker-plot', 'figure'),
        Input('update-marker-plot-button', 'n_clicks'),
        [State('visualize-marker-gene', 'value'),
         State('marker-plot-type', 'value'),
         State('adata-store', 'data')],
        prevent_initial_call=True
    )
    def update_selected_marker_plot(n_clicks, gene, plot_type, adata_metadata):
        """Update the selected marker gene plot."""
        if n_clicks is None or gene is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            if plot_type == 'feature':
                fig = create_feature_plot(
                    adata, gene, reduction='umap', title=f"{gene} Expression")
            elif plot_type == 'violin':
                fig = create_violin_plot(adata, gene, cluster_key)
            elif plot_type == 'ridge':
                fig = create_ridge_plot(adata, gene, cluster_key)
            else:
                fig = {}

            return fig

        except Exception as e:
            print(f"Error creating marker plot: {e}")
            fig = go.Figure()
            fig.add_annotation(
                text=f"Error: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="red")
            )
            return fig
# Simple reload trigger using interval

    @app.callback(
        [Output('marker-table-container', 'children', allow_duplicate=True),
         Output('marker-quick-actions', 'style', allow_duplicate=True),
         Output('marker-heatmap', 'figure', allow_duplicate=True),
         Output('marker-dotplot', 'figure', allow_duplicate=True),
         Output('visualize-marker-gene', 'options', allow_duplicate=True),
         Output('marker-heatmap-info', 'children', allow_duplicate=True),
         Output('markers-reload-interval', 'disabled')],
        Input('markers-reload-interval', 'n_intervals'),
        [State('markers-store', 'data'),
         State('adata-store', 'data')],
        prevent_initial_call=True
    )
    def reload_on_interval(n, markers_data, adata_metadata):
        """Reload markers when interval fires (happens when tab is created)."""
        print(f"DEBUG: Interval fired, n={n}")

        if not markers_data or not adata_metadata:
            print("DEBUG: No data to reload")
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            print("DEBUG: Reloading markers via interval")
            result = _rehydrate_markers_view(adata, markers_data)

            # Disable interval after first fire
            return (*result, True)

        except Exception as e:
            print(f"Error in interval reload: {e}")
            import traceback
            traceback.print_exc()
            raise PreventUpdate

    # Callback to clear markers

    @app.callback(
        [Output('marker-calculation-status', 'children', allow_duplicate=True),
         Output('marker-quick-actions', 'style', allow_duplicate=True),
         Output('marker-table-container', 'children', allow_duplicate=True)],
        Input('clear-markers-button', 'n_clicks'),
        State('adata-store', 'data'),
        prevent_initial_call=True
    )
    def clear_markers(n_clicks, adata_metadata):
        """Clear marker results."""
        if n_clicks is None:
            raise PreventUpdate

        if adata_metadata:
            session_id = adata_metadata['session_id']
            data_store.store_additional_data(session_id, 'markers', None)
            data_store.store_additional_data(session_id, 'top_markers', None)

        status = dbc.Alert("Marker results cleared",
                           color="info", duration=3000)
        empty_msg = html.H5(
            "No markers calculated yet. Click 'Find Marker Genes' to start analysis.",
            style={'color': '#999', 'textAlign': 'center', 'padding': '50px'}
        )

        return status, {'display': 'none'}, empty_msg

    # Download callbacks

    @app.callback(
        Output('download-all-markers', 'data'),
        Input('download-all-markers-button', 'n_clicks'),
        State('adata-store', 'data'),
        prevent_initial_call=True
    )
    def download_all_markers(n_clicks, adata_metadata):
        """Download all markers as Excel."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            markers_data = data_store.load_additional_data(
                session_id, 'markers')

            if markers_data is None:
                raise PreventUpdate

            markers_df = pd.DataFrame(markers_data)

            # Create Excel file in memory
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                # Write all markers
                markers_df.to_excel(
                    writer, sheet_name='All_Markers', index=False)

                # Write per cluster if one_vs_all
                if 'cluster' in markers_df.columns:
                    for cluster in sorted(markers_df['cluster'].unique()):
                        cluster_markers = markers_df[markers_df['cluster'] == cluster]
                        sheet_name = f"Cluster_{cluster}"[:31]  # Excel limit
                        cluster_markers.to_excel(
                            writer, sheet_name=sheet_name, index=False)

            output.seek(0)

            return dcc.send_bytes(output.getvalue(),
                                  f"marker_genes_all_{datetime.now().strftime('%Y%m%d')}.xlsx")

        except Exception as e:
            print(f"Error downloading markers: {e}")
            raise PreventUpdate

    @app.callback(
        Output('download-top-markers', 'data'),
        Input('download-top-markers-button', 'n_clicks'),
        State('adata-store', 'data'),
        prevent_initial_call=True
    )
    def download_top_markers(n_clicks, adata_metadata):
        """Download top markers as CSV."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            top_markers_data = data_store.load_additional_data(
                session_id, 'top_markers')

            if top_markers_data is None:
                raise PreventUpdate

            top_markers_df = pd.DataFrame(top_markers_data)

            return dcc.send_data_frame(
                top_markers_df.to_csv,
                f"top_marker_genes_{datetime.now().strftime('%Y%m%d')}.csv",
                index=False
            )

        except Exception as e:
            print(f"Error downloading top markers: {e}")
            raise PreventUpdate

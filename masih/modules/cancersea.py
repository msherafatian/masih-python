"""
CancerSEA Pathways Module for MASIH application.

This module provides functionality to score cancer-related functional states
using the CancerSEA gene sets.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context
from dash.exceptions import PreventUpdate
import pandas as pd
import plotly.graph_objects as go

from masih.config.settings import AppConfig
from masih.utils.data_store import get_adata, data_store
from masih.utils.cancersea_utils import (
    calculate_pathway_score,
    calculate_multiple_pathways,
    calculate_pathway_averages,
    get_pathway_statistics
)
from masih.utils.cancersea_plotting import (
    create_pathway_feature_plot,
    create_pathway_violin_plot,
    create_pathway_heatmap,
    create_pathway_comparison_plot
)


def create_cancersea_layout():
    """
    Create the layout for the CancerSEA Pathways module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval to trigger reload when tab becomes visible
        dcc.Interval(
            id='cancersea-reload-interval',
            interval=500,
            max_intervals=1,
            disabled=False
        ),

        dbc.Row([
            dbc.Col([
                html.H3("CancerSEA Functional State Analysis", className="mb-3"),
                html.P(
                    "Score cancer-related functional states using CancerSEA gene signatures.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Pathway Selection Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Functional State Selection", className="mb-0")),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Select CancerSEA pathway:"),
                                dcc.Dropdown(
                                    id='cancersea-pathway-select',
                                    options=[
                                        {'label': 'Angiogenesis',
                                            'value': 'Angiogenesis'},
                                        {'label': 'Apoptosis',
                                            'value': 'Apoptosis'},
                                        {'label': 'Cell Cycle',
                                            'value': 'Cell_cycle'},
                                        {'label': 'Differentiation',
                                            'value': 'Differentiation'},
                                        {'label': 'DNA Damage',
                                            'value': 'DNA_damage'},
                                        {'label': 'DNA Repair',
                                            'value': 'DNA_repair'},
                                        {'label': 'EMT', 'value': 'EMT'},
                                        {'label': 'Hypoxia', 'value': 'Hypoxia'},
                                        {'label': 'Inflammation',
                                            'value': 'Inflammation'},
                                        {'label': 'Invasion', 'value': 'Invasion'},
                                        {'label': 'Metastasis',
                                            'value': 'Metastasis'},
                                        {'label': 'Proliferation',
                                            'value': 'Proliferation'},
                                        {'label': 'Quiescence',
                                            'value': 'Quiescence'},
                                        {'label': 'Stemness', 'value': 'Stemness'}
                                    ],
                                    value='Proliferation',
                                    clearable=False
                                )
                            ], md=8),

                            dbc.Col([
                                html.Br(),
                                dbc.Button(
                                    "Calculate/Update Score",
                                    id='calculate-cancersea-button',
                                    color='info',
                                    className="w-100"
                                )
                            ], md=4)
                        ]),

                        html.Div(id='cancersea-calculation-status',
                                 className="mt-3"),

                        html.Hr(),

                        # Batch calculation section
                        html.Div([
                            html.H6("Batch Calculate Multiple Pathways:",
                                    className="mb-2"),
                            dbc.Checklist(
                                id='cancersea-batch-select',
                                options=[
                                    {'label': 'Select All', 'value': 'all'}
                                ],
                                value=[],
                                inline=True,
                                className="mb-2"
                            ),
                            dbc.Button(
                                "Calculate Selected Pathways",
                                id='calculate-batch-cancersea-button',
                                color='success',
                                size='sm'
                            )
                        ], id='cancersea-batch-area')
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Visualization Cards
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Functional State Visualization", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-cancersea-feature",
                            children=dcc.Graph(
                                id='cancersea-feature-plot',
                                config={'displayModeBar': True},
                                style={'height': '500px'}
                            )
                        )
                    ])
                ])
            ], md=6),

            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Violin Plot by Cluster", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-cancersea-violin",
                            children=dcc.Graph(
                                id='cancersea-violin-plot',
                                config={'displayModeBar': True},
                                style={'height': '500px'}
                            )
                        )
                    ])
                ])
            ], md=6)
        ], className="mb-4"),

        # Heatmap Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Functional State Heatmap", className="mb-0")),
                    dbc.CardBody([
                        html.P(
                            "Heatmap will appear after calculating at least 2 pathways",
                            className="text-muted mb-3"
                        ),
                        dcc.Loading(
                            id="loading-cancersea-heatmap",
                            children=dcc.Graph(
                                id='cancersea-heatmap',
                                config={'displayModeBar': True},
                                figure=go.Figure()
                            )
                        )
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Statistics Table
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Pathway Score Statistics", className="mb-0")),
                    dbc.CardBody([
                        html.Div(id='cancersea-stats-table')
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


def register_cancersea_callbacks(app):
    """
    Register all callbacks for the CancerSEA pathways module.

    Args:
        app: Dash application instance
    """

    from datetime import datetime
    import traceback

    # Callback to calculate single pathway score

    @app.callback(
        [Output('cancersea-calculation-status', 'children'),
         Output('cancersea-feature-plot', 'figure'),
         Output('cancersea-violin-plot', 'figure'),
         Output('cancersea-heatmap', 'figure'),
         Output('cancersea-stats-table', 'children'),
         Output('cancersea-store', 'data')],
        Input('calculate-cancersea-button', 'n_clicks'),
        [State('cancersea-pathway-select', 'value'),
         State('adata-store', 'data'),
         State('cancersea-store', 'data')],
        prevent_initial_call=True
    )
    def calculate_single_pathway(n_clicks, pathway, adata_metadata, cancersea_data):
        """Calculate score for a single pathway."""
        if n_clicks is None or pathway is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get pathway genes from config
            pathway_genes = AppConfig.CANCERSEA_PATHWAYS.get(pathway, [])

            if len(pathway_genes) == 0:
                status = dbc.Alert(
                    f"No genes found for pathway: {pathway}",
                    color="warning"
                )
                return status, {}, {}, {}, "", cancersea_data

            # Calculate score
            score_name = calculate_pathway_score(
                adata,
                pathway,
                pathway_genes,
                ctrl_size=50
            )

            # Update adata in storage
            data_store.update_adata(session_id, adata)

            # Initialize or update cancersea_data
            if cancersea_data is None:
                cancersea_data = {'scores': {}}

            cancersea_data['scores'][pathway] = score_name

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Create visualizations
            feature_fig = create_pathway_feature_plot(
                adata, score_name, pathway)
            violin_fig = create_pathway_violin_plot(
                adata, score_name, pathway, cluster_key)

            # Create heatmap if multiple pathways calculated
            if len(cancersea_data['scores']) >= 2:
                heatmap_fig = create_pathway_heatmap(
                    adata, cancersea_data['scores'], cluster_key)
            else:
                heatmap_fig = go.Figure()
                heatmap_fig.add_annotation(
                    text="Calculate at least 2 pathways to show heatmap",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=14, color="gray")
                )

            # Get statistics
            stats = get_pathway_statistics(adata, score_name, cluster_key)
            stats_table = dbc.Table.from_dataframe(
                stats.round(3),
                striped=True,
                bordered=True,
                hover=True,
                size='sm'
            )

            # Status message
            status = dbc.Alert([
                html.I(className="fas fa-check-circle me-2"),
                f"Successfully calculated {pathway} score!"
            ], color="success")

            return status, feature_fig, violin_fig, heatmap_fig, stats_table, cancersea_data

        except Exception as e:
            print(f"Error calculating pathway: {e}")
            traceback.print_exc()

            status = dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error: {str(e)}"
            ], color="danger")

            return status, {}, {}, {}, "", cancersea_data

    # Callback to reload stored pathway scores when tab is opened

    @app.callback(
        [Output('cancersea-feature-plot', 'figure', allow_duplicate=True),
         Output('cancersea-violin-plot', 'figure', allow_duplicate=True),
         Output('cancersea-heatmap', 'figure', allow_duplicate=True),
         Output('cancersea-stats-table', 'children', allow_duplicate=True),
         Output('cancersea-reload-interval', 'disabled')],
        Input('cancersea-reload-interval', 'n_intervals'),
        [State('cancersea-store', 'data'),
         State('adata-store', 'data'),
         State('cancersea-pathway-select', 'value')],
        prevent_initial_call=True
    )
    def reload_cancersea_on_mount(n, cancersea_data, adata_metadata, pathway):
        """Reload pathway visualizations when tab is opened."""
        print(f"DEBUG: CancerSEA interval fired, n={n}")

        if cancersea_data is None or adata_metadata is None:
            print("DEBUG: No CancerSEA data to reload")
            raise PreventUpdate

        if 'scores' not in cancersea_data or len(cancersea_data['scores']) == 0:
            print("DEBUG: No pathway scores calculated yet")
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            print("DEBUG: Reloading CancerSEA visualizations")

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Get current pathway score
            if pathway in cancersea_data['scores']:
                score_name = cancersea_data['scores'][pathway]

                # Create visualizations
                feature_fig = create_pathway_feature_plot(
                    adata, score_name, pathway)
                violin_fig = create_pathway_violin_plot(
                    adata, score_name, pathway, cluster_key)

                # Get statistics
                stats = get_pathway_statistics(adata, score_name, cluster_key)
                stats_table = dbc.Table.from_dataframe(
                    stats.round(3),
                    striped=True,
                    bordered=True,
                    hover=True,
                    size='sm'
                )
            else:
                feature_fig = go.Figure()
                violin_fig = go.Figure()
                stats_table = ""

            # Create heatmap if multiple pathways
            if len(cancersea_data['scores']) >= 2:
                heatmap_fig = create_pathway_heatmap(
                    adata, cancersea_data['scores'], cluster_key)
            else:
                heatmap_fig = go.Figure()

            print("DEBUG: Successfully reloaded CancerSEA visualizations")

            return feature_fig, violin_fig, heatmap_fig, stats_table, True

        except Exception as e:
            print(f"ERROR reloading CancerSEA: {e}")
            traceback.print_exc()
            raise PreventUpdate

    # Callback to update visualizations when pathway selection changes

    @app.callback(
        [Output('cancersea-feature-plot', 'figure', allow_duplicate=True),
         Output('cancersea-violin-plot', 'figure', allow_duplicate=True),
         Output('cancersea-stats-table', 'children', allow_duplicate=True)],
        Input('cancersea-pathway-select', 'value'),
        [State('cancersea-store', 'data'),
         State('adata-store', 'data')],
        prevent_initial_call=True
    )
    def update_pathway_display(pathway, cancersea_data, adata_metadata):
        """Update display when pathway selection changes."""
        if pathway is None or cancersea_data is None or adata_metadata is None:
            raise PreventUpdate

        if 'scores' not in cancersea_data or pathway not in cancersea_data['scores']:
            # Pathway not yet calculated
            empty_fig = go.Figure()
            empty_fig.add_annotation(
                text=f"Click 'Calculate/Update Score' to score {pathway}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray")
            )
            return empty_fig, empty_fig, ""

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            score_name = cancersea_data['scores'][pathway]
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Create visualizations
            feature_fig = create_pathway_feature_plot(
                adata, score_name, pathway)
            violin_fig = create_pathway_violin_plot(
                adata, score_name, pathway, cluster_key)

            # Get statistics
            stats = get_pathway_statistics(adata, score_name, cluster_key)
            stats_table = dbc.Table.from_dataframe(
                stats.round(3),
                striped=True,
                bordered=True,
                hover=True,
                size='sm'
            )

            return feature_fig, violin_fig, stats_table

        except Exception as e:
            print(f"Error updating pathway display: {e}")
            raise PreventUpdate

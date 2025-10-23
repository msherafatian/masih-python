"""
Cell Cycle Module for MASIH application.

This module provides cell cycle phase analysis and visualization.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go

from masih.utils.data_store import get_adata
from masih.utils.cellcycle_plotting import (
    create_cellcycle_umap,
    create_cycle_score_plot,
    create_phase_by_cluster_plot
)
from masih.analysis.cellcycle_utils import get_cell_cycle_stats


def create_cellcycle_layout():
    """
    Create the layout for the Cell Cycle module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval for reload
        dcc.Interval(
            id='cellcycle-reload-interval',
            interval=500,
            max_intervals=1,
            disabled=False
        ),

        dbc.Row([
            dbc.Col([
                html.H3("Cell Cycle Analysis", className="mb-3"),
                html.P(
                    "Visualize cell cycle phase distribution and scores across your dataset.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Status alert
        html.Div(id='cellcycle-status-alert'),

        # Top row - Phase distribution and scores
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Cell Cycle Phase Distribution", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-cellcycle-plot",
                            children=dcc.Graph(
                                id='cellcycle-phase-plot',
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
                        html.H5("Cell Cycle Scores", className="mb-0")),
                    dbc.CardBody([
                        dbc.RadioItems(
                            id='cellcycle-score-select',
                            options=[
                                {'label': 'S Phase Score', 'value': 'S_score'},
                                {'label': 'G2M Phase Score', 'value': 'G2M_score'}
                            ],
                            value='S_score',
                            inline=True,
                            className="mb-3"
                        ),
                        dcc.Loading(
                            id="loading-cellcycle-score",
                            children=dcc.Graph(
                                id='cellcycle-score-plot',
                                config={'displayModeBar': True},
                                style={'height': '450px'}
                            )
                        )
                    ])
                ])
            ], md=6)
        ], className="mb-4"),

        # Bottom row - Phase by cluster
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Phase Distribution by Cluster", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            id="loading-phase-by-cluster",
                            children=dcc.Graph(
                                id='phase-by-cluster-plot',
                                config={'displayModeBar': True},
                                style={'height': '400px'}
                            )
                        )
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Statistics table
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Cell Cycle Statistics", className="mb-0")),
                    dbc.CardBody([
                        html.Div(id='cellcycle-stats-table')
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


def register_cellcycle_callbacks(app):
    """
    Register all callbacks for the cell cycle module.

    Args:
        app: Dash application instance
    """

    import traceback

    # Helper function to check if cell cycle data is available

    def has_cellcycle_data(adata):
        """Check if cell cycle scores are present."""
        required_cols = ['S_score', 'G2M_score', 'phase']
        return all(col in adata.obs.columns for col in required_cols)

    # Callback to update all plots

    @app.callback(
        [Output('cellcycle-phase-plot', 'figure'),
         Output('cellcycle-score-plot', 'figure'),
         Output('phase-by-cluster-plot', 'figure'),
         Output('cellcycle-stats-table', 'children'),
         Output('cellcycle-status-alert', 'children')],
        [Input('cellcycle-reload-interval', 'n_intervals'),
         Input('cellcycle-score-select', 'value')],
        State('adata-store', 'data'),
        prevent_initial_call=False
    )
    def update_cellcycle_plots(n, score_type, adata_metadata):
        """Update all cell cycle plots."""
        print(f"DEBUG: Cell cycle callback fired, n={n}")

        if adata_metadata is None:
            empty_fig = go.Figure()
            empty_fig.add_annotation(
                text="No data loaded",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray")
            )
            return empty_fig, empty_fig, empty_fig, "", ""

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Check if cell cycle data is available
            if not has_cellcycle_data(adata):
                print("DEBUG: Cell cycle data not available")

                # Create empty figures with message
                empty_fig = go.Figure()
                empty_fig.add_annotation(
                    text="Cell cycle scores not calculated.\nPlease enable cell cycle scoring during data processing.",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=14, color="gray")
                )

                # Warning alert
                alert = dbc.Alert([
                    html.I(className="fas fa-exclamation-triangle me-2"),
                    "Cell cycle scores not found. Please run cell cycle scoring during data processing."
                ], color="warning", className="mb-4")

                return empty_fig, empty_fig, empty_fig, "", alert

            print("DEBUG: Creating cell cycle plots")

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Create plots
            phase_fig = create_cellcycle_umap(adata)
            score_fig = create_cycle_score_plot(adata, score_type)
            cluster_fig = create_phase_by_cluster_plot(adata, cluster_key)

            # Get statistics
            stats = get_cell_cycle_stats(adata, cluster_key)
            stats_table = dbc.Table.from_dataframe(
                stats.round(2),
                striped=True,
                bordered=True,
                hover=True,
                size='sm'
            )

            # Success alert
            alert = dbc.Alert([
                html.I(className="fas fa-check-circle me-2"),
                f"Showing cell cycle analysis for {adata.n_obs:,} cells"
            ], color="success", dismissable=True, className="mb-4")

            print("DEBUG: Cell cycle plots created successfully")

            return phase_fig, score_fig, cluster_fig, stats_table, alert

        except Exception as e:
            print(f"ERROR in cell cycle callback: {e}")
            traceback.print_exc()

            empty_fig = go.Figure()
            empty_fig.add_annotation(
                text=f"Error: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="red")
            )

            alert = dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error loading cell cycle data: {str(e)}"
            ], color="danger", className="mb-4")

            return empty_fig, empty_fig, empty_fig, "", alert

    # Callback to disable interval after first load

    @app.callback(
        Output('cellcycle-reload-interval', 'disabled'),
        Input('cellcycle-reload-interval', 'n_intervals')
    )
    def disable_cellcycle_interval(n):
        """Disable interval after first fire."""
        return True

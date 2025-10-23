"""
Comparative Pathway Analysis Module for MASIH application.

This module provides functionality to compare multiple CancerSEA pathways
across clusters and calculate correlations.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, dash_table
from dash.exceptions import PreventUpdate
import pandas as pd
import plotly.graph_objects as go

from masih.config.settings import AppConfig
from masih.utils.data_store import get_adata
from masih.analysis.compare_utils import (
    calculate_pathway_correlation,
    calculate_pathway_by_cluster,
    calculate_pathway_statistics,
    create_pathway_profile
)
from masih.utils.compare_plotting import (
    create_correlation_matrix_plot,
    create_pathway_by_cluster_plot,
    create_comparison_heatmap
)


def create_compare_layout():
    """
    Create the layout for the Comparative Analysis module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval for reload
        dcc.Interval(
            id='compare-reload-interval',
            interval=500,
            max_intervals=1,
            disabled=False
        ),

        dbc.Row([
            dbc.Col([
                html.H3("Comparative Pathway Analysis", className="mb-3"),
                html.P(
                    "Compare multiple CancerSEA pathways to identify correlations and patterns.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Setup Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Pathway Comparison Setup", className="mb-0")),
                    dbc.CardBody([
                        dbc.Label(
                            "Select pathways to compare (2-5 pathways):"),
                        dcc.Dropdown(
                            id='compare-pathways-select',
                            options=[],
                            multi=True,
                            placeholder="Select pathways...",
                            className="mb-3"
                        ),

                        dbc.Button(
                            "Run Comparison",
                            id='run-comparison-button',
                            color='success',
                            className="w-100"
                        ),

                        html.Div(id='compare-status', className="mt-3")
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Visualizations
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Correlation Matrix", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            children=dcc.Graph(
                                id='compare-correlation-matrix',
                                config={'displayModeBar': True},
                                # Added width constraint
                                style={'height': '500px', 'width': '100%'}
                            )
                        )
                    ], style={'overflow': 'hidden'})  # Prevent overflow
                ])
            ], md=6),

            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Pathway Expression by Cluster", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            children=dcc.Graph(
                                id='compare-pathway-by-cluster',
                                config={'displayModeBar': True},
                                # Added width constraint
                                style={'height': '500px', 'width': '100%'}
                            )
                        )
                    ], style={'overflow': 'hidden'})  # Prevent overflow
                ])
            ], md=6)
        ], className="mb-4"),

        # Heatmap
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Pathway Comparison Heatmap", className="mb-0")),
                    dbc.CardBody([
                        dcc.Loading(
                            children=dcc.Graph(
                                id='compare-heatmap',
                                config={'displayModeBar': True},
                                # Let height be dynamic
                                style={'width': '100%'}
                            )
                        )
                    ], style={'overflow': 'auto'})  # Allow scrolling if needed
                ])
            ], width=12)
        ], className="mb-4"),

        # Statistics Table
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Statistical Summary", className="mb-0")),
                    dbc.CardBody([
                        html.Div(id='compare-stats-table')
                    ])
                ])
            ], width=12)
        ])
    ], fluid=True)

    return layout


def register_compare_callbacks(app):
    """
    Register all callbacks for the comparative analysis module.

    Args:
        app: Dash application instance
    """

    import traceback

    # Callback to update pathway options from calculated CancerSEA scores

    @app.callback(
        Output('compare-pathways-select', 'options'),
        Input('cancersea-store', 'data')
    )
    def update_pathway_options(cancersea_data):
        """Update pathway options based on calculated CancerSEA scores."""
        if cancersea_data is None or 'scores' not in cancersea_data:
            # Return all available pathways even if none calculated yet
            all_pathways = list(AppConfig.CANCERSEA_PATHWAYS.keys())
            return [{'label': p, 'value': p} for p in all_pathways]

        # Show calculated pathways at the top
        calculated = list(cancersea_data['scores'].keys())
        all_pathways = list(AppConfig.CANCERSEA_PATHWAYS.keys())

        options = []

        # Add calculated pathways first (with indicator)
        for pathway in calculated:
            options.append({
                'label': f"✓ {pathway}",
                'value': pathway
            })

        # Add remaining pathways
        for pathway in all_pathways:
            if pathway not in calculated:
                options.append({
                    'label': pathway,
                    'value': pathway
                })

        return options

    # Main callback to run comparison
    @app.callback(
        [Output('compare-status', 'children'),
         Output('compare-store', 'data')],
        Input('run-comparison-button', 'n_clicks'),
        [State('compare-pathways-select', 'value'),
         State('adata-store', 'data'),
         State('cancersea-store', 'data')],
        prevent_initial_call=True
    )
    def run_comparison(n_clicks, selected_pathways, adata_metadata, cancersea_data):
        """Run comparative pathway analysis."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        if not selected_pathways or len(selected_pathways) < 2:
            return dbc.Alert(
                "Please select at least 2 pathways for comparison",
                color="warning"
            ), None

        if len(selected_pathways) > 5:
            return dbc.Alert(
                "Please select no more than 5 pathways for optimal visualization",
                color="warning"
            ), None

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get existing CancerSEA scores
            cancersea_scores = cancersea_data.get(
                'scores', {}) if cancersea_data else {}

            print(f"DEBUG: Selected pathways: {selected_pathways}")
            print(
                f"DEBUG: Existing CancerSEA scores: {list(cancersea_scores.keys())}")

            # Check which pathways are already calculated
            available_scores = {}
            pathways_to_calculate = []

            for pathway in selected_pathways:
                if pathway in cancersea_scores:
                    score_col = cancersea_scores[pathway]
                    # Verify the score column actually exists in adata
                    if score_col in adata.obs.columns:
                        available_scores[pathway] = score_col
                        print(
                            f"DEBUG: Found existing score for {pathway}: {score_col}")
                    else:
                        print(
                            f"DEBUG: Score column {score_col} not found in adata, will recalculate {pathway}")
                        pathways_to_calculate.append(pathway)
                else:
                    pathways_to_calculate.append(pathway)

            print(f"DEBUG: Available scores: {available_scores}")
            print(f"DEBUG: Pathways to calculate: {pathways_to_calculate}")

            # Calculate missing pathways
            if pathways_to_calculate:
                from masih.utils.cancersea_utils import calculate_pathway_score
                from masih.utils.data_store import data_store

                for pathway in pathways_to_calculate:
                    try:
                        print(f"Calculating {pathway} for comparison...")
                        pathway_genes = AppConfig.CANCERSEA_PATHWAYS.get(
                            pathway, [])

                        if pathway_genes:
                            score_name = calculate_pathway_score(
                                adata,
                                pathway,
                                pathway_genes,
                                ctrl_size=50
                            )

                            # Verify score was added
                            if score_name in adata.obs.columns:
                                available_scores[pathway] = score_name
                                cancersea_scores[pathway] = score_name
                                print(
                                    f"✓ Calculated {pathway}, score column: {score_name}")
                            else:
                                print(
                                    f"ERROR: Score column {score_name} not found after calculation")
                        else:
                            print(f"ERROR: No genes found for {pathway}")

                    except Exception as e:
                        print(f"Warning: Could not calculate {pathway}: {e}")
                        traceback.print_exc()

                # Update adata in storage
                data_store.update_adata(session_id, adata)

                # Update cancersea store
                if cancersea_data is None:
                    cancersea_data = {'scores': {}}
                cancersea_data['scores'] = cancersea_scores
                data_store.store_additional_data(
                    session_id, 'cancersea_scores', cancersea_data)

            # Final check - make sure we have at least 2 pathways
            if len(available_scores) < 2:
                return dbc.Alert(
                    f"Only {len(available_scores)} pathway(s) available. Need at least 2 for comparison. Available: {list(available_scores.keys())}",
                    color="danger"
                ), None

            print(f"DEBUG: Final available scores: {available_scores}")

            # Store comparison results
            comparison_data = {
                'pathways': list(available_scores.keys()),
                'scores': available_scores,
                'session_id': session_id  # Store session ID for later retrieval
            }

            status = dbc.Alert([
                html.I(className="fas fa-check-circle me-2"),
                f"Comparison complete for {len(available_scores)} pathways: {', '.join(available_scores.keys())}"
            ], color="success")

            return status, comparison_data

        except Exception as e:
            print(f"Error in comparison: {e}")
            traceback.print_exc()

            status = dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error: {str(e)}"
            ], color="danger")

            return status, None

# Callback to update all visualizations
    @app.callback(
        [Output('compare-correlation-matrix', 'figure'),
         Output('compare-pathway-by-cluster', 'figure'),
         Output('compare-heatmap', 'figure'),
         Output('compare-stats-table', 'children')],
        [Input('compare-store', 'data'),
         Input('compare-reload-interval', 'n_intervals')],
        State('adata-store', 'data'),
        prevent_initial_call=False
    )
    def update_comparison_visualizations(comparison_data, n, adata_metadata):
        """Update all comparison visualizations."""
        print(f"DEBUG: Compare viz callback fired, n={n}")

        # Empty figures
        empty_fig = go.Figure()
        empty_fig.add_annotation(
            text="No comparison data available. Select pathways and click 'Run Comparison'.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=12, color="gray")
        )

        if comparison_data is None:
            print("DEBUG: No comparison data")
            return empty_fig, empty_fig, empty_fig, ""

        try:
            # Use session_id from comparison_data if available, otherwise from adata_metadata
            session_id = comparison_data.get('session_id') or (
                adata_metadata.get('session_id') if adata_metadata else None)

            if session_id is None:
                print("DEBUG: No session ID available")
                return empty_fig, empty_fig, empty_fig, ""

            adata = get_adata(session_id)

            if adata is None:
                print("DEBUG: Could not load adata")
                raise PreventUpdate

            score_columns = comparison_data['scores']
            print(f"DEBUG: Score columns: {score_columns}")

            # Verify all score columns exist
            missing_cols = [
                col for col in score_columns.values() if col not in adata.obs.columns]
            if missing_cols:
                print(f"ERROR: Missing score columns: {missing_cols}")
                print(f"Available columns: {adata.obs.columns.tolist()}")

            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Calculate correlation matrix
            print("Calculating correlations...")
            corr_results = calculate_pathway_correlation(adata, score_columns)
            corr_fig = create_correlation_matrix_plot(
                corr_results['correlation'],
                corr_results['p_values']
            )

            # Calculate pathway by cluster
            print("Calculating pathway by cluster...")
            pathway_by_cluster = calculate_pathway_by_cluster(
                adata, score_columns, cluster_key)
            print(f"DEBUG: Pathway by cluster data:\n{pathway_by_cluster}")
            cluster_fig = create_pathway_by_cluster_plot(pathway_by_cluster)

            # Create comparison heatmap
            print("Creating comparison heatmap...")
            heatmap_fig = create_comparison_heatmap(
                adata, score_columns, cluster_key)

            # Calculate statistics
            print("Calculating statistics...")
            stats = calculate_pathway_statistics(adata, score_columns)
            stats_table = dbc.Table.from_dataframe(
                stats.round(3),
                striped=True,
                bordered=True,
                hover=True,
                size='sm'
            )

            print("DEBUG: Comparison visualizations updated successfully")

            return corr_fig, cluster_fig, heatmap_fig, stats_table

        except Exception as e:
            print(f"Error updating comparison visualizations: {e}")
            traceback.print_exc()

            error_fig = go.Figure()
            error_fig.add_annotation(
                text=f"Error: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=12, color="red")
            )

            return error_fig, error_fig, error_fig, html.P(f"Error: {str(e)}", className="text-danger")

    # Callback to disable reload interval

    @app.callback(
        Output('compare-reload-interval', 'disabled'),
        Input('compare-reload-interval', 'n_intervals')
    )
    def disable_compare_interval(n):
        """Disable interval after first fire."""
        return True

"""
Trajectory Analysis Module for MASIH application.

This module provides trajectory inference using PAGA and diffusion pseudotime.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context
from dash.exceptions import PreventUpdate
import pandas as pd
import plotly.graph_objects as go

from masih.utils.data_store import get_adata, data_store
from masih.analysis.trajectory_utils import (
    prepare_trajectory_subset,
    run_trajectory_analysis,
    calculate_trajectory_stats,
    get_pseudotime_by_cluster,
    export_trajectory_data
)
from masih.utils.trajectory_plotting import (
    create_paga_trajectory_plot,
    create_pseudotime_plot,
    create_stratified_trajectory_plot
)


def create_trajectory_layout():
    """
    Create the layout for the Trajectory Analysis module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval for reload
        dcc.Interval(
            id='trajectory-reload-interval',
            interval=500,
            max_intervals=1,
            disabled=False
        ),

        dbc.Row([
            dbc.Col([
                html.H3("Trajectory Analysis", className="mb-3"),
                html.P(
                    "Infer developmental trajectories and pseudotime using PAGA.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Settings Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Trajectory Analysis Settings", className="mb-0")),
                    dbc.CardBody([
                        html.H6("Select Cell Types for Trajectory",
                                className="mb-2"),
                        html.P("Choose which cell types to include in the trajectory analysis.",
                               className="text-muted mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Annotation column:"),
                                dcc.Dropdown(
                                    id='trajectory-annotation-col',
                                    options=[],
                                    value='leiden'
                                )
                            ], md=6),

                            dbc.Col([
                                dbc.Label("Select cell types:"),
                                dcc.Dropdown(
                                    id='trajectory-celltypes',
                                    options=[],
                                    multi=True,
                                    placeholder="Select at least 2 cell types..."
                                )
                            ], md=6)
                        ], className="mb-3"),

                        html.Hr(),

                        html.H6("Trajectory Parameters", className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Number of neighbors:"),
                                dbc.Input(
                                    id='trajectory-n-neighbors',
                                    type='number',
                                    value=15,
                                    min=5,
                                    max=50,
                                    step=5
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Starting cluster (optional):"),
                                dcc.Dropdown(
                                    id='trajectory-start-cluster',
                                    options=[
                                        {'label': 'Auto-detect', 'value': 'auto'}],
                                    value='auto'
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Ending cluster (optional):"),
                                dcc.Dropdown(
                                    id='trajectory-end-cluster',
                                    options=[
                                        {'label': 'Auto-detect', 'value': 'auto'}],
                                    value='auto'
                                )
                            ], md=4)
                        ], className="mb-3"),

                        dbc.Alert([
                            html.I(className="fas fa-lightbulb me-2"),
                            "Tip: Leave start/end as 'Auto-detect' to let PAGA determine the trajectory automatically."
                        ], color="info", className="mb-3"),

                        html.Div(id='trajectory-filter-preview',
                                 className="mb-3"),

                        html.Hr(),

                        dbc.Button(
                            "Run Trajectory Analysis",
                            id='run-trajectory-button',
                            color='success',
                            size='lg',
                            className="w-100"
                        ),

                        html.Div(id='trajectory-status', className="mt-3")
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Visualization Cards
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        dbc.Tabs([
                            dbc.Tab(label="Main Trajectory", tab_id="main"),
                            dbc.Tab(label="Stratified Trajectory",
                                    tab_id="stratified"),
                            dbc.Tab(label="Pseudotime", tab_id="pseudotime")
                        ], id="trajectory-viz-tabs", active_tab="main")
                    ),
                    dbc.CardBody([
                        html.Div(id='trajectory-viz-content')
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Statistics
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Trajectory Statistics", className="mb-0")),
                    dbc.CardBody([
                        html.Div(id='trajectory-stats-display')
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


def register_trajectory_callbacks(app):
    """
    Register all callbacks for the trajectory analysis module.

    Args:
        app: Dash application instance
    """

    import traceback

    # Callback to update annotation column options

    @app.callback(
        [Output('trajectory-annotation-col', 'options'),
         Output('trajectory-annotation-col', 'value')],
        Input('adata-store', 'data')
    )
    def update_annotation_options(adata_metadata):
        """Update available annotation columns."""
        if adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get available columns
            meta_cols = list(adata.obs.columns)

            # Prioritize useful columns
            priority_cols = ['leiden', 'louvain',
                             'seurat_clusters', 'annotate', 'cell_type']
            available_priority = [c for c in priority_cols if c in meta_cols]
            other_cols = [c for c in meta_cols if c not in priority_cols]

            all_cols = available_priority + other_cols

            options = [{'label': col, 'value': col} for col in all_cols]

            # Default value
            default = 'leiden' if 'leiden' in meta_cols else all_cols[0]

            return options, default

        except Exception as e:
            print(f"Error updating annotation options: {e}")
            raise PreventUpdate

    # Callback to update cell type options based on annotation column

    @app.callback(
        [Output('trajectory-celltypes', 'options'),
         Output('trajectory-start-cluster', 'options'),
         Output('trajectory-end-cluster', 'options')],
        Input('trajectory-annotation-col', 'value'),
        State('adata-store', 'data')
    )
    def update_celltype_options(annotation_col, adata_metadata):
        """Update cell type options."""
        if annotation_col is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None or annotation_col not in adata.obs.columns:
                raise PreventUpdate

            # Get unique cell types
            cell_types = sorted(adata.obs[annotation_col].unique())

            celltype_options = [
                {'label': str(ct), 'value': str(ct)} for ct in cell_types]

            # Start/end options include auto-detect
            cluster_options = [{'label': 'Auto-detect',
                                'value': 'auto'}] + celltype_options

            return celltype_options, cluster_options, cluster_options

        except Exception as e:
            print(f"Error updating celltype options: {e}")
            raise PreventUpdate

    # Callback to show filter preview

    @app.callback(
        Output('trajectory-filter-preview', 'children'),
        [Input('trajectory-celltypes', 'value'),
         Input('trajectory-annotation-col', 'value')],
        State('adata-store', 'data')
    )
    def update_filter_preview(selected_types, annotation_col, adata_metadata):
        """Show preview of selected filters."""
        if not selected_types or adata_metadata is None:
            return html.P("No cell types selected", className="text-muted")

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Count cells
            mask = adata.obs[annotation_col].isin(selected_types)
            n_selected = mask.sum()
            n_total = adata.n_obs

            return dbc.Alert([
                html.Strong(f"Selected: {n_selected:,} cells "),
                f"out of {n_total:,} total cells",
                html.Br(),
                html.Small(f"Cell types: {', '.join(selected_types)}")
            ], color="light", className="mb-0")

        except Exception as e:
            return html.P(f"Error: {str(e)}", className="text-danger")

    # Main callback to run trajectory analysis

 # Main callback to run trajectory analysis
    @app.callback(
        [Output('trajectory-status', 'children'),
         Output('trajectory-store', 'data')],
        Input('run-trajectory-button', 'n_clicks'),
        [State('adata-store', 'data'),
         State('trajectory-celltypes', 'value'),
         State('trajectory-annotation-col', 'value'),
         State('trajectory-n-neighbors', 'value'),
         State('trajectory-start-cluster', 'value'),
         State('trajectory-end-cluster', 'value')],
        prevent_initial_call=True
    )
    def run_trajectory(n_clicks, adata_metadata, cell_types, annotation_col,
                       n_neighbors, start_cluster, end_cluster):
        """Run trajectory analysis."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        if not cell_types or len(cell_types) < 2:
            return dbc.Alert(
                "Please select at least 2 cell types for trajectory analysis",
                color="warning"
            ), None

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Prepare subset
            print("Preparing subset for trajectory...")
            adata_subset = prepare_trajectory_subset(
                adata,
                cell_types=cell_types,
                annotation_col=annotation_col
            )

            # Run trajectory analysis
            print("Running trajectory analysis...")
            start = None if start_cluster == 'auto' else start_cluster
            end = None if end_cluster == 'auto' else end_cluster

            trajectory_results = run_trajectory_analysis(
                adata_subset,
                cluster_key=annotation_col,
                n_neighbors=n_neighbors,
                root_cluster=start,
                end_cluster=end
            )

            # Calculate statistics
            stats = calculate_trajectory_stats(
                adata_subset,
                cluster_key=annotation_col,
                root_cluster=start,
                end_cluster=end
            )

            # Store the subset AnnData object with a new session ID
            import uuid
            subset_session_id = f"{session_id}_trajectory_{uuid.uuid4().hex[:8]}"
            data_store.store_adata(adata_subset, subset_session_id)

            # Store results (just the session ID string, not the whole object)
            trajectory_data = {
                'subset_session_id': subset_session_id,  # Changed from 'subset_session'
                'results': trajectory_results,
                'stats': stats,
                'parameters': {
                    'cell_types': cell_types,
                    'annotation_col': annotation_col,
                    'n_neighbors': n_neighbors,
                    'start_cluster': start,
                    'end_cluster': end
                }
            }

            # Also store in main data store
            data_store.store_additional_data(
                session_id, 'trajectory', trajectory_data)

            status = dbc.Alert([
                html.I(className="fas fa-check-circle me-2"),
                f"Trajectory analysis complete! Analyzed {stats['n_cells']:,} cells across {stats['n_clusters']} cell types."
            ], color="success")

            return status, trajectory_data

        except Exception as e:
            print(f"Error in trajectory analysis: {e}")
            traceback.print_exc()

            status = dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error: {str(e)}"
            ], color="danger")

            return status, None

    # Callback to update visualizations based on selected tab

    # Callback to update visualizations based on selected tab
    @app.callback(
        Output('trajectory-viz-content', 'children'),
        [Input('trajectory-viz-tabs', 'active_tab'),
         Input('trajectory-store', 'data')],
        prevent_initial_call=False
    )
    def update_trajectory_viz(active_tab, trajectory_data):
        """Update trajectory visualization based on selected tab."""
        if trajectory_data is None:
            return html.Div([
                html.H5(
                    "No trajectory calculated yet.",
                    className="text-muted text-center",
                    style={'padding': '50px'}
                ),
                html.P(
                    "Select cell types and click 'Run Trajectory Analysis' to start.",
                    className="text-muted text-center"
                )
            ])

        try:
            # Get subset AnnData using the session ID string
            # Changed key
            subset_session_id = trajectory_data['subset_session_id']
            adata_subset = get_adata(subset_session_id)

            if adata_subset is None:
                raise ValueError("Could not load trajectory data")

            params = trajectory_data['parameters']
            annotation_col = params['annotation_col']

            # Create appropriate plot based on tab
            if active_tab == "main":
                fig = create_paga_trajectory_plot(adata_subset, annotation_col)
            elif active_tab == "stratified":
                fig = create_stratified_trajectory_plot(
                    adata_subset, annotation_col)
            elif active_tab == "pseudotime":
                fig = create_pseudotime_plot(adata_subset)
            else:
                fig = go.Figure()

            return dcc.Loading(
                children=dcc.Graph(
                    figure=fig,
                    config={'displayModeBar': True},
                    style={'height': '600px'}
                )
            )

        except Exception as e:
            print(f"Error creating trajectory plot: {e}")
            traceback.print_exc()

            return dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error creating visualization: {str(e)}"
            ], color="danger")

    # Callback to show statistics

    @app.callback(
        Output('trajectory-stats-display', 'children'),
        Input('trajectory-store', 'data')
    )
    def update_trajectory_stats(trajectory_data):
        """Display trajectory statistics."""
        if trajectory_data is None:
            return html.P("No trajectory calculated", className="text-muted")

        try:
            stats = trajectory_data['stats']
            params = trajectory_data['parameters']

            stats_text = f"""
═══════════════════════════════════════
Trajectory Analysis Summary
═══════════════════════════════════════

Number of cells: {stats['n_cells']:,}
Number of clusters: {stats['n_clusters']}
Clusters included: {', '.join(map(str, stats['clusters']))}

Starting cluster: {params['start_cluster'] or 'Auto-detected'}
Ending cluster: {params['end_cluster'] or 'Auto-detected'}
"""

            if 'pseudotime_range' in stats:
                pt_min, pt_max = stats['pseudotime_range']
                stats_text += f"""
Pseudotime range: {pt_min:.3f} - {pt_max:.3f}
Mean pseudotime: {stats['pseudotime_mean']:.3f}
"""

            if 'n_connections' in stats:
                stats_text += f"""
PAGA connections: {stats['n_connections']}
Mean connectivity: {stats['mean_connectivity']:.3f}
"""

            stats_text += f"""
═══════════════════════════════════════
Method: PAGA + Diffusion Pseudotime
"""

            return html.Pre(stats_text, style={'fontFamily': 'monospace'})

        except Exception as e:
            return html.P(f"Error: {str(e)}", className="text-danger")

    # Callback to reload trajectory on tab open

    @app.callback(
        Output('trajectory-reload-interval', 'disabled'),
        Input('trajectory-reload-interval', 'n_intervals'),
        State('trajectory-store', 'data')
    )
    def reload_trajectory_on_mount(n, trajectory_data):
        """Reload trajectory when tab opens."""
        print(f"DEBUG: Trajectory interval fired, n={n}")

        if trajectory_data:
            print("DEBUG: Trajectory data available")

        return True  # Disable after first fire

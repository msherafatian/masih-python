"""
Export Module for MASIH application.

This module provides functionality to export plots, data tables,
and generate methods text for publications.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context
from dash.exceptions import PreventUpdate
import pandas as pd
import io
import zipfile
import traceback
import tempfile
import os
from datetime import datetime

from masih.config.settings import AppConfig
from masih.utils.data_store import get_adata, data_store
from masih.utils.export_utils import (
    calculate_cluster_stats,
    generate_methods_text,
    generate_citations,
    generate_figure_legend,
    print_analysis_summary,
    prepare_download_data
)
from masih.utils.plot_export import (
    create_export_plot,
    fig_to_base64,
    get_available_plot_types
)


def create_export_layout():
    """
    Create the layout for the Export module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval for reload
        dcc.Interval(
            id='export-reload-interval',
            interval=500,
            max_intervals=1,
            disabled=False
        ),

        dbc.Row([
            dbc.Col([
                html.H3("Export & Documentation", className="mb-3"),
                html.P(
                    "Export your analysis results, plots, and generate publication-ready documentation.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Plot Export Section
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Plot Export Settings", className="mb-0")),
                    dbc.CardBody([
                        html.H6("Export Configuration", className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Select plot to export:"),
                                dcc.Dropdown(
                                    id='export-plot-type',
                                    options=[
                                        {'label': 'Current Cluster Plot',
                                            'value': 'cluster'},
                                        {'label': 'Cell Cycle Plot',
                                            'value': 'cellcycle'},
                                        {'label': 'Cell Cycle Distribution',
                                            'value': 'cellcycle_distribution'},
                                        {'label': 'CancerSEA Feature Plot',
                                            'value': 'cancersea_feature'},
                                        {'label': 'CancerSEA Heatmap',
                                            'value': 'cancersea_heatmap'},
                                        {'label': 'Pathway Correlation Matrix',
                                            'value': 'correlation_matrix'},
                                        {'label': 'Cluster Tree', 'value': 'tree'},
                                        {'label': 'Custom Gene Expression',
                                            'value': 'gene'},
                                        {'label': 'Marker Heatmap',
                                            'value': 'marker_heatmap'}
                                    ],
                                    value='cluster'
                                )
                            ], md=12)
                        ], className="mb-3"),

                        # Conditional inputs
                        html.Div(id='export-plot-options'),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Format:"),
                                dcc.Dropdown(
                                    id='export-format',
                                    options=[
                                        {'label': 'PNG', 'value': 'png'},
                                        {'label': 'JPEG', 'value': 'jpeg'},
                                        {'label': 'SVG', 'value': 'svg'},
                                        {'label': 'PDF', 'value': 'pdf'}
                                    ],
                                    value='png',
                                    clearable=False
                                )
                            ], md=6),

                            dbc.Col([
                                dbc.Label("Resolution (DPI):"),
                                dcc.Dropdown(
                                    id='export-dpi',
                                    options=[
                                        {'label': 'Screen (72)', 'value': 1},
                                        {'label': 'Print (300)',
                                         'value': 4.17},
                                        {'label': 'Publication (600)',
                                         'value': 8.33}
                                    ],
                                    value=4.17,
                                    clearable=False
                                )
                            ], md=6)
                        ], className="mb-3"),

                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Width (inches):"),
                                dbc.Input(
                                    id='export-width',
                                    type='number',
                                    value=8,
                                    min=2,
                                    max=20,
                                    step=0.5
                                )
                            ], md=6),

                            dbc.Col([
                                dbc.Label("Height (inches):"),
                                dbc.Input(
                                    id='export-height',
                                    type='number',
                                    value=6,
                                    min=2,
                                    max=20,
                                    step=0.5
                                )
                            ], md=6)
                        ], className="mb-3"),

                        dbc.Checklist(
                            id='export-white-bg',
                            options=[
                                {'label': ' White background', 'value': 'white'}],
                            value=['white'],
                            inline=True,
                            className="mb-3"
                        ),

                        dbc.Input(
                            id='export-custom-title',
                            placeholder='Custom title (optional)',
                            type='text',
                            className="mb-3"
                        ),

                        dbc.Button(
                            "Download Plot",
                            id='download-plot-button',
                            color='success',
                            className="w-100 mb-3"
                        ),
                        dcc.Download(id='download-plot'),

                        html.Hr(),

                        html.H6("Batch Export", className="mb-2"),
                        html.P(
                            "Export all available plots with current settings:", className="text-muted"),
                        dbc.Button(
                            "Download All Plots (ZIP)",
                            id='download-all-plots-button',
                            color='warning',
                            className="w-100"
                        ),
                        dcc.Download(id='download-all-plots')
                    ])
                ])
            ], md=6),

            # Data Export Section
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(html.H5("Data Export", className="mb-0")),
                    dbc.CardBody([
                        html.H6("Export Data Tables", className="mb-3"),

                        dbc.Checklist(
                            id='export-data-options',
                            options=[
                                {'label': 'Cell metadata', 'value': 'metadata'},
                                {'label': 'Cluster statistics',
                                 'value': 'cluster_stats'},
                                {'label': 'CancerSEA scores',
                                 'value': 'cancersea_scores'},
                                {'label': 'Pathway averages per cluster',
                                 'value': 'pathway_avg'},
                                {'label': 'Cell cycle proportions',
                                 'value': 'cellcycle_props'},
                                {'label': 'Marker genes', 'value': 'markers'},
                                {'label': 'Trajectory summary',
                                 'value': 'trajectory'}
                            ],
                            value=['metadata', 'cluster_stats'],
                            className="mb-3"
                        ),

                        dbc.Button(
                            "Download Data (Excel)",
                            id='download-data-button',
                            color='primary',
                            className="w-100 mb-3"
                        ),
                        dcc.Download(id='download-data'),

                        html.Hr(),

                        html.H6("Export Full Object", className="mb-2"),
                        html.P(
                            "Download the complete AnnData object for use in Python/R:",
                            className="text-muted"
                        ),
                        dbc.Button(
                            "Download AnnData (.h5ad)",
                            id='download-anndata-button',
                            color='info',
                            className="w-100"
                        ),
                        dcc.Download(id='download-anndata')
                    ])
                ])
            ], md=6)
        ], className="mb-4"),

        # Documentation Section
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Documentation & Citations", className="mb-0")),
                    dbc.CardBody([
                        dbc.Tabs(
                            id="documentation-tabs",
                            active_tab="methods",
                            children=[
                                dbc.Tab(label="Methods", tab_id="methods"),
                                dbc.Tab(label="Citations", tab_id="citations"),
                                dbc.Tab(label="Figure Legend",
                                        tab_id="legend"),
                                dbc.Tab(label="Analysis Summary",
                                        tab_id="summary")
                            ]
                        ),
                        html.Div(id='documentation-content', className="mt-3")
                    ])
                ])
            ])
        ])
    ], fluid=True)

    return layout


def register_export_callbacks(app):
    """
    Register all callbacks for the Export module.

    Args:
        app: Dash application instance
    """

    # Callback for conditional plot options
    @app.callback(
        Output('export-plot-options', 'children'),
        Input('export-plot-type', 'value'),
        State('cancersea-store', 'data')
    )
    def update_plot_options(plot_type, cancersea_data):
        """Show additional options based on plot type."""
        if plot_type == 'gene':
            return dbc.Row([
                dbc.Col([
                    dbc.Label("Gene name:"),
                    dbc.Input(
                        id='export-gene-name',
                        type='text',
                        placeholder='Enter gene name',
                        className="mb-3"
                    )
                ])
            ])

        elif plot_type in ['cancersea_feature', 'cancersea_heatmap']:
            if cancersea_data and 'scores' in cancersea_data:
                pathways = list(cancersea_data['scores'].keys())
                if pathways:
                    return dbc.Row([
                        dbc.Col([
                            dbc.Label("Select pathway:"),
                            dcc.Dropdown(
                                id='export-pathway-name',
                                options=[{'label': p, 'value': p}
                                         for p in pathways],
                                value=pathways[0],
                                className="mb-3"
                            )
                        ])
                    ])

        return html.Div()

    # Callback to download single plot
    @app.callback(
        Output('download-plot', 'data'),
        Input('download-plot-button', 'n_clicks'),
        [State('export-plot-type', 'value'),
         State('export-format', 'value'),
         State('export-width', 'value'),
         State('export-height', 'value'),
         State('export-dpi', 'value'),
         State('export-white-bg', 'value'),
         State('export-custom-title', 'value'),
         State('adata-store', 'data'),
         State('cancersea-store', 'data'),
         State('markers-store', 'data')],
        prevent_initial_call=True
    )
    def download_plot(n_clicks, plot_type, format_type, width, height, dpi,
                      white_bg, custom_title, adata_metadata, cancersea_data, markers_data):
        """Download a single plot."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Get additional data in the format expected by create_export_plot
            cancersea_scores = cancersea_data.get(
                'scores', {}) if cancersea_data else {}

            markers_genes = None
            if markers_data:
                top_markers = pd.DataFrame(markers_data.get('top_markers', []))
                if not top_markers.empty:
                    markers_genes = top_markers['gene'].unique()[:50].tolist()

            # Create plot using the correct function signature
            fig = create_export_plot(
                adata,
                plot_type=plot_type,
                cluster_key=cluster_key,
                cancersea_scores=cancersea_scores,
                markers_genes=markers_genes,
                custom_title=custom_title,
                white_bg='white' in white_bg if white_bg else True
            )

            if fig is None:
                raise ValueError("Could not create plot")

            # Convert inches to pixels (96 DPI base)
            width_px = int(width * 96)
            height_px = int(height * 96)

            img_bytes = fig.to_image(
                format=format_type,
                width=width_px,
                height=height_px,
                scale=dpi
            )

            filename = f"masih_{plot_type}_{datetime.now().strftime('%Y%m%d')}.{format_type}"

            return dcc.send_bytes(img_bytes, filename)

        except Exception as e:
            print(f"Error downloading plot: {e}")
            traceback.print_exc()
            raise PreventUpdate

    # Callback to download all plots as ZIP
    @app.callback(
        Output('download-all-plots', 'data'),
        Input('download-all-plots-button', 'n_clicks'),
        [State('export-format', 'value'),
         State('export-width', 'value'),
         State('export-height', 'value'),
         State('export-dpi', 'value'),
         State('export-white-bg', 'value'),
         State('adata-store', 'data'),
         State('cancersea-store', 'data'),
         State('markers-store', 'data')],
        prevent_initial_call=True
    )
    def download_all_plots(n_clicks, format_type, width, height, dpi, white_bg,
                           adata_metadata, cancersea_data, markers_data):
        """Download all available plots as a ZIP file."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get cluster key
            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Get additional data
            cancersea_scores = cancersea_data.get(
                'scores', {}) if cancersea_data else {}

            markers_genes = None
            if markers_data:
                top_markers = pd.DataFrame(markers_data.get('top_markers', []))
                if not top_markers.empty:
                    markers_genes = top_markers['gene'].unique()[:50].tolist()

            # Get available plot types based on data
            plot_types = get_available_plot_types(
                adata, cancersea_data, markers_data)

            zip_buffer = io.BytesIO()

            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for plot_type in plot_types:
                    try:
                        # Create plot using correct signature
                        fig = create_export_plot(
                            adata,
                            plot_type=plot_type,
                            cluster_key=cluster_key,
                            cancersea_scores=cancersea_scores,
                            markers_genes=markers_genes,
                            white_bg='white' in white_bg if white_bg else True
                        )

                        if fig is not None:
                            # Convert inches to pixels (96 DPI base)
                            width_px = int(width * 96)
                            height_px = int(height * 96)

                            img_bytes = fig.to_image(
                                format=format_type,
                                width=width_px,
                                height=height_px,
                                scale=dpi
                            )

                            # Add to ZIP
                            filename = f"{plot_type}.{format_type}"
                            zip_file.writestr(filename, img_bytes)

                    except Exception as e:
                        print(f"Error creating plot {plot_type}: {e}")
                        continue

            zip_buffer.seek(0)
            filename = f"masih_all_plots_{datetime.now().strftime('%Y%m%d')}.zip"

            return dcc.send_bytes(zip_buffer.getvalue(), filename)

        except Exception as e:
            print(f"Error creating ZIP: {e}")
            traceback.print_exc()
            raise PreventUpdate

    # Callback to download data as Excel
    @app.callback(
        Output('download-data', 'data'),
        Input('download-data-button', 'n_clicks'),
        [State('export-data-options', 'value'),
         State('adata-store', 'data'),
         State('cancersea-store', 'data'),
         State('markers-store', 'data'),
         State('trajectory-store', 'data')],
        prevent_initial_call=True
    )
    def download_data(n_clicks, options, adata_metadata, cancersea_data, markers_data, trajectory_data):
        """Download data tables as Excel file."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'

            # Get additional data with improved handling
            cancersea_scores = cancersea_data.get(
                'scores', {}) if cancersea_data else None

            # Handle markers data - try multiple possible structures
            markers_df = None
            if markers_data:
                print(
                    f"DEBUG: markers_data keys: {markers_data.keys() if isinstance(markers_data, dict) else 'not a dict'}")

                # Try 'markers' key first
                if isinstance(markers_data, dict) and 'markers' in markers_data:
                    markers_list = markers_data.get('markers', [])
                    if markers_list:
                        markers_df = pd.DataFrame(markers_list)
                        print(
                            f"DEBUG: Created markers_df from 'markers' key, shape: {markers_df.shape}")

                # Try 'top_markers' key as fallback
                elif isinstance(markers_data, dict) and 'top_markers' in markers_data:
                    top_markers = markers_data.get('top_markers', [])
                    if top_markers:
                        markers_df = pd.DataFrame(top_markers)
                        print(
                            f"DEBUG: Created markers_df from 'top_markers' key, shape: {markers_df.shape}")

                # If markers_data is already a list
                elif isinstance(markers_data, list) and markers_data:
                    markers_df = pd.DataFrame(markers_data)
                    print(
                        f"DEBUG: Created markers_df from list, shape: {markers_df.shape}")

            if markers_df is None or (isinstance(markers_df, pd.DataFrame) and markers_df.empty):
                print("WARNING: No markers data available for export")

            # Handle trajectory data with better checking
            traj_data = None
            if trajectory_data:
                print(f"DEBUG: trajectory_data type: {type(trajectory_data)}")
                if isinstance(trajectory_data, dict):
                    print(
                        f"DEBUG: trajectory_data keys: {trajectory_data.keys()}")
                    # Only pass if it has actual data
                    if trajectory_data:  # Non-empty dict
                        traj_data = trajectory_data
                        print(
                            "DEBUG: Passing trajectory data to prepare_download_data")

            # Debug: Show what options were selected
            print(f"DEBUG: Selected export options: {options}")

            # Prepare data sheets
            sheets = prepare_download_data(
                adata,
                include_options=options or [],
                cluster_key=cluster_key,
                cancersea_scores=cancersea_scores,
                markers_df=markers_df,
                trajectory_data=traj_data
            )

            # Debug: Show which sheets were created
            print(f"DEBUG: Created sheets: {list(sheets.keys())}")
            for sheet_name, df in sheets.items():
                print(
                    f"  - {sheet_name}: {df.shape[0]} rows, {df.shape[1]} columns")

            if not sheets:
                print("WARNING: No data sheets to export!")
                raise PreventUpdate

            # Create Excel file in memory
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='openpyxl') as writer:
                for sheet_name, df in sheets.items():
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

            output.seek(0)
            filename = f"masih_data_{datetime.now().strftime('%Y%m%d')}.xlsx"

            print(f"SUCCESS: Excel file created with {len(sheets)} sheets")
            return dcc.send_bytes(output.getvalue(), filename)

        except Exception as e:
            print(f"Error creating Excel file: {e}")
            traceback.print_exc()
            raise PreventUpdate

    # Callback to download AnnData object
    @app.callback(
        Output('download-anndata', 'data'),
        Input('download-anndata-button', 'n_clicks'),
        State('adata-store', 'data'),
        prevent_initial_call=True
    )
    def download_anndata(n_clicks, adata_metadata):
        """Download AnnData object."""
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Create a temporary file (Python 3.13 compatible)
            with tempfile.NamedTemporaryFile(mode='wb', suffix='.h5ad', delete=False) as tmp_file:
                tmp_path = tmp_file.name

            try:
                # Write to temporary file
                adata.write_h5ad(tmp_path)

                # Read the file content
                with open(tmp_path, 'rb') as f:
                    file_content = f.read()

                filename = f"masih_adata_{datetime.now().strftime('%Y%m%d')}.h5ad"

                return dcc.send_bytes(file_content, filename)

            finally:
                # Clean up temporary file
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)

        except Exception as e:
            print(f"Error saving AnnData: {e}")
            traceback.print_exc()
            raise PreventUpdate

    # Callback to update documentation content based on selected tab
    @app.callback(
        Output('documentation-content', 'children'),
        Input('documentation-tabs', 'active_tab'),
        [State('adata-store', 'data'),
         State('processed-flags', 'data'),
         State('cancersea-store', 'data'),
         State('markers-store', 'data'),
         State('trajectory-store', 'data')]
    )
    def update_documentation_content(active_tab, adata_metadata, flags,
                                     cancersea_data, markers_data, trajectory_data):
        """Update documentation content based on selected tab."""
        if adata_metadata is None:
            return html.P("No data loaded", className="text-muted")

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            if active_tab == "methods":
                # Generate methods text
                processing_params = {
                    'version': AppConfig.VERSION,
                    'qc_performed': flags.get('qc_done', False) if flags else False,
                    'min_genes': 200,
                    'max_genes': 2500,
                    'max_mt_percent': 5,
                    'target_sum': 10000,
                    'n_top_genes': 2000,
                    'n_pcs': 50,
                    'resolution': 0.5,
                    'n_neighbors': 15
                }

                cancersea_pathways = list(cancersea_data.get(
                    'scores', {}).keys()) if cancersea_data else []
                has_markers = markers_data is not None
                has_trajectory = trajectory_data is not None and bool(
                    trajectory_data)

                methods_text = generate_methods_text(
                    processing_params,
                    cancersea_pathways,
                    has_markers,
                    has_trajectory
                )

                return html.Div([
                    html.Pre(methods_text, style={'whiteSpace': 'pre-wrap'}),
                    dbc.Button(
                        [html.I(className="fas fa-copy me-2"),
                         "Copy to Clipboard"],
                        id='copy-methods-button',
                        color='primary',
                        size='sm',
                        className="mt-2"
                    )
                ])

            elif active_tab == "citations":
                citations = generate_citations()
                return html.Pre(citations, style={'whiteSpace': 'pre-wrap'})

            elif active_tab == "legend":
                cluster_key = 'leiden' if 'leiden' in adata.obs.columns else 'louvain'
                legend = generate_figure_legend(adata, cluster_key)
                return html.Pre(legend, style={'whiteSpace': 'pre-wrap'})

            elif active_tab == "summary":
                processing_params = {
                    'target_sum': 10000,
                    'n_top_genes': 2000,
                    'n_pcs': 50,
                    'resolution': 0.5,
                    'n_neighbors': 15
                }

                cancersea_scores = cancersea_data.get(
                    'scores', {}) if cancersea_data else {}
                has_markers = markers_data is not None
                has_trajectory = trajectory_data is not None and bool(
                    trajectory_data)

                summary = print_analysis_summary(
                    adata,
                    processing_params,
                    cancersea_scores,
                    has_markers,
                    has_trajectory
                )

                return html.Pre(summary, style={'whiteSpace': 'pre-wrap', 'fontFamily': 'monospace'})

        except Exception as e:
            print(f"Error generating documentation: {e}")
            traceback.print_exc()
            return html.P(f"Error: {str(e)}", className="text-danger")

    # Callback to reload documentation when tab opens
    @app.callback(
        Output('export-reload-interval', 'disabled'),
        Input('export-reload-interval', 'n_intervals')
    )
    def reload_export_on_mount(n):
        """Disable interval after first fire."""
        print(f"DEBUG: Export interval fired, n={n}")
        return True

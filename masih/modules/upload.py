"""
Upload & QC Module for MASIH application.

This module handles data upload, quality control, normalization,
dimensionality reduction, and clustering - matching the R/Shiny implementation.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly.graph_objects as go
import base64
import io

from masih.config.settings import AppConfig


def create_upload_layout():
    """
    Create the layout for the Upload & QC module.
    Matches the R mod_upload_ui layout.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H3("Data Upload & Processing", className="mb-3"),
                html.P(
                    "Upload your single-cell RNA-seq data and perform quality control and preprocessing.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Main Upload Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(html.H5("Data Upload", className="mb-0")),
                    dbc.CardBody([
                        # Input type selection
                        dbc.Label("Select input type:", className="fw-bold"),
                        dbc.RadioItems(
                            id='upload-input-type',
                            options=[
                                {'label': ' AnnData Object (.h5ad)',
                                 'value': 'h5ad'},
                                {'label': ' 10X HDF5 (.h5)', 'value': 'h5'},
                                {'label': ' 10X Directory (filtered_feature_bc_matrix)',
                                 'value': '10x_dir'}
                            ],
                            value='h5ad',
                            inline=False,
                            className="mb-3"
                        ),

                        html.Hr(),

                        # File upload (shown for h5ad and h5)
                        html.Div(id='upload-file-area', children=[
                            dcc.Upload(
                                id='upload-data-file',
                                children=html.Div([
                                    html.I(
                                        className="fas fa-cloud-upload-alt fa-3x mb-3"),
                                    html.H5(
                                        'Drag and Drop or Click to Select File'),
                                    html.P(
                                        'Supported formats: .h5ad, .h5',
                                        className="text-muted small"
                                    )
                                ]),
                                style={
                                    'width': '100%',
                                    'height': '150px',
                                    'lineHeight': '150px',
                                    'borderWidth': '2px',
                                    'borderStyle': 'dashed',
                                    'borderRadius': '10px',
                                    'textAlign': 'center',
                                    'cursor': 'pointer',
                                    'backgroundColor': '#fafafa'
                                },
                                multiple=False
                            )
                        ]),

                        # Upload status
                        html.Div(id='upload-status', className="mt-3"),

                        # Processing options (shown after upload)
                        html.Div(id='processing-options-area', style={'display': 'none'}, children=[
                            html.Hr(),
                            html.H5("Processing Options",
                                    className="mt-3 mb-3"),

                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("Number of variable features:"),
                                    dbc.Input(
                                        id='process-nfeatures',
                                        type='number',
                                        value=AppConfig.DEFAULT_NFEATURES,
                                        min=1000,
                                        max=5000,
                                        step=500
                                    )
                                ], md=3),

                                dbc.Col([
                                    dbc.Label("Number of PCs for clustering:"),
                                    dbc.Input(
                                        id='process-dims',
                                        type='number',
                                        value=AppConfig.DEFAULT_DIMS,
                                        min=10,
                                        max=50,
                                        step=5
                                    )
                                ], md=3),

                                dbc.Col([
                                    dbc.Label("Clustering resolution:"),
                                    dbc.Input(
                                        id='process-resolution',
                                        type='number',
                                        value=AppConfig.DEFAULT_RESOLUTION,
                                        min=0.1,
                                        max=2.0,
                                        step=0.1
                                    )
                                ], md=3),

                                dbc.Col([
                                    dbc.Label("Number of neighbors:"),
                                    dbc.Input(
                                        id='process-n-neighbors',
                                        type='number',
                                        value=AppConfig.DEFAULT_N_NEIGHBORS,
                                        min=5,
                                        max=100,
                                        step=5
                                    )
                                ], md=3)
                            ], className="mb-3"),

                            dbc.Row([
                                dbc.Col([
                                    dbc.Checklist(
                                        id='process-recluster',
                                        options=[
                                            {'label': ' Force re-clustering', 'value': 'recluster'}],
                                        value=[],
                                        inline=True
                                    )
                                ], md=4),

                                dbc.Col([
                                    dbc.Checklist(
                                        id='process-run-qc',
                                        options=[
                                            {'label': ' Run QC filtering', 'value': 'qc'}],
                                        value=[],
                                        inline=True
                                    )
                                ], md=4)
                            ], className="mb-3"),

                            # QC parameters (conditionally shown)
                            html.Div(id='qc-parameters-area', style={'display': 'none'}, children=[
                                html.H6("QC Parameters:",
                                        className="mt-2 mb-2"),
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Label("Max mitochondrial %:"),
                                        dbc.Input(
                                            id='qc-max-mt',
                                            type='number',
                                            value=AppConfig.DEFAULT_MAX_MT,
                                            min=1,
                                            max=100,
                                            step=1
                                        )
                                    ], md=4),

                                    dbc.Col([
                                        dbc.Label("Min nCount_RNA:"),
                                        dbc.Input(
                                            id='qc-min-ncount',
                                            type='number',
                                            value=AppConfig.DEFAULT_MIN_NCOUNT,
                                            min=0,
                                            max=10000,
                                            step=100
                                        )
                                    ], md=4),

                                    dbc.Col([
                                        dbc.Label("Max nCount_RNA:"),
                                        dbc.Input(
                                            id='qc-max-ncount',
                                            type='number',
                                            value=AppConfig.DEFAULT_MAX_NCOUNT,
                                            min=1000,
                                            max=100000,
                                            step=1000
                                        )
                                    ], md=4)
                                ])
                            ]),

                            dbc.Button(
                                "Process Data",
                                id='run-processing-button',
                                color='success',
                                className="mt-3",
                                size="lg"
                            ),

                            html.Hr(),
                            html.H5("Data Summary", className="mt-3"),
                            html.Pre(id='data-summary',
                                     className="p-3 bg-light border rounded")
                        ])
                    ])
                ])
            ], width=12)
        ], className="mb-4"),

        # Processing Log
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Processing Log", className="mb-0")),
                    dbc.CardBody([
                        html.Pre(
                            id='processing-log',
                            className="p-3 bg-dark text-light",
                            style={
                                'maxHeight': '400px',
                                'overflowY': 'auto',
                                'fontFamily': 'monospace',
                                'fontSize': '0.9rem'
                            }
                        )
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


# Callback registration function (to be called from callback_manager)
def register_upload_callbacks(app):
    """
    Register all callbacks for the upload module.

    Args:
        app: Dash application instance
    """
    import base64
    import io
    import traceback
    from datetime import datetime

    import anndata as ad
    import scanpy as sc

    from masih.utils.data_store import save_adata, get_adata, update_adata
    from masih.utils.scanpy_utils import (
        setup_scanpy_settings,
        check_adata_analyses,
        calculate_qc_metrics,
        filter_cells_qc,
        normalize_data,
        find_variable_genes,
        scale_data,
        run_pca,
        compute_neighbors,
        run_umap,
        run_tsne,
        find_clusters,
        score_cell_cycle,
        CELL_CYCLE_GENES
    )

    # Setup Scanpy settings
    setup_scanpy_settings()

    # Callback to show/hide QC parameters

    @app.callback(
        Output('qc-parameters-area', 'style'),
        Input('process-run-qc', 'value')
    )
    def toggle_qc_parameters(run_qc):
        """Show QC parameters when QC checkbox is checked."""
        if run_qc and 'qc' in run_qc:
            return {'display': 'block'}
        return {'display': 'none'}

    # Callback to handle file upload

    @app.callback(
        [Output('processing-options-area', 'style'),
         Output('upload-status', 'children'),
         Output('adata-store', 'data'),
         Output('processing-log', 'children'),
         Output('data-summary', 'children')],
        Input('upload-data-file', 'contents'),
        [State('upload-data-file', 'filename'),
         State('upload-input-type', 'value')]
    )
    def handle_file_upload(contents, filename, input_type):
        """Handle file upload and initial data loading."""
        if contents is None:
            raise PreventUpdate

        log = f"[{datetime.now().strftime('%H:%M:%S')}] Loading file: {filename}\n"

        try:
            # Decode the file
            content_type, content_string = contents.split(',')
            decoded = base64.b64decode(content_string)

            # Read the file based on type
            if input_type == 'h5ad':
                log += "Reading .h5ad file...\n"
                adata = ad.read_h5ad(io.BytesIO(decoded))
            elif input_type == 'h5':
                log += "Reading 10X HDF5 file...\n"
                adata = sc.read_10x_h5(io.BytesIO(decoded))
            else:
                return {'display': 'none'}, dbc.Alert("Unsupported file type", color="danger"), None, log, ""

            log += f"✓ File loaded successfully\n"
            log += f"  - Cells: {adata.n_obs}\n"
            log += f"  - Genes: {adata.n_vars}\n\n"

            # Check existing analyses
            log += "Checking existing analyses:\n"
            existing_analyses = check_adata_analyses(adata)

            for analysis, status in existing_analyses.items():
                if status:
                    log += f"  ✓ {analysis}\n"

            if not any(existing_analyses.values()):
                log += "  ⚠ No preprocessing detected. Full processing required.\n"

            log += "\n⚠ Click 'Process Data' to run analyses.\n"

            # Save to data store
            metadata = save_adata(adata)

            # Create data summary
            summary = f"""Input type: {input_type.upper()}
Number of cells: {adata.n_obs:,}
Number of genes: {adata.n_vars:,}

Existing analyses:
"""
            for analysis, status in existing_analyses.items():
                summary += f"  - {analysis}: {'Yes' if status else 'No'}\n"

            # Show success message
            status = dbc.Alert([
                html.I(className="fas fa-check-circle me-2"),
                f"File uploaded successfully: {filename}"
            ], color="success", className="mb-0")

            return {'display': 'block'}, status, metadata, log, summary

        except Exception as e:
            log += f"\n✗ Error loading file: {str(e)}\n"
            log += f"{traceback.format_exc()}\n"

            status = dbc.Alert([
                html.I(className="fas fa-exclamation-circle me-2"),
                f"Error loading file: {str(e)}"
            ], color="danger")

            return {'display': 'none'}, status, None, log, ""

    # Main processing callback

    @app.callback(
        [Output('processing-log', 'children', allow_duplicate=True),
         Output('data-summary', 'children', allow_duplicate=True),
         Output('processed-flags', 'data', allow_duplicate=True),
         Output('adata-store', 'data', allow_duplicate=True)],
        Input('run-processing-button', 'n_clicks'),
        [State('adata-store', 'data'),
         State('process-nfeatures', 'value'),
         State('process-dims', 'value'),
         State('process-resolution', 'value'),
         State('process-n-neighbors', 'value'),
         State('process-run-qc', 'value'),
         State('process-recluster', 'value'),
         State('qc-max-mt', 'value'),
         State('qc-min-ncount', 'value'),
         State('qc-max-ncount', 'value'),
         State('processed-flags', 'data'),
         State('processing-log', 'children')],
        prevent_initial_call=True
    )
    def process_data(n_clicks, adata_metadata, nfeatures, dims, resolution,
                     n_neighbors, run_qc, recluster, max_mt, min_ncount,
                     max_ncount, flags, existing_log):
        """
        Main data processing pipeline.
        Matches the R version's processing workflow.
        """
        if n_clicks is None or adata_metadata is None:
            raise PreventUpdate

        # Initialize log
        log = existing_log if existing_log else ""
        log += f"\n{'='*60}\n"
        log += f"[{datetime.now().strftime('%H:%M:%S')}] PROCESSING STARTED\n"
        log += f"{'='*60}\n\n"

        try:
            # Load data
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise ValueError("Could not load data from store")

            log += f"Parameters:\n"
            log += f"  - Variable features: {nfeatures}\n"
            log += f"  - PCs: {dims}\n"
            log += f"  - Resolution: {resolution}\n"
            log += f"  - Neighbors: {n_neighbors}\n"
            log += f"  - QC enabled: {'Yes' if run_qc and 'qc' in run_qc else 'No'}\n"
            log += f"  - Force recluster: {'Yes' if recluster and 'recluster' in recluster else 'No'}\n\n"

            # Check existing analyses
            existing = check_adata_analyses(adata)

            # QC Filtering
            if run_qc and 'qc' in run_qc:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Running QC filtering...\n"

                # Calculate QC metrics if not present
                if 'pct_counts_mt' not in adata.obs.columns:
                    calculate_qc_metrics(adata)

                # Filter cells
                adata, n_before, n_after = filter_cells_qc(
                    adata,
                    min_genes=200,
                    max_genes=2500,
                    min_counts=min_ncount,
                    max_counts=max_ncount,
                    max_mt_percent=max_mt
                )

                log += f"  ✓ QC complete. Cells: {n_before} → {n_after} "
                log += f"(removed {n_before - n_after})\n"
                log += f"    - Max MT%: {max_mt}\n"
                log += f"    - Min counts: {min_ncount}\n"
                log += f"    - Max counts: {max_ncount}\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping QC (not requested)\n\n"

            # Normalization
            if not existing['normalized']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Normalizing data...\n"
                normalize_data(adata)
                log += f"  ✓ Normalization complete\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping normalization (already done)\n\n"

            # Find variable genes
            if not existing['highly_variable_genes']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Finding variable features...\n"
                find_variable_genes(adata, n_top_genes=nfeatures)
                log += f"  ✓ {nfeatures} variable features identified\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping variable features (already done)\n\n"

            # Scale data
            if not existing['scaled']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Scaling data...\n"
                scale_data(adata)
                log += f"  ✓ Scaling complete\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping scaling (already done)\n\n"

            # PCA
            if not existing['pca']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Running PCA...\n"
                run_pca(adata, n_comps=50)
                log += f"  ✓ PCA complete\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping PCA (already done)\n\n"

            # Neighbors
            if not existing['neighbors']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Computing neighborhood graph...\n"
                compute_neighbors(adata, n_neighbors=n_neighbors, n_pcs=dims)
                log += f"  ✓ Neighbors computed\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping neighbors (already done)\n\n"

            # UMAP
            if not existing['umap']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Running UMAP...\n"
                run_umap(adata)
                log += f"  ✓ UMAP complete\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping UMAP (already done)\n\n"

            # t-SNE
            if not existing['tsne']:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] Running t-SNE...\n"
                run_tsne(adata, n_pcs=dims)
                log += f"  ✓ t-SNE complete\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping t-SNE (already done)\n\n"

            # Clustering
            should_cluster = (not existing['clusters']) or (
                recluster and 'recluster' in recluster)

            if should_cluster:
                action = "Re-clustering" if (
                    recluster and 'recluster' in recluster) else "Clustering"
                log += f"[{datetime.now().strftime('%H:%M:%S')}] {action} (resolution={resolution})...\n"
                find_clusters(adata, resolution=resolution)
                n_clusters = len(adata.obs['leiden'].unique())
                log += f"  ✓ {action} complete ({n_clusters} clusters)\n\n"
            else:
                log += f"[{datetime.now().strftime('%H:%M:%S')}] → Skipping clustering (already done)\n\n"

            # Cell cycle scoring
            log += f"[{datetime.now().strftime('%H:%M:%S')}] Cell cycle scoring...\n"
            score_cell_cycle(
                adata,
                s_genes=CELL_CYCLE_GENES['s_genes'],
                g2m_genes=CELL_CYCLE_GENES['g2m_genes']
            )
            log += f"  ✓ Cell cycle scoring complete\n\n"

            # Save updated data
            log += f"[{datetime.now().strftime('%H:%M:%S')}] Saving processed data...\n"
            updated_metadata = update_adata(session_id, adata)

            # Update flags
            if flags is None:
                flags = {}
            flags['uploaded'] = True
            flags['qc_done'] = run_qc and 'qc' in run_qc
            flags['processed'] = True

            # Final log
            log += f"\n{'='*60}\n"
            log += f"✓ ALL PROCESSING COMPLETE!\n"
            log += f"{'='*60}\n"

            # Create updated summary
            summary = f"""Number of cells: {adata.n_obs:,}
Number of genes: {adata.n_vars:,}
Number of clusters: {len(adata.obs['leiden'].unique())}

Reductions available: {', '.join(adata.obsm.keys())}
"""

            # Only show QC metrics if they exist
            if 'pct_counts_mt' in adata.obs.columns:
                summary += f"""
QC Metrics:
  - Mean MT%: {adata.obs['pct_counts_mt'].mean():.2f}
  - Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.0f}
  - Mean counts per cell: {adata.obs['total_counts'].mean():.0f}
"""

            # Only show cell cycle if it exists
            if 'phase' in adata.obs.columns:
                summary += f"""
Cell Cycle:
  - G1: {(adata.obs['phase'] == 'G1').sum()}
  - S: {(adata.obs['phase'] == 'S').sum()}
  - G2M: {(adata.obs['phase'] == 'G2M').sum()}
"""

            return log, summary, flags, updated_metadata

        except Exception as e:
            log += f"\n✗ ERROR during processing:\n"
            log += f"{str(e)}\n\n"
            log += f"Traceback:\n{traceback.format_exc()}\n"

            return log, "Processing failed - see log for details", flags, adata_metadata

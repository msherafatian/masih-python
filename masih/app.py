"""
Main application entry point for MASIH.

This module initializes the Dash application, sets up the layout,
and registers all callbacks.
"""

import dash
import dash_bootstrap_components as dbc
from dash import html, dcc

from masih.config.settings import AppConfig
from masih.layout.main_layout import create_layout, create_data_stores, create_tabs
from masih.callbacks.callback_manager import register_all_callbacks


# Initialize Dash app with Bootstrap theme
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    suppress_callback_exceptions=True,  # CRITICAL for dynamic layouts
    title=AppConfig.APP_NAME,
    update_title="Loading...",
    meta_tags=[
        {
            "name": "viewport",
            "content": "width=device-width, initial-scale=1.0"
        },
        {
            "name": "description",
            "content": AppConfig.APP_FULL_NAME
        }
    ]
)

# Server instance for deployment
server = app.server

# CRITICAL: Validation layout contains all IDs that callbacks will reference
# Even if components aren't currently visible, Dash needs to know they exist
app.validation_layout = html.Div([
    # Stores (always exist at top level)
    create_data_stores(),

    # Tabs (always exist)
    create_tabs(),

    # Tab content container
    html.Div(id='tab-content'),
    # CancerSEA module IDs
    html.Div([
        dcc.Interval(id='cancersea-reload-interval'),
        dcc.Dropdown(id='cancersea-pathway-select'),
        dbc.Button(id='calculate-cancersea-button'),
        dbc.Button(id='calculate-batch-cancersea-button'),
        dbc.Checklist(id='cancersea-batch-select'),
        html.Div(id='cancersea-calculation-status'),
        html.Div(id='cancersea-batch-area'),
        dcc.Graph(id='cancersea-feature-plot'),
        dcc.Graph(id='cancersea-violin-plot'),
        dcc.Graph(id='cancersea-heatmap'),
        html.Div(id='cancersea-stats-table'),
    ]),
    # Upload module IDs
    html.Div([
        dcc.Upload(id='upload-data-file'),
        html.Div(id='upload-status'),
        html.Div(id='processing-options-area'),
        html.Div(id='qc-parameters-area'),
        html.Div(id='upload-file-area'),
        html.Div(id='processing-log'),
        html.Div(id='data-summary'),
        dbc.Input(id='process-nfeatures'),
        dbc.Input(id='process-dims'),
        dbc.Input(id='process-resolution'),
        dbc.Input(id='process-n-neighbors'),
        dbc.Checklist(id='process-recluster'),
        dbc.Checklist(id='process-run-qc'),
        dbc.Input(id='qc-max-mt'),
        dbc.Input(id='qc-min-ncount'),
        dbc.Input(id='qc-max-ncount'),
        dbc.Button(id='run-processing-button'),
        dcc.Dropdown(id='upload-input-type'),
    ]),

    # Cluster Analysis module IDs
    html.Div([
        dcc.Dropdown(id='cluster-reduction-select'),
        dcc.Dropdown(id='cluster-color-by-select'),
        dcc.Slider(id='cluster-point-size'),
        dcc.Graph(id='cluster-plot'),
        dcc.Graph(id='cluster-tree-plot'),
        html.Div(id='cluster-stats-table'),
    ]),

    # Markers module IDs
    html.Div([
        dcc.Dropdown(id='marker-test-method'),
        dbc.RadioItems(id='marker-comparison-type'),
        dbc.Input(id='marker-min-logfc'),
        dbc.Input(id='marker-min-pct'),
        dbc.Checklist(id='marker-only-pos'),
        dbc.Input(id='marker-top-n'),
        dbc.Input(id='marker-max-display'),
        dcc.Dropdown(id='marker-cluster1'),
        dcc.Dropdown(id='marker-cluster2'),
        html.Div(id='marker-pairwise-selectors'),
        dbc.Button(id='calculate-markers-button'),
        html.Div(id='marker-calculation-status'),
        html.Div(id='marker-quick-actions'),
        html.Div(id='marker-table-container'),
        dcc.Graph(id='marker-heatmap'),
        dcc.Graph(id='marker-dotplot'),
        dcc.Dropdown(id='visualize-marker-gene'),
        html.Div(id='marker-heatmap-info'),
        dbc.RadioItems(id='marker-plot-type'),
        dbc.Button(id='update-marker-plot-button'),
        dcc.Graph(id='selected-marker-plot'),
        dbc.Button(id='download-all-markers-button'),
        dbc.Button(id='download-top-markers-button'),
        dbc.Button(id='clear-markers-button'),
        dcc.Download(id='download-all-markers'),
        dcc.Download(id='download-top-markers'),
    ]),
    # Cell Cycle module IDs
    html.Div([
        dcc.Interval(id='cellcycle-reload-interval'),
        dcc.Graph(id='cellcycle-phase-plot'),
        dcc.Graph(id='cellcycle-score-plot'),
        dcc.Graph(id='phase-by-cluster-plot'),
        dbc.RadioItems(id='cellcycle-score-select'),
        html.Div(id='cellcycle-stats-table'),
        html.Div(id='cellcycle-status-alert'),
    ]),
    # Placeholder IDs for future modules
    html.Div([
        html.Div(id='pathways-placeholder'),
        html.Div(id='trajectory-placeholder'),
        html.Div(id='compare-placeholder'),
        html.Div(id='explorer-placeholder'),
        html.Div(id='export-placeholder'),
    ]),
    # Trajectory module IDs
    html.Div([
        dcc.Interval(id='trajectory-reload-interval'),
        dcc.Dropdown(id='trajectory-annotation-col'),
        dcc.Dropdown(id='trajectory-celltypes'),
        dcc.Dropdown(id='trajectory-start-cluster'),
        dcc.Dropdown(id='trajectory-end-cluster'),
        dbc.Input(id='trajectory-n-neighbors'),
        dbc.Button(id='run-trajectory-button'),
        html.Div(id='trajectory-filter-preview'),
        html.Div(id='trajectory-status'),
        dbc.Tabs(id='trajectory-viz-tabs'),
        html.Div(id='trajectory-viz-content'),
        html.Div(id='trajectory-stats-display'),
    ]),
    # Compare module IDs
    html.Div([
        dcc.Interval(id='compare-reload-interval'),
        dcc.Dropdown(id='compare-pathways-select'),
        dbc.Button(id='run-comparison-button'),
        html.Div(id='compare-status'),
        dcc.Graph(id='compare-correlation-matrix'),
        dcc.Graph(id='compare-pathway-by-cluster'),
        dcc.Graph(id='compare-heatmap'),
        html.Div(id='compare-stats-table'),
    ]),
    # Export module IDs
    html.Div([
        dcc.Interval(id='export-reload-interval'),
        dcc.Dropdown(id='export-plot-type'),
        dcc.Dropdown(id='export-format'),
        dcc.Dropdown(id='export-dpi'),
        dbc.Input(id='export-width'),
        dbc.Input(id='export-height'),
        dbc.Checklist(id='export-white-bg'),
        dbc.Input(id='export-custom-title'),
        html.Div(id='export-plot-options'),
        dbc.Button(id='download-plot-button'),
        dbc.Button(id='download-all-plots-button'),
        dbc.Button(id='download-data-button'),
        dbc.Button(id='download-anndata-button'),
        dcc.Download(id='download-plot'),
        dcc.Download(id='download-all-plots'),
        dcc.Download(id='download-data'),
        dcc.Download(id='download-anndata'),
        dbc.Checklist(id='export-data-options'),
        dbc.Tabs(id='documentation-tabs'),
        html.Div(id='documentation-content'),
        dbc.Button(id='copy-methods-button'),
    ]),
    # Loading overlay
    html.Div(id='loading-output')
])

# Set the actual layout
app.layout = create_layout()

# Register all callbacks
register_all_callbacks(app)


def run_app(debug=False, host="127.0.0.1", port=8050):
    """
    Run the MASIH application.

    Args:
        debug: Enable debug mode (default: False)
        host: Host address (default: "127.0.0.1")
        port: Port number (default: 8050)
    """
    print(f"""
    ╔══════════════════════════════════════════════════════════════╗
    ║                          MASIH                               ║
    ║  Modular Analysis Suite for Interactive Heterogeneity        ║
    ║                      Version {AppConfig.VERSION}             ║
    ╚══════════════════════════════════════════════════════════════╝
    
    Starting application...
    Running on: http://{host}:{port}
    
    Press CTRL+C to quit
    """)

    app.run(
        debug=debug,
        host=host,
        port=port,
        dev_tools_hot_reload=debug
    )


if __name__ == "__main__":
    run_app(debug=True)

"""
Main layout structure for MASIH application.

This module defines the overall UI structure including header, navigation tabs,
data stores for state management, and the main content area.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html

from masih.config.settings import AppConfig


def create_header():
    """Create the application header with title and links."""
    return dbc.Navbar(
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.Div([
                        html.H2(
                            AppConfig.APP_NAME,
                            className="mb-0 text-white",
                            style={"fontWeight": "bold"}
                        ),
                        html.P(
                            AppConfig.APP_FULL_NAME,
                            className="mb-0 text-white-50",
                            style={"fontSize": "0.9rem"}
                        )
                    ])
                ], width="auto"),

                dbc.Col([
                    html.Div([
                        dbc.Button(
                            [html.I(className="fab fa-github me-2"), "GitHub"],
                            color="light",
                            outline=True,
                            href=AppConfig.GITHUB_URL,
                            target="_blank",
                            className="me-2"
                        ),
                        dbc.Button(
                            [html.I(className="fab fa-r-project me-2"),
                             "R Version"],
                            color="light",
                            outline=True,
                            href=AppConfig.R_VERSION_URL,
                            target="_blank"
                        )
                    ], className="d-flex justify-content-end")
                ], width=True)
            ], className="w-100 align-items-center justify-content-between")
        ], fluid=True),
        color="primary",
        dark=True,
        className="mb-4"
    )


def create_data_stores():
    """
    Create dcc.Store components for state management.

    These stores maintain application state across callbacks.
    Data is stored as JSON and persists during the session.
    """
    return html.Div([
        # Main data storage - stores session ID and metadata
        dcc.Store(
            id='adata-store',
            storage_type='session',
            data=None
        ),

        # Processing flags - tracks which analyses have been completed
        dcc.Store(
            id='processed-flags',
            storage_type='session',
            data={
                'uploaded': False,
                'qc_done': False,
                'processed': False,
                'markers_calculated': False,
                'pathways_calculated': [],
                'trajectory_calculated': False
            }
        ),

        # Marker genes results
        dcc.Store(
            id='markers-store',
            storage_type='session',
            data=None
        ),

        # Pathway scores
        dcc.Store(
            id='pathway-scores-store',
            storage_type='session',
            data=None
        ),

        # Trajectory results
        dcc.Store(
            id='trajectory-store',
            storage_type='session',
            data=None
        ),

        # UI state preferences
        dcc.Store(
            id='ui-state',
            storage_type='session',
            data={
                'active_tab': 'upload',
                'color_by': 'cluster',
                'reduction': 'umap'
            }
        )
    ])


def create_tabs():
    """Create the main navigation tabs."""
    tabs = dbc.Tabs(
        id="main-tabs",
        active_tab="upload",
        children=[
            dbc.Tab(
                label="Upload & QC",
                tab_id="upload",
                label_style={"cursor": "pointer"}
            ),
            dbc.Tab(
                label="Cluster Analysis",
                tab_id="cluster",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="cluster-tab"
            ),
            dbc.Tab(
                label="Marker Genes",
                tab_id="markers",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="markers-tab"
            ),
            dbc.Tab(
                label="CancerSEA Pathways",
                tab_id="pathways",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="pathways-tab"
            ),
            dbc.Tab(
                label="Trajectory",
                tab_id="trajectory",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="trajectory-tab"
            ),
            dbc.Tab(
                label="Compare",
                tab_id="compare",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="compare-tab"
            ),
            dbc.Tab(
                label="Explorer",
                tab_id="explorer",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="explorer-tab"
            ),
            dbc.Tab(
                label="Export",
                tab_id="export",
                label_style={"cursor": "pointer"},
                disabled=True,
                id="export-tab"
            )
        ],
        className="mb-3"
    )

    return tabs


def create_content_area():
    """Create the main content area that will be dynamically updated."""
    return html.Div(
        id="tab-content",
        children=[
            dbc.Container([
                dbc.Alert(
                    [
                        html.H4("Welcome to MASIH!",
                                className="alert-heading"),
                        html.P(
                            "Please upload your data to begin analysis.",
                            className="mb-0"
                        )
                    ],
                    color="info"
                )
            ])
        ]
    )


def create_footer():
    """Create the application footer."""
    return dbc.Container([
        html.Hr(),
        html.Div([
            html.P([
                f"MASIH v{AppConfig.VERSION} | ",
                html.A(
                    "Documentation",
                    href=f"{AppConfig.GITHUB_URL}#readme",
                    target="_blank"
                ),
                " | ",
                html.A(
                    "Report Issue",
                    href=f"{AppConfig.GITHUB_URL}/issues",
                    target="_blank"
                ),
                " | Built with ",
                html.A("Dash", href="https://dash.plotly.com", target="_blank"),
                " and ",
                html.A("Scanpy", href="https://scanpy.readthedocs.io",
                       target="_blank")
            ], className="text-center text-muted small mb-2"),
            html.P([
                "If you use MASIH in your research, please cite our paper."
            ], className="text-center text-muted small")
        ])
    ], fluid=True, className="mt-5 mb-3")


def create_layout():
    """
    Create the complete application layout.

    Returns:
        Dash layout component tree
    """
    return html.Div([
        # Data stores (hidden state management)
        create_data_stores(),

        # Header
        create_header(),

        # Main content
        dbc.Container([
            # Navigation tabs
            create_tabs(),

            # Dynamic content area
            create_content_area()
        ], fluid=True),

        # Footer
        create_footer(),

        # Loading overlay (for future use)
        dcc.Loading(
            id="loading-overlay",
            type="default",
            children=html.Div(id="loading-output")
        )
    ])

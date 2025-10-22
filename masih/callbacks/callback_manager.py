"""
Callback manager for MASIH application.

This module centralizes callback registration and manages the interaction
between different modules. It registers callbacks for tab switching and
coordinates module-specific callbacks.
"""

from dash import Input, Output, State, callback_context
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
from dash import html


def register_tab_switching(app):
    """
    Register callbacks for tab switching and enabling/disabling tabs.

    Tabs are enabled based on processing flags to ensure proper workflow.
    """

    @app.callback(
        Output("tab-content", "children"),
        Input("main-tabs", "active_tab"),
        State("processed-flags", "data")
    )
    def render_tab_content(active_tab, flags):
        """
        Render content for the active tab.

        Args:
            active_tab: ID of the currently active tab
            flags: Processing flags dictionary

        Returns:
            Layout for the active tab
        """
        if active_tab is None:
            raise PreventUpdate

        # Import module layouts dynamically to avoid circular imports
        if active_tab == "upload":
            from masih.modules.upload import create_upload_layout
            return create_upload_layout()

        elif active_tab == "cluster":
            from masih.modules.cluster_analysis import create_cluster_analysis_layout
            return create_cluster_analysis_layout()

        elif active_tab == "markers":
            return dbc.Container([
                dbc.Alert(
                    [
                        html.H4("Marker Genes Module",
                                className="alert-heading"),
                        html.P("Marker detection module will be implemented here."),
                    ],
                    color="info"
                )
            ])

        elif active_tab == "pathways":
            return dbc.Container([
                dbc.Alert(
                    [
                        html.H4("CancerSEA Pathways Module",
                                className="alert-heading"),
                        html.P("Pathway scoring module will be implemented here."),
                    ],
                    color="info"
                )
            ])

        elif active_tab == "trajectory":
            return dbc.Container([
                dbc.Alert(
                    [
                        html.H4("Trajectory Analysis Module",
                                className="alert-heading"),
                        html.P(
                            "Trajectory inference module will be implemented here."),
                    ],
                    color="info"
                )
            ])

        elif active_tab == "compare":
            return dbc.Container([
                dbc.Alert(
                    [
                        html.H4("Comparative Analysis Module",
                                className="alert-heading"),
                        html.P(
                            "Pathway comparison module will be implemented here."),
                    ],
                    color="info"
                )
            ])

        elif active_tab == "explorer":
            return dbc.Container([
                dbc.Alert(
                    [
                        html.H4("Interactive Explorer Module",
                                className="alert-heading"),
                        html.P("Cell explorer module will be implemented here."),
                    ],
                    color="info"
                )
            ])

        elif active_tab == "export":
            return dbc.Container([
                dbc.Alert(
                    [
                        html.H4("Export Module", className="alert-heading"),
                        html.P("Export functionality will be implemented here."),
                    ],
                    color="info"
                )
            ])

        return html.Div("Unknown tab")

    @app.callback(
        [
            Output("cluster-tab", "disabled"),
            Output("markers-tab", "disabled"),
            Output("pathways-tab", "disabled"),
            Output("trajectory-tab", "disabled"),
            Output("compare-tab", "disabled"),
            Output("explorer-tab", "disabled"),
            Output("export-tab", "disabled")
        ],
        Input("processed-flags", "data")
    )
    def update_tab_states(flags):
        """
        Enable or disable tabs based on processing state.

        Args:
            flags: Processing flags dictionary

        Returns:
            Tuple of boolean values for each tab's disabled state
        """
        if flags is None:
            # All tabs disabled if no data
            return True, True, True, True, True, True, True

        # Cluster tab: enabled after processing is done
        cluster_disabled = not flags.get('processed', False)

        # Markers tab: enabled after processing is done
        markers_disabled = not flags.get('processed', False)

        # Pathways tab: enabled after processing is done
        pathways_disabled = not flags.get('processed', False)

        # Trajectory tab: enabled after processing is done
        trajectory_disabled = not flags.get('processed', False)

        # Compare tab: enabled after processing is done
        compare_disabled = not flags.get('processed', False)

        # Explorer tab: enabled after processing is done
        explorer_disabled = not flags.get('processed', False)

        # Export tab: enabled after processing is done
        export_disabled = not flags.get('processed', False)

        return (
            cluster_disabled,
            markers_disabled,
            pathways_disabled,
            trajectory_disabled,
            compare_disabled,
            explorer_disabled,
            export_disabled
        )


def register_all_callbacks(app):
    """
    Register all application callbacks.

    This function is called from app.py to set up all callbacks.
    As we add modules, we'll import and register their callbacks here.

    Args:
        app: Dash application instance
    """

    # Register core tab switching callbacks
    register_tab_switching(app)

    # Register upload module callbacks
    from masih.modules.upload import register_upload_callbacks
    register_upload_callbacks(app)

    # Register cluster analysis module callbacks
    from masih.modules.cluster_analysis import register_cluster_analysis_callbacks
    register_cluster_analysis_callbacks(app)

    # Module-specific callbacks will be registered here as we build them
    # Example:
    # from masih.modules.cluster_analysis import register_cluster_callbacks
    # register_cluster_callbacks(app)

    print("✓ All callbacks registered successfully")

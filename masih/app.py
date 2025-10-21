"""
Main application entry point for MASIH.

This module initializes the Dash application, sets up the layout,
and registers all callbacks.
"""

import dash
import dash_bootstrap_components as dbc
from dash import html

from masih.config.settings import AppConfig
from masih.layout.main_layout import create_layout
from masih.callbacks.callback_manager import register_all_callbacks


# Initialize Dash app with Bootstrap theme
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    suppress_callback_exceptions=True,  # Required for dynamic callbacks
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

# Set the layout
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

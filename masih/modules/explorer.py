"""
Interactive Cell Explorer Module for MASIH application.

This module provides an interactive scatter plot for exploring individual cells
and their properties.
"""

import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State, callback_context
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import pandas as pd
import numpy as np

from masih.utils.data_store import get_adata
from masih.utils.plotting import get_color_palette


def create_explorer_layout():
    """
    Create the layout for the Interactive Cell Explorer module.

    Returns:
        Dash layout component
    """

    layout = dbc.Container([
        # Hidden interval for reload
        dcc.Interval(
            id='explorer-reload-interval',
            interval=500,
            max_intervals=1,
            disabled=False
        ),

        dbc.Row([
            dbc.Col([
                html.H3("Interactive Cell Explorer", className="mb-3"),
                html.P(
                    "Click on individual cells to explore their properties and gene expression.",
                    className="text-muted"
                )
            ])
        ], className="mb-4"),

        # Main Explorer Card
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Interactive Cell Explorer", className="mb-0")),
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                dbc.Label("Color by:"),
                                dcc.Dropdown(
                                    id='explorer-color-by',
                                    options=[],
                                    value='leiden',
                                    clearable=False
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Reduction:"),
                                dcc.Dropdown(
                                    id='explorer-reduction',
                                    options=[
                                        {'label': 'UMAP', 'value': 'umap'},
                                        {'label': 'tSNE', 'value': 'tsne'},
                                        {'label': 'PCA', 'value': 'pca'}
                                    ],
                                    value='umap',
                                    clearable=False
                                )
                            ], md=4),

                            dbc.Col([
                                dbc.Label("Gene expression:"),
                                dcc.Dropdown(
                                    id='explorer-gene',
                                    options=[
                                        {'label': 'None', 'value': 'none'}],
                                    value='none',
                                    placeholder="Type to search genes..."
                                )
                            ], md=4)
                        ], className="mb-3"),

                        dcc.Loading(
                            children=dcc.Graph(
                                id='explorer-plot',
                                config={
                                    'displayModeBar': True,
                                    'displaylogo': False,
                                    'modeBarButtonsToRemove': ['lasso2d', 'select2d']
                                },
                                style={'height': '600px'}
                            )
                        ),

                        html.Hr(),

                        html.H5("Selected Cell Information:"),
                        html.Div(id='explorer-cell-info', children=[
                            html.P("Click on a cell to see its information",
                                   className="text-muted",
                                   style={'padding': '20px'})
                        ])
                    ])
                ])
            ], width=12)
        ])

    ], fluid=True)

    return layout


def register_explorer_callbacks(app):
    """
    Register all callbacks for the explorer module.

    Args:
        app: Dash application instance
    """

    import traceback

    # Callback to update dropdown options

    @app.callback(
        [Output('explorer-color-by', 'options'),
         Output('explorer-color-by', 'value'),
         Output('explorer-gene', 'options')],
        Input('adata-store', 'data')
    )
    def update_explorer_options(adata_metadata):
        """Update available options for coloring and gene selection."""
        if adata_metadata is None:
            raise PreventUpdate

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Color by options - all obs columns
            color_options = []

            # Add common ones first
            priority_cols = ['leiden', 'louvain', 'phase', 'cell_type']
            for col in priority_cols:
                if col in adata.obs.columns:
                    color_options.append(
                        {'label': col.replace('_', ' ').title(), 'value': col})

            # Add remaining obs columns
            for col in adata.obs.columns:
                if col not in priority_cols:
                    color_options.append({'label': col, 'value': col})

            # Default color
            if 'leiden' in adata.obs.columns:
                default_color = 'leiden'
            elif 'louvain' in adata.obs.columns:
                default_color = 'louvain'
            else:
                default_color = adata.obs.columns[0]

            # Gene options - first 100 genes for performance
            gene_options = [{'label': 'None', 'value': 'none'}]
            gene_options.extend([
                {'label': gene, 'value': gene}
                for gene in adata.var_names[:100]
            ])

            return color_options, default_color, gene_options

        except Exception as e:
            print(f"Error updating explorer options: {e}")
            raise PreventUpdate

# Main callback to create explorer plot
    @app.callback(
        Output('explorer-plot', 'figure'),
        [Input('explorer-color-by', 'value'),
         Input('explorer-reduction', 'value'),
         Input('explorer-gene', 'value'),
         Input('explorer-reload-interval', 'n_intervals')],
        State('adata-store', 'data'),
        prevent_initial_call=False
    )
    def update_explorer_plot(color_by, reduction, gene, n, adata_metadata):
        """Update the explorer plot."""
        print(f"DEBUG: Explorer plot callback fired, n={n}")

        if adata_metadata is None:
            fig = go.Figure()
            fig.add_annotation(
                text="No data loaded",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False
            )
            return fig

        try:
            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get coordinates
            coord_key = f'X_{reduction}'
            if coord_key not in adata.obsm:
                fig = go.Figure()
                fig.add_annotation(
                    text=f"{reduction.upper()} not available. Please run {reduction.upper()} in processing.",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=12, color="red")
                )
                return fig

            coords = adata.obsm[coord_key][:, :2]

            # Get cell IDs
            cell_ids = adata.obs_names.tolist()

            # Create hover text
            hover_texts = []

            # Determine what to color by
            if gene and gene != 'none' and gene in adata.var_names:
                # Color by gene expression
                gene_idx = adata.var_names.get_loc(gene)
                if hasattr(adata.X, 'toarray'):
                    expression = adata.X[:, gene_idx].toarray().flatten()
                else:
                    expression = adata.X[:, gene_idx].flatten()

                # Create hover text
                for i, cell_id in enumerate(cell_ids):
                    text = f"Cell: {cell_id}<br>"
                    if color_by and color_by in adata.obs.columns:
                        text += f"{color_by}: {adata.obs.iloc[i][color_by]}<br>"
                    text += f"{gene}: {expression[i]:.2f}"
                    hover_texts.append(text)

                # Create scatter plot with continuous color scale
                fig = go.Figure(data=go.Scattergl(
                    x=coords[:, 0],
                    y=coords[:, 1],
                    mode='markers',
                    marker=dict(
                        size=3,
                        color=expression,
                        colorscale='Viridis',
                        showscale=True,
                        colorbar=dict(title="Expression"),
                        line=dict(width=0)
                    ),
                    text=hover_texts,
                    hovertemplate='%{text}<extra></extra>',
                    customdata=cell_ids
                ))

                title = f"{gene} Expression"

            else:
                # Color by metadata
                if color_by and color_by in adata.obs.columns:
                    color_data = adata.obs[color_by]

                    # Check if numeric (continuous) or categorical
                    is_numeric = pd.api.types.is_numeric_dtype(color_data)

                    print(
                        f"DEBUG: Coloring by {color_by}, is_numeric: {is_numeric}")

                    if is_numeric:
                        # CONTINUOUS DATA (like differentiation_score)
                        print(f"DEBUG: Treating {color_by} as continuous")

                        color_values = color_data.values

                        # Create hover text
                        for i, cell_id in enumerate(cell_ids):
                            text = f"Cell: {cell_id}<br>{color_by}: {color_values[i]:.3f}"
                            hover_texts.append(text)

                        # Create scatter plot with continuous color scale
                        fig = go.Figure(data=go.Scattergl(
                            x=coords[:, 0],
                            y=coords[:, 1],
                            mode='markers',
                            marker=dict(
                                size=3,
                                color=color_values,  # Use numpy array directly
                                colorscale='Viridis',
                                showscale=True,
                                colorbar=dict(title=color_by),
                                line=dict(width=0)
                            ),
                            text=hover_texts,
                            hovertemplate='%{text}<extra></extra>',
                            customdata=cell_ids
                        ))

                        title = f"Colored by {color_by}"

                    else:
                        # CATEGORICAL DATA
                        color_data = color_data.astype(str)
                        unique_vals = sorted(color_data.unique())

                        print(
                            f"DEBUG: Treating {color_by} as categorical with {len(unique_vals)} categories")

                        # Create hover text
                        for i, cell_id in enumerate(cell_ids):
                            text = f"Cell: {cell_id}<br>{color_by}: {color_data.iloc[i]}"
                            hover_texts.append(text)

                        # Check if too many categories
                        if len(unique_vals) > 50:
                            # Too many categories - use continuous color mapping
                            print(
                                f"DEBUG: Too many categories ({len(unique_vals)}), using continuous coloring")

                            # Map categories to numbers
                            category_map = {cat: i for i,
                                            cat in enumerate(unique_vals)}
                            color_numeric = [category_map[cat]
                                             for cat in color_data]

                            fig = go.Figure(data=go.Scattergl(
                                x=coords[:, 0],
                                y=coords[:, 1],
                                mode='markers',
                                marker=dict(
                                    size=3,
                                    color=color_numeric,  # Use list, not range
                                    colorscale='Viridis',
                                    showscale=False,
                                    line=dict(width=0)
                                ),
                                text=hover_texts,
                                hovertemplate='%{text}<extra></extra>',
                                customdata=cell_ids
                            ))

                        else:
                            # Reasonable number of categories - create separate traces
                            colors = get_color_palette(
                                len(unique_vals), "discrete")
                            color_map = dict(zip(unique_vals, colors))

                            # Create traces for each category
                            fig = go.Figure()

                            for val in unique_vals:
                                mask = color_data == val
                                indices = np.where(mask)[0]

                                fig.add_trace(go.Scattergl(
                                    x=coords[mask, 0],
                                    y=coords[mask, 1],
                                    mode='markers',
                                    name=str(val),
                                    marker=dict(
                                        size=3,
                                        color=color_map[val],
                                        line=dict(width=0)
                                    ),
                                    text=[hover_texts[i] for i in indices],
                                    hovertemplate='%{text}<extra></extra>',
                                    customdata=[cell_ids[i] for i in indices]
                                ))

                        title = f"Colored by {color_by}"

                else:
                    # Default: all same color
                    fig = go.Figure(data=go.Scattergl(
                        x=coords[:, 0],
                        y=coords[:, 1],
                        mode='markers',
                        marker=dict(size=3, color='blue', line=dict(width=0)),
                        text=[f"Cell: {cid}" for cid in cell_ids],
                        hovertemplate='%{text}<extra></extra>',
                        customdata=cell_ids
                    ))

                    title = f"{reduction.upper()} Plot"

            # Update layout
            fig.update_layout(
                title=title,
                xaxis_title=f"{reduction.upper()}_1",
                yaxis_title=f"{reduction.upper()}_2",
                height=600,
                template='plotly_white',
                hovermode='closest',
                clickmode='event+select'
            )

            return fig

        except Exception as e:
            print(f"Error creating explorer plot: {e}")
            traceback.print_exc()

            fig = go.Figure()
            fig.add_annotation(
                text=f"Error: {str(e)}",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=12, color="red")
            )
            return fig

    # Callback to show cell information on click

    @app.callback(
        Output('explorer-cell-info', 'children'),
        Input('explorer-plot', 'clickData'),
        [State('adata-store', 'data'),
         State('explorer-gene', 'value')],
        prevent_initial_call=True
    )
    def display_cell_info(click_data, adata_metadata, selected_gene):
        """Display information about clicked cell."""
        if click_data is None or adata_metadata is None:
            raise PreventUpdate

        try:
            # Get cell ID from click data
            cell_id = click_data['points'][0].get('customdata')

            if cell_id is None:
                return html.P("Could not identify cell", className="text-muted")

            session_id = adata_metadata['session_id']
            adata = get_adata(session_id)

            if adata is None:
                raise PreventUpdate

            # Get cell index
            try:
                cell_idx = adata.obs_names.get_loc(cell_id)
            except:
                return html.P(f"Cell {cell_id} not found", className="text-danger")

            # Get cell metadata
            cell_meta = adata.obs.iloc[cell_idx]

            # Create display
            info_content = []

            # Cell ID
            info_content.append(
                html.H6(f"Cell ID: {cell_id}", className="mb-3"))

            # Metadata
            info_content.append(html.H6("Cell Metadata:", className="mt-3"))
            meta_items = []
            for col in cell_meta.index[:10]:  # Show first 10 metadata columns
                meta_items.append(html.Li(f"{col}: {cell_meta[col]}"))
            info_content.append(html.Ul(meta_items))

            # Selected gene expression
            if selected_gene and selected_gene != 'none' and selected_gene in adata.var_names:
                gene_idx = adata.var_names.get_loc(selected_gene)
                if hasattr(adata.X, 'toarray'):
                    expr_val = adata.X[cell_idx, gene_idx].toarray()[0, 0]
                else:
                    expr_val = adata.X[cell_idx, gene_idx]

                info_content.append(
                    html.H6("Gene Expression:", className="mt-3"))
                info_content.append(html.P(f"{selected_gene}: {expr_val:.3f}"))

            # Top expressed genes
            info_content.append(
                html.H6("Top 10 Expressed Genes:", className="mt-3"))

            if hasattr(adata.X, 'toarray'):
                cell_expr = adata.X[cell_idx, :].toarray().flatten()
            else:
                cell_expr = adata.X[cell_idx, :].flatten()

            top_indices = np.argsort(cell_expr)[-10:][::-1]
            top_genes = [(adata.var_names[i], cell_expr[i])
                         for i in top_indices]

            gene_items = []
            for gene, expr in top_genes:
                gene_items.append(html.Li(f"{gene}: {expr:.2f}"))
            info_content.append(html.Ul(gene_items))

            return html.Div(info_content, style={'maxHeight': '400px', 'overflowY': 'auto'})

        except Exception as e:
            print(f"Error displaying cell info: {e}")
            traceback.print_exc()
            return html.P(f"Error: {str(e)}", className="text-danger")

    # Callback to disable reload interval

    @app.callback(
        Output('explorer-reload-interval', 'disabled'),
        Input('explorer-reload-interval', 'n_intervals')
    )
    def disable_explorer_interval(n):
        """Disable interval after first fire."""
        return True

import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import os

# Define shared results folder
RESULTS_FOLDER = "seurat_outputs"

# Initialize Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])

# Layout of the dashboard
app.layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Navbar(
            dbc.Container([
                html.A(
                    dbc.Row([
                        dbc.Col(html.Img(src="https://via.placeholder.com/50", height="40px")),
                        dbc.Col(dbc.NavbarBrand("Neocellomics Dashboard", className="ms-2")),
                    ], align="center", className="g-0"),
                    href="#",
                    style={"textDecoration": "none"}
                )
            ]),
            color="dark",
            dark=True,
            className="mb-4"
        ),
        dbc.Tabs([
            dbc.Tab(label="UMAP Results", tab_id="tab-umap"),
            dbc.Tab(label="Differential Genes", tab_id="tab-genes"),
            dbc.Tab(label="Cluster Summary", tab_id="tab-clusters")
        ], id="tabs", active_tab="tab-umap", className="mb-4"),
        html.Div(id="tab-content", className="p-4")
    ]
)

# Tab content callback
@app.callback(
    Output("tab-content", "children"),
    Input("tabs", "active_tab")
)
def render_tab_content(active_tab):
    if active_tab == "tab-umap":
        umap_plot_path = os.path.join(RESULTS_FOLDER, "umap_singleR_plot.png")
        if os.path.exists(umap_plot_path):
            return html.Div([
                html.H4("UMAP Visualization", className="text-primary"),
                html.Img(src=umap_plot_path, style={"width": "100%", "height": "auto"}),
            ])
        else:
            return html.P("UMAP plot not available. Please run the Streamlit analysis.", className="text-danger")

    elif active_tab == "tab-genes":
        genes_file = os.path.join(RESULTS_FOLDER, "differential_genes.csv")
        if os.path.exists(genes_file):
            df = pd.read_csv(genes_file)
            fig = px.scatter(
                df,
                x="avg_log2FC",
                y="p_val_adj",
                text="gene",
                title="Volcano Plot of Differential Genes",
                labels={"avg_log2FC": "Log2 Fold Change", "p_val_adj": "-Log10 P-value"}
            )
            return html.Div([
                html.H4("Differential Gene Expression", className="text-primary"),
                dcc.Graph(figure=fig)
            ])
        else:
            return html.P("Differential gene data not available. Please run the Streamlit analysis.", className="text-danger")

    elif active_tab == "tab-clusters":
        cell_type_summary_file = os.path.join(RESULTS_FOLDER, "cell_type_summary.csv")
        if os.path.exists(cell_type_summary_file):
            df = pd.read_csv(cell_type_summary_file)
            fig = px.pie(df, names="Cell Type", values="Count", title="Cell-Type Distribution")
            return html.Div([
                html.H4("Cluster Summary", className="text-primary"),
                dcc.Graph(figure=fig)
            ])
        else:
            return html.P("Cell type summary not available. Please run the Streamlit analysis.", className="text-danger")

    return html.P("Select a tab to view content.")

if __name__ == "__main__":
    app.run_server(debug=True, port=8050)

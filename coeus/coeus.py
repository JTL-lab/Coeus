import os
import time
import json

import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template
from dash.dependencies import Input, Output, State
from IPython.display import Markdown as md
from flask import jsonify

def get_gene_options():
    """
    Retrieves all gene names in alphabetical order for user selector dropdown component of dashboard.
    """
    all_filenames = [filename.split('.')[0] for filename in os.listdir('assets/clustermap/JSON')]
    files = [filename for filename in all_filenames if 'surrogates' not in filename]
    return sorted(list(set(files)), key=str.casefold)


def get_gene_names():
    """
    Retrieves all gene names in alphabetical order for fetching surrogates.
    """
    files = [filename.split('.')[0] for filename in os.listdir('assets/clustermap/JSON')]
    genes = [file for file in files if not 'surrogates' in file]
    return sorted(list(set(genes)), key=str.casefold)


def get_gene_surrogates(gene):
    """
    Retrieves data from surrogate textfile to display if 'Representative neighborhoods only' view mode is selected.
    """
    with open('assets/clustermap/JSON/surrogates/' + gene + '_surrogates.txt', 'r') as infile:
        data = infile.readlines()
    return data


def load_gene_json(gene):
    """
    Updates index.html for rendering gene order visualization using clustermap.js to load the JSON file corresponding
    to the gene the user chose from the dropdown.
    """
    with open('assets/clustermap/JSON/index.html', 'r') as infile:
        data = infile.readlines()

    # Modify which gene's JSON is being loaded
    data[89] = '\t\td3.json("' + gene + '.json"' + ')\n'

    # Edit index file
    with open('assets/clustermap/JSON/index.html', 'w') as outfile:
        outfile.writelines(data)


def filter_incomplete_neighborhoods(gene, neighborhood_size, surrogates=False):
    """
    Modifies data being rendered in gene order visualization: shows only complete neighborhoods
    """
    if surrogates:
        filename = 'assets/clustermap/JSON/' + gene + '_surrogates.json'
    else:
        filename = 'assets/clustermap/JSON/' + gene + '.json'

    with open('assets/clustermap/JSON/' + gene + '_surrogates.json', 'r') as infile:
        if len(infile.readlines()) != 0:
            infile.seek(0)
            json_data = json.load(infile)

    cluster_data = json_data["clusters"].copy()
    for cluster_id in cluster_data:
        genes = cluster_data[cluster_id]["loci"]["genes"].keys()
        # Remove clusters/genomes/neighborhoods below the size threshold that is specified
        if len(genes) < neighborhood_size:
            cluster_data.pop(cluster_id, None)

    return


# --------------------------------------------------- DASHBOARD --------------------------------------------------------
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.LUX])
server = app.server
app.config.suppress_callback_exceptions = True
app.title = 'Gene neighborhoods visualizer'

app.layout = html.Div(children=[
    html.Div(className='row',
             children=[

                 # Left hand pane for user controls
                 html.Div(className='three columns div--user-controls',
                          children=[

                              # Pane header: ARETE logo + dashboard title!
                              dbc.Row([
                                  html.Div(className='two columns',
                                           children=[
                                               html.Img(src=app.get_asset_url('arete_logo.png'),
                                                        style={'max-width': '400px', 'height': '90px', 'float': 'left',
                                                               'padding-left': '10px'})
                                           ]),
                                  html.Div(className='one column',
                                           children=[
                                                html.Img(src=app.get_asset_url('coeus.png'),
                                                         style={'max-width': '400px', 'height': '90px', 'float': 'left',
                                                               'padding-left': '170px', 'padding-top': '5px'})
                                               #html.H1('Coeus',
                                               #        style={'color': '#FFFFFF', 'text-align': 'right',
                                               #               'padding-left': '175px', 'padding-top': '25px'}),
                                           ])
                              ], style={'padding-top': '20px'}),

                              html.P("Coeus allows for easy comparison of a given "
                                     "gene's neighborhoods across multiple "
                                     "genomes, along with the neighborhoods' clusterings using three algorithms "
                                     "(UPGMA, MCL, and DBSCAN) to provide additional context on their similarities "
                                     "and differences.", style={'font-size': '10pt', 'padding-top': '20px',
                                                                'padding-left': '25px', 'padding-bottom': '10px'}),

                              # Pane content
                              dbc.Row([
                                  html.P('Select an AMR gene to analyze:', style={'color': '#FFFFFF', 'font-size': '13pt'}),
                                  dcc.Dropdown(id='gene-selector', options=get_gene_options(),
                                               value=get_gene_options()[0], className='gene-selector'),
                                  html.P('View gene neighborhoods according to:', style={'color': '#FFFFFF',
                                                                                         'font-size': '13pt',
                                                                                         'padding-top': '25px'}),
                                  dcc.RadioItems(id='clustermap-mode',
                                                 options=['All genomes', 'Unique neighborhoods only',
                                                          'Representative UPGMA cluster'],
                                                 labelStyle={'display': 'block'},
                                                 value='All genomes', style={'color': '#dc2284',
                                                                             'font-size': '11pt',
                                                                             'padding-right': '25px',
                                                                             'margin-right': '25px'}),
                                  html.Div([dcc.Slider(0, 1, 0.1, id='upgma-slider', value=0.5)], id='slider-div')
                                  #html.P('Color neighborhoods by: ', style={'color': '#FFFFFF',
                                  #                                        'font-size': '13pt',
                                  #                                        'padding-top': '25px'}),
                                  #dcc.RadioItems(id='clustermap-color-mode',
                                  #               options=['Gene Cluster', 'RGI Hit'],
                                  #               labelStyle={'display': 'block'},
                                  #               value='Gene Cluster', style={'color':'#FFFFFF',
                                  #                                            'font-size':'11pt'})
                              ], style={'padding-left': '20px', 'padding-right': '20px', 'padding-top': '20px'})

                          ]
                          ),

                 # Gene order visualizations pane
                 html.Div(className='six columns side-by-side bg-grey', children=[
                     dbc.Row([
                         dcc.Loading(
                             id='clustermap-loading',
                             children=[html.Div([html.Div(id='clustermap-loading-output')])],
                             type='circle',
                             style={'padding-top': '900px'}
                         ),
                         html.Iframe(src=app.get_asset_url('assets/clustermap/JSON/' + get_gene_names()[0] + '.html'),
                                     id='clustermap', title='clustermap'),

                         dcc.Markdown(get_gene_surrogates(get_gene_names()[0]), id='surrogates-md',
                                      style={'color': 'black', 'font-weight': 700, 'overflow-y': 'scroll'})
                     ])
                 ]),

                 # Neighborhoods clustering visualizations pane
                 html.Div(className='three columns side-by-side bg-grey', children=[
                     dcc.Loading(
                         id='clustering-loading',
                         children=[html.Div([html.Div(id='clustering-loading-output')])],
                         type='circle',
                         style={'padding-top': '900px'}
                     ),
                     dbc.Row([
                         html.Iframe(src=app.get_asset_url('clustering/UPGMA/acrB.html'),
                                     id='UPGMA', title='UPGMA', width='419px', height='316px')
                     ]),
                     dbc.Row([
                         html.Iframe(src=app.get_asset_url('clustering/MCL/acrB.html'),
                                     id='MCL', title='MCL', width='419px', height='316px')
                     ]),
                     dbc.Row([
                         html.Iframe(src=app.get_asset_url('clustering/DBSCAN/acrB.html'),
                                     id='DBSCAN', title='DBSCAN', width='419px', height='316px')
                     ])
                 ])
             ])
])


# ------------------------------- Update all graphs based on AMR gene dropdown selection -------------------------------
@app.callback(Output(component_id='clustermap', component_property='src'),
              [Input(component_id='gene-selector', component_property='value'),
               Input(component_id='clustermap-mode', component_property='value')])
def clustermap_callback(gene_value, mode_value):
    if mode_value == 'Unique neighborhoods only':
        return app.get_asset_url('clustermap/JSON/' + gene_value + '_surrogates.html')
    else:
        return app.get_asset_url('clustermap/JSON/' + gene_value + '.html')


@app.callback(Output(component_id='UPGMA', component_property='src'),
              [Input(component_id='gene-selector', component_property='value')])
def UPGMA_callback(input_value):
    return app.get_asset_url('clustering/UPGMA/' + input_value + '.html')


@app.callback(Output(component_id='MCL', component_property='src'),
              [Input(component_id='gene-selector', component_property='value')])
def MCL_callback(input_value):
    return app.get_asset_url('clustering/MCL/' + input_value + '.html')


@app.callback(Output(component_id='DBSCAN', component_property='src'),
              [Input(component_id='gene-selector', component_property='value')])
def DBSCAN_callback(input_value):
    return app.get_asset_url('clustering/DBSCAN/' + input_value + '.html')


# ------------------------------------ Show circle spinner while graphs are loading ------------------------------------
@app.callback(Output(component_id='clustermap-loading', component_property='children'),
              [Input(component_id='gene-selector', component_property='value')])
def clustermap_loading_callback(input_value):
    time.sleep(1)


# Show loading animation for UPGMA graph window when user is selecting gene from dropdown
@app.callback(Output(component_id='clustering-loading', component_property='children'),
              [Input(component_id='gene-selector', component_property='value')])
def clustering_loading_callback(input_value):
    time.sleep(1)


# ----------------- Hide slider for representative cluster height threshold only if radio button selected --------------
@app.callback(Output(component_id='slider-div', component_property='style'),
              [Input(component_id='clustermap-mode', component_property='value')])
def clustermap_loading_callback(input_value):
    if input_value == 'Representative UPGMA cluster':
        return {'display': 'inline'}
    else:
        return {'display': 'none'}


# ---------- Hide secondary window for showing genome surrogates unless 'Unique neighborhoods only' selected -----------
@app.callback(Output(component_id='surrogates-md', component_property='style'),
              [Input(component_id='clustermap-mode', component_property='value'),
               Input(component_id='gene-selector', component_property='value')])
def clustermap_loading_callback(input_mode, input_value):
    get_gene_surrogates(input_value)
    if input_mode == 'Unique neighborhoods only':
        return {'height': 200, 'width': 2400, 'display': 'inline'}
    else:
        return {'height': 200, 'width': 2400, 'display': 'none'}


# ------- Show surrogates list in markdown format below clustermap viz if 'Unique neighborhoods only' is selected ------
@app.callback(Output(component_id='clustermap', component_property='style'),
              [Input(component_id='clustermap-mode', component_property='value')])
def clustermap_loading_callback(input_value):
    if input_value == 'Unique neighborhoods only':
        return {'height': 750, 'width': 2400}
    else:
        return {'height': 950, 'width': 2400}


# -------------------- Load proper surrogates file into markdown section according to gene selected --------------------
@app.callback(Output(component_id='surrogates-md', component_property='children'),
              [Input(component_id='gene-selector', component_property='value')])
def surrogates_loading_callback(input_value):
    return get_gene_surrogates(input_value)


if __name__ == '__main__':
    app.run_server(debug=True)

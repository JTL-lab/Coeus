import os
import time

import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output


def get_AMR_gene_names():
    genes = [filename.split('.')[0] for filename in os.listdir('assets/clustermap/JSON')]
    return sorted(list(set(genes)), key=str.casefold)


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
                                               html.H1('Neighborhoods visualizer',
                                                       style={'color': '#FFFFFF', 'text-align': 'right',
                                                              'padding-left': '175px', 'padding-top': '25px'}),
                                           ])
                              ], style={'padding-top': '20px'}),

                              html.P("The neighborhood visualizer allows for easy comparison of a given "
                                     "antimicrobial resistance (AMR) gene's neighborhoods across multiple "
                                     "genomes, along with the neighborhoods' clusterings using three algorithms "
                                     "(UPGMA, MCL, and DBSCAN) to provide additional context on their similarities "
                                     "and differences.", style={'font-size': '10pt', 'padding-top': '20px',
                                                                'padding-left': '25px', 'padding-bottom': '10px'}),

                              # Pane content
                              dbc.Row([
                                  html.P('Select an AMR gene to analyze:', style={'color': '#FFFFFF', 'font-size': '13pt'}),
                                  dcc.Dropdown(id='gene-selector', options=get_AMR_gene_names(),
                                               value='acrB', className='gene-selector'),
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
                         html.Iframe(src=app.get_asset_url('clustermap/JSON/acrB.html'),
                                     id='clustermap', title='clustermap',
                                     height='950vh', width='2400vw'),
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
              [Input(component_id='gene-selector', component_property='value')])
def clustermap_callback(input_value):
    return app.get_asset_url('clustermap/JSON/' + input_value + '.html')


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


if __name__ == '__main__':
    app.run_server(debug=True)

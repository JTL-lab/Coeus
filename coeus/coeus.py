import os
import sys
import time
import json
import numpy as np
import pandas as pd

import dash
from dash import dcc, DiskcacheManager
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State, ClientsideFunction
from dash import dash_table
import functools

from scipy.cluster import hierarchy
from scipy import sparse
import networkx as nx
import plotly.express as px
import plotly.graph_objects as go
from plotly.figure_factory import create_dendrogram
import markov_clustering as mcl
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from skbio.stats.ordination import pcoa
import random

from uuid import uuid4
import diskcache

# ------------------------------------------------ LOCAL CACHING SETUP  ------------------------------------------------
# Initialize caching: assumes app is launched locally
launch_uid = uuid4()
cache = diskcache.Cache("./cache")
background_callback_manager = DiskcacheManager(
    cache, cache_by=[lambda: launch_uid]
)

# ------------------------------------------- GENERAL APPLICATION METHODS ----------------------------------------------
def get_colors(num):
    """
    Color mapping for visualizations.
    """
    colors = ['#6e40aa', '#b83cb0', '#c33dad', '#ff4f7c', '#f6478d', '#ff6956', '#f59f30', '#c4d93e',
              '#83f557', '#38f17a', '#22e599', '#19d3b5', '#29a0dd', '#5069d9', '#5f56c9', '#bbbbbb']

    if num <= len(colors):
        return colors
    else:
        for i in range(num - len(colors)):
            rand_color = lambda: random.randint(0, 255)
            colors.append('#%02X%02X%02X' % (rand_color(), rand_color(), rand_color()))
        return colors


def get_gene_options():
    """
    Retrieves all gene names in alphabetical order for user selector dropdown component of dashboard.
    """
    try:
        json_files = os.listdir('assets/clustermap/JSON')
    except FileNotFoundError:
        print("Error: JSON directory missing. Please follow documentation instructions for JSON directory placement.")
        sys.exit(1)
    all_filenames = [filename.split('.')[0] for filename in os.listdir('assets/clustermap/JSON')]
    files = [filename for filename in all_filenames if 'surrogates' not in filename and 'upgma' not in filename]
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


def gene_surrogates_to_df(gene):
    """
    Retrieves data from surrogate textfile as a Pandas dataframe for easy loading into a Dash Datatable.
    """
    try:
        df = pd.read_csv('assets/clustermap/JSON/surrogates/' + gene + '_surrogates.txt', delimiter=':',
                     names=['Surrogate Genome', 'Represented Genomes'])
        df['Represented Genomes'] = df['Represented Genomes'].str.strip(' []')

    except AttributeError:
        df = pd.read_csv('assets/clustermap/JSON/surrogates/' + gene + '_surrogates.txt', delimiter=':',
                         names=['Surrogate Genome'])

    return df


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

    return cluster_data


def load_similarity_matrix(gene):
    """
    Loads similarity matrix required for rendering MCL clustering.
    """
    df = pd.read_csv('assets/clustering/similarity_matrices/' + gene + '_similarity_matrix.csv', sep='\t')
    return df


def load_distance_matrix(gene):
    """
    Loads distance matrix required for rendering UPGMA and DBSCAN clustering.
    """
    df = pd.read_csv('assets/clustering/distance_matrices/' + gene + '_distance_matrix.csv', sep='\t')
    return df


def get_UPGMA_heights(gene):
    """
    Obtains all unique node heights from the UPGMA dendrogram.
    """
    distance_matrix_df = load_distance_matrix(gene)
    genome_names = list(distance_matrix_df.columns)
    linkage_matrix = hierarchy.linkage(distance_matrix_df.values, method='average')
    dendrogram = hierarchy.dendrogram(linkage_matrix, labels=genome_names)

    node_dcoords = dendrogram['dcoord']
    node_heights = []
    for genome_coords in dendrogram['dcoord']:
        node_heights.append(genome_coords[1])

    unique_node_heights = list(set(node_heights))

    return unique_node_heights


def get_UPGMA_filtering_params(gene):
    node_heights = get_UPGMA_heights(gene)
    sorted_node_heights = sorted(node_heights)

    # Get min and max values
    min_height = round(min(node_heights), 2)
    max_height = round(max(node_heights), 2)

    return dcc.RadioItems(sorted_node_heights, value=max_height, id='upgma-filter')


# ----------------------------------- FUNCTIONS FOR APPLYING CLUSTERING ALGORITHMS  ------------------------------------
@functools.lru_cache(maxsize=32)
def render_UPGMA(gene):
    """
    Generates a UPGMA graph by applying the figure factory dendrogram to the selected gene's distance matrix.
    """
    distance_matrix_df = load_distance_matrix(gene)
    genome_names = list(distance_matrix_df.columns)
    dendrogram_fig = plotly_dendrogram(distance_matrix_df.values, genome_names, gene)

    linkage_matrix = hierarchy.linkage(distance_matrix_df.values, method='average')
    dendrogram = hierarchy.dendrogram(linkage_matrix, labels=genome_names)
    return dendrogram_fig


@functools.lru_cache(maxsize=32)
def render_MCL(gene, inflation=2):
    """
    Generates an MCL network using markov clustering, networkx, and Plotly go.
    """
    similarity_matrix_df = load_similarity_matrix(gene)
    genome_names = list(similarity_matrix_df.columns)

    np_matrix = np.array(similarity_matrix_df.values)
    sparse_matrix = sparse.csr_matrix(np_matrix)

    result = mcl.run_mcl(sparse_matrix, inflation=inflation)
    clusters = mcl.get_clusters(result)
    fig = plotly_mcl_network(sparse_matrix, clusters, genome_names, gene)
    return fig


@functools.lru_cache(maxsize=32)
def render_DBSCAN(gene, min_samples=5, eps=0.5):
    """
    Applies DBSCAN clustering to a neighborhood similarity or symmetric distance matrix.
    """
    distance_matrix_df = load_distance_matrix(gene)
    genome_names = list(distance_matrix_df.columns)

    np_distance_matrix = np.array(distance_matrix_df.values)
    distance_matrix = StandardScaler().fit_transform(np_distance_matrix)

    dbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(distance_matrix)
    fig = plotly_pcoa(distance_matrix_df, genome_names, dbscan.labels_, gene)
    return fig


# ----------------------------------- FUNCTIONS FOR GENERATING CLUSTERING GRAPHS ---------------------------------------
def plotly_dendrogram(linkage_matrix, genome_names, AMR_gene):
    """
    Generates an interactive dendrogram visualization using Plotly figure factory.
    """
    title = "UPGMA dendrogram for {g}".format(g=AMR_gene)
    # The create_dendrogram function wraps the scipy.hierarchy.dendrogram function. It automatically calculates the
    # condensed distance matrix using pdist assuming Euclidean distances and obtains the UPGMA linkage matrix.
    fig = create_dendrogram(linkage_matrix,
                            linkagefun=lambda x: hierarchy.linkage(x, "average"),
                            labels=genome_names,
                            colorscale=get_colors(len(genome_names)))
    fig.update_layout(autosize=True, title=title, paper_bgcolor='white', template='plotly_white', width=419, height=316)
    return fig


def plotly_mcl_network(matrix, clusters, genome_names, AMR_gene):
    # make a networkx graph from the adjacency matrix
    graph = nx.Graph(matrix)
    pos = nx.spring_layout(graph)

    # Get node cluster assignments from MCL
    cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}

    # Determine edge lines between nodes: only add edges between nodes in the same cluster
    edge_x = []
    edge_y = []

    for cluster_group in clusters:
        cluster_nodes = list(cluster_group)
        for edge in graph.edges(cluster_nodes):
            if edge[0] in cluster_nodes and edge[1] in cluster_nodes:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.append(x0)
                edge_x.append(x1)
                edge_x.append(None)
                edge_y.append(y0)
                edge_y.append(y1)
                edge_y.append(None)

    # Assign a random color to every cluster
    hex_colors = []
    cluster_colors_dict = {}

    colors = get_colors(len(genome_names))
    for hex in colors:
        hex_colors.append(hex)

    for cluster in set(cluster_map.values()):
        cluster_colors_dict[str(cluster)] = hex_colors[cluster]

    # Get graph nodes
    node_x = []
    node_y = []
    for node in graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    # Add genome_ID as hover text to allow user to see which node comes from which genome neighborhood
    node_cluster = []
    node_text = {}

    node_to_genome_id_map = {}
    for node, genome_id in zip(cluster_map.keys(), genome_names):
        node_to_genome_id_map[node] = genome_id

    for i, cluster in enumerate(clusters):
        for node in cluster:
            node_cluster.append(i)
            node_text[node] = "{g}: Cluster {c}".format(g=node_to_genome_id_map[node], c=cluster_map[node])

    # Draw edges
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        showlegend=False,
        hoverinfo='none',
        mode='lines'
    )

    # Baseline figure with graph edges
    fig = go.Figure(data=edge_trace,
                    layout=go.Layout(
                        title='MCL network clusters for {}'.format(AMR_gene),
                        autosize=False,
                        width=419,
                        height=316,
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    # Draw graph nodes: overlay on top of edge traces at respective positions
    legend_clusters = []
    for node in graph.nodes():
        x_pos, y_pos = pos[node]
        if cluster_map[node] not in legend_clusters:
            fig.add_trace(
                go.Scatter(x=[x_pos], y=[y_pos], name=str(cluster_map[node]), mode='markers',
                           hoverinfo='text', text=[node_text[node]], line=dict(width=2, color='black'),
                           legendgroup=str(cluster_map[node]), legendrank=cluster_map[node],
                           marker=dict(color=cluster_colors_dict[str(cluster_map[node])], size=15)))
            legend_clusters.append(cluster_map[node])

        else:
            fig.add_trace(
                go.Scatter(x=[x_pos], y=[y_pos], mode='markers',
                           hoverinfo='text', text=[node_text[node]],
                           legendgroup=str(cluster_map[node]), showlegend=False, legendrank=cluster_map[node],
                           marker=dict(color=cluster_colors_dict[str(cluster_map[node])], size=15)))

    fig.update_layout(legend_title='MCL Cluster', font_family='"Open Sans", verdana, arial, sans-serif',
                      paper_bgcolor='white', template='plotly_white')

    return fig


def plotly_pcoa(distance_matrix_df, genome_ids, labels, AMR_gene):
    """
    Make Plotly dash interactive scatterplot of PCoA visualization of DBSCAN clusters
    """
    # Get PCoA dimension values from distance matrix: preserve their relative relationships
    pcoa_vals = pcoa(distance_matrix_df)

    # Create new dataframe indexed by genome name containing data on PC1, PC2, and cluster label
    df_data = {'GenomeID': genome_ids,
               'PC1': pcoa_vals.samples['PC1'].values,
               'PC2': pcoa_vals.samples['PC2'].values,
               'Cluster': labels}

    df = pd.DataFrame(data=df_data, index=genome_ids)
    df['Cluster'] = df['Cluster'].astype(str)

    # Make noise cluster black by default
    colors = get_colors(len(genome_ids))
    colors.insert(0, '#111111')

    fig = px.scatter(df, x='PC1', y='PC2',
                     color='Cluster',
                     color_discrete_sequence=colors,
                     hover_name='GenomeID',
                     title='PCoA DBSCAN clusters for {g}'.format(g=AMR_gene))
    fig.update_traces(marker_size=5, line=dict(width=2, color='black'))
    fig.update_layout(paper_bgcolor='white', font_family='"Open Sans", verdana, arial, sans-serif',
                      template='plotly_white', width=419, height=316)

    return fig


# --------------------------------------------------- DASHBOARD --------------------------------------------------------
app = dash.Dash(__name__, background_callback_manager=background_callback_manager,
                external_stylesheets=[dbc.themes.LUX],
                external_scripts=['https://cdn.plot.ly/plotly-2.3.0.min.js']
                )
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
                              html.Div([
                                  dbc.Row([
                                      html.Div(className='one column',
                                      children=[html.Img(src=app.get_asset_url('coeus-logo-dark.svg'),
                                                         style={'float': 'left'}
                                                         )
                                  ])
                              ]),
                              ], style={'padding-top': '20px', 'align-items': 'center', 'display': 'flex', 'justify-content': 'center'}),


                              html.P("Coeus allows for easy comparison of a given "
                                     "gene's neighborhoods across multiple "
                                     "genomes, along with the neighborhoods' clusterings using three algorithms "
                                     "(UPGMA, MCL, and DBSCAN) to provide additional context on their similarities "
                                     "and differences.", style={'font-size': '10pt', 'padding-top': '10px',
                                                                'padding-left': '25px', 'padding-bottom': '11px'}),

                              # Pane content
                              dbc.Row([
                                  html.P('Select a gene to analyze:', style={'color': '#FFFFFF',
                                                                             'font-size': '13pt'}),
                                  dcc.Dropdown(id='gene-selector', options=get_gene_options(),
                                               value=get_gene_options()[0], className='gene-selector'),
                                  html.P('View gene neighborhoods according to:', style={'color': '#FFFFFF',
                                                                                         'font-size': '13pt',
                                                                                         'padding-top': '25px'}),
                                  dcc.RadioItems(id='clustermap-mode',
                                                 options=['All genomes', 'Unique neighborhoods only',
                                                          'Representative UPGMA cluster'],
                                                 labelStyle={'display': 'block'},
                                                 value='Unique neighborhoods only', style={'color': '#dc2284',
                                                                                           'font-size': '11pt',
                                                                                           'padding-right': '25px',
                                                                                           'padding-bottom': '25px',
                                                                                           'margin-right': '25px'}),
                                  html.Button('Toggle clustering hyperparameters',
                                              id='hyperparameter-toggle-btn',
                                              n_clicks=0,
                                              style={'background-color': '#7571B4', 'color': '#393939'}),
                                  html.P(
                                      'Note: For larger datasets, hyperparameter changes may take some time to update.',
                                      style={'font-size': '11pt', 'padding-top': '25px'}),
                                  html.Div([
                                      html.P('Toggle hyperparameters:', style={'color': '#FFFFFF',
                                                                               'font-size': '13pt',
                                                                               'padding-top': '25px'}),
                                      html.P('MCL Inflation:', style={'color': '#dc2284', 'font-size': '13pt',
                                                                      'padding-top': '10px'}),
                                      dcc.Slider(2, 20, 1, value=2, id='inflation-slider'),
                                      html.P('DBSCAN Minimum Points:', style={'color': '#dc2284', 'font-size': '13pt',
                                                                              'padding-top': '10px'}),
                                      dcc.Slider(1, 20, 1, value=5, id='min-pts-slider'),
                                      html.P('DBSCAN Epsilon:', style={'color': '#dc2284', 'font-size': '13pt',
                                                                       'padding-top': '10px'}),
                                      dcc.Slider(0.1, 1, 0.1, value=0.5, id='epsilon-slider')
                                  ], id='clustering-hyperparameters-div'),
                              ], style={'padding-left': '20px', 'padding-right': '20px', 'padding-top': '20px'})
                          ],
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
                         dbc.Table.from_dataframe(df=gene_surrogates_to_df(get_gene_names()[0]),
                                   id='surrogates-table', bordered=True, hover=True,
                                   style={'color': 'black', 'font-size': 12, 'height': 200, 'width': 2400,
                                          'display': 'inline', 'overflow': 'auto'})
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
                         html.Iframe(src=app.get_asset_url('clustering/UPGMA/' + get_gene_names()[0] + '.html'),
                                     id='UPGMA-iframe', title='UPGMA', width='419px', height='316px')
                     ]),
                     dbc.Row([
                         html.Iframe(src=app.get_asset_url('clustering/MCL/' + get_gene_names()[0] + '.html'),
                                     id='MCL-iframe', title='MCL', width='419px', height='316px'),
                         dcc.Graph(figure={}, id='MCL-fig', style={'display': 'none'})
                     ]),
                     dbc.Row([
                         html.Iframe(src=app.get_asset_url('clustering/DBSCAN/' + get_gene_names()[0] + '.html'),
                                     id='DBSCAN-iframe', title='DBSCAN', width='419px', height='316px'),
                         dcc.Graph(figure={}, id='DBSCAN-fig', style={'display': 'none'})
                     ])
                 ])
             ])
])


# ----------- If user clicks toggle button, allow them to see controls for hyperparameter tuning -----------------------
@app.callback(Output(component_id='clustering-hyperparameters-div', component_property='style'),
              [Input(component_id='hyperparameter-toggle-btn', component_property='n_clicks')])
def clustering_hyperparameters_callback(num_clicks):
    if num_clicks % 2 == 0:
        return {'display': 'none'}
    else:
        return {'display': 'inline'}


# ------------------------------- Update all graphs based on AMR gene dropdown selection -------------------------------
@app.callback(Output(component_id='clustermap', component_property='src'),
              [Input(component_id='gene-selector', component_property='value'),
               Input(component_id='clustermap-mode', component_property='value')])
def clustermap_callback(gene_value, mode_value):
    if mode_value == 'All genomes':
        return app.get_asset_url('clustermap/JSON/' + gene_value + '.html')
    elif mode_value == 'Representative UPGMA cluster':
        return app.get_asset_url('clustermap/JSON/' + gene_value + '_upgma.html')
    else:
        return app.get_asset_url('clustermap/JSON/' + gene_value + '_surrogates.html')


@app.callback(Output(component_id='UPGMA-iframe', component_property='src'),
              [Input(component_id='gene-selector', component_property='value')])
def UPGMA_callback(gene):
    return app.get_asset_url('clustering/UPGMA/' + gene + '.html')


@app.callback(Output(component_id='MCL-fig', component_property='figure'),
              Output(component_id='MCL-fig', component_property='style'),
              Output(component_id='MCL-iframe', component_property='style'),
              Output(component_id='MCL-iframe', component_property='src'),
              [Input(component_id='gene-selector', component_property='value'),
               Input(component_id='hyperparameter-toggle-btn', component_property='n_clicks'),
               Input(component_id='inflation-slider', component_property='value')])
def MCL_callback(gene, num_clicks, inflation_value=2):
    if num_clicks % 2 == 0:
        return {}, {'display': 'none'}, {}, app.get_asset_url('clustering/MCL/' + gene + '.html')
    elif num_clicks % 2 != 0:
        return render_MCL(gene, inflation_value), {}, {'display': 'none'}, \
               app.get_asset_url('clustering/MCL/' + gene + '.html')


# ----------------------- Redo clustering with user selected hyperparameters if sliders are used -----------------------
@app.callback(Output(component_id='DBSCAN-fig', component_property='figure'),
              Output(component_id='DBSCAN-fig', component_property='style'),
              Output(component_id='DBSCAN-iframe', component_property='style'),
              Output(component_id='DBSCAN-iframe', component_property='src'),
              [Input(component_id='gene-selector', component_property='value'),
               Input(component_id='hyperparameter-toggle-btn', component_property='n_clicks'),
               Input(component_id='min-pts-slider', component_property='value'),
               Input(component_id='epsilon-slider', component_property='value')])
def DBSCAN_callback(gene, num_clicks, min_samples=5, epsilon=0.5):
    if num_clicks % 2 == 0:
        return {}, {'display': 'none'}, {}, app.get_asset_url('clustering/DBSCAN/' + gene + '.html')
    elif num_clicks % 2 != 0:
        return render_DBSCAN(gene, min_samples=min_samples, eps=epsilon), {}, {'display': 'none'}, app.get_asset_url(
            'clustering/DBSCAN/' + gene + '.html')


# ------------------------------------ Show circle spinner while graphs are loading ------------------------------------
@app.callback(Output(component_id='clustermap-loading', component_property='children'),
              [Input(component_id='gene-selector', component_property='value')])
def clustermap_loading_callback(input_value):
    time.sleep(1)


@app.callback(Output(component_id='clustering-loading', component_property='children'),
              [Input(component_id='MCL', component_property='value')])
def clustering_loading_callback(input_value):
    time.sleep(1)


# ---------- Hide secondary window for showing genome surrogates unless 'Unique neighborhoods only' selected -----------
@app.callback(Output(component_id='surrogates-table', component_property='children'),
              Output(component_id='surrogates-table', component_property='style'),
              [Input(component_id='clustermap-mode', component_property='value'),
               Input(component_id='gene-selector', component_property='value')])
def clustermap_loading_callback(input_mode, gene):
    if input_mode == 'Unique neighborhoods only' and gene is not None:
        df = gene_surrogates_to_df(gene)
        return dbc.Table.from_dataframe(df, bordered=True, hover=True), \
               {'color': 'black', 'font-size': 12, 'height': 200, 'width': 2400,
                'display': 'inline', 'overflow': 'auto'}
    else:
        return '', {'height': 200, 'width': 2400, 'display': 'none'}


# ------- Show surrogates list in markdown format below clustermap viz if 'Unique neighborhoods only' is selected ------
@app.callback(Output(component_id='clustermap', component_property='style'),
              [Input(component_id='clustermap-mode', component_property='value')])
def clustermap_loading_callback(input_value):
    if input_value == 'Unique neighborhoods only':
        return {'height': 750, 'width': 2400}
    else:
        return {'height': 950, 'width': 2400}


if __name__ == '__main__':
    app.run_server(debug=True)

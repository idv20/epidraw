import pandas as pd
from pathlib import Path, PurePath
import igraph
import plotly.express as px
import math
import random
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import colorcet as cc
matplotlib.use('Agg')

def get_graph(exp, out_dir):
    """Reads graph and creates igraph object.

    Args:
        exp (string): Experiment name.
        out_dir (string): Output directory.

    Returns:
        igraph: Igraph object containing embedding graph.
    """
    # Read graph if it exists
    my_graph = Path(out_dir)/f'{exp}_graph_leiden.pkl'

    if my_graph.is_file():

        graph = pd.read_pickle(my_graph)

        # Construct igraph object
        Gm = igraph.Graph.TupleList(
                                    graph.itertuples(index=False),
                                    directed = False,
                                    edge_attrs= ['weight']
                                    )

    else:
        print(
            'Could NOT load GRAPH.\
            Make sure GRAPH was generated!')

    return Gm

def cluster_leiden(Gm, df_peaks, res, fig_dir, exp):
    """Clusters the data using Leiden algorithm. Option to use user-defined resolution or
    calculate a suggested one based on stability or improvement of clustering.

    Args:
        Gm (igraph): Igraph object containing embedding graph.
        df_peaks (dataframe): Dataframe of peaks with embedding.
        res (string): 'auto' or user-defined value.
        fig_dir (path): Path to figures directory.
        exp (string): Name of experiment.

    Returns:
        dataframe: Dataframe of peaks with embedding and assigned cluster.
        float: User-defined resolution or suggested resolution.
    """
    # Provides a possible resolution parameter of the most 'stable' clustering in a range of res
    if res == 'auto':
        # List of tuples (res, max_clust_number)
        results = []
        # List of resolution parameters
        res_list = []
        # List of max_clust_number
        nr_clust_list = []
        # Previous number of cluster
        previous_n_clust = 0

        # Perform Leiden clustering
        for i in range(30):
            random.seed(1234)
            res_var = 0.00001 * i
            clusters = Gm.community_leiden(objective_function='CPM',
                                            weights='weight',
                                            resolution=res_var, # for py=3.8+
                                            #resolution_parameter=res_var,
                                            n_iterations=-1)
            clusters = pd.Series(clusters.membership, index=Gm.vs["name"]).sort_index().values
            # Current number of clusters
            curr_n_clust = len(set(clusters))
            print(f'RES: {res_var:.5f} | PREVIOUS (max): {previous_n_clust} | CURRENT: {curr_n_clust}')

            # Appends tuple (res, maximum value either current number or previous)
            results.append((res_var, max([curr_n_clust, previous_n_clust])))
            res_list.append(res_var)
            # Appends value of maximum either current number or previous
            nr_clust_list.append(max([curr_n_clust, previous_n_clust]))
            # Updates value of maximum either current or previous
            previous_n_clust = max([curr_n_clust, previous_n_clust])

        # Plot variation of res vs number clusters
        df_results = pd.DataFrame(results)
        fig = sns.scatterplot(x=df_results[0], y=df_results[1])
        plt.xlabel('Resolution')
        plt.ylabel('Max number of clusters')
        plt.savefig(f'{fig_dir}/{exp}_clust_res_var.jpeg')

        # Find a possible resolution
        # Extract all max cluster numbers found
        set_nr_clust = set(nr_clust_list)

        # Count occurrence of each max cluster number
        counts = [nr_clust_list.count(number) for number in set_nr_clust]

        # Make dict of max clust number: counts of occurrences of max
        dict_counts = {}
        for nr_clust in set_nr_clust:
            dict_counts[nr_clust] = nr_clust_list.count(nr_clust)

        # Find max streak
        max_streak = max(counts)

        # Index of max streak initiation
        pos = counts.index(max_streak)

        # Find resolution of same index
        optimum_clust_nr = list(dict_counts.keys())[pos]
        first_occur_clust_nr = nr_clust_list.index(optimum_clust_nr)
        optimum_res = res_list[first_occur_clust_nr]
        res = optimum_res
        print(f'Suggested resolution parameter: {res:.5f}')

    else:
        res = float(res)

    random.seed(1234)
    clusters = Gm.community_leiden(objective_function='CPM',
                                    weights='weight',
                                    resolution=res, # for py=3.8+
                                    #resolution_parameter=res,
                                    n_iterations=-1)

    # Extract series of cluster_id
    clusters = pd.Series(clusters.membership, index=Gm.vs['name']).sort_index().values

    # Add +1 to cluster numbers since count starts with 0
    clusters = clusters + 1
    print(f'TOTAL clusters identified = {len(set(clusters))}')

    # Add cluster_id column to pd peaks df
    df_peaks['cluster_id'] = clusters
    df_peaks = df_peaks.sort_values(by=['cluster_id'])

    return df_peaks, res

def plot_clusters(df_peaks, res, fig_dir, exp):
    """Plots clusters.

    Args:
        df_peaks (dataframe): Dataframe of peaks with embedding and assigned cluster.
        res (float): Resolution parameter.
        fig_dir (path): Path to figures directory.
        exp (string): Name of experiment.
    """
    # Read cluster_id as strings to create categorical coloring
    df_peaks['cluster_id'] = df_peaks['cluster_id'].astype(int)
    df_peaks = df_peaks.sort_values(by='cluster_id')
    df_peaks['cluster_id'] = df_peaks['cluster_id'].astype(str)
    df_peaks = df_peaks.rename(columns={'cluster_id': 'Cluster'})

    colors = cc.glasbey_light

    # Identify ranges for plotting the axes
    min_x = min(df_peaks.x)
    min_y = min(df_peaks.y)
    max_x = max(df_peaks.x)
    max_y = max(df_peaks.y)

    # Make Plotly Express interactive plot
    fig = px.scatter(
                    df_peaks, x='x', y='y',
                    color='Cluster',
                    color_discrete_sequence=colors,
                    hover_data=['name','element','Cluster'],
                    width=940,
                    height=940,
                    labels={'x': 'UMAP 1',
                            'y': 'UMAP 2',},
                    )
    # Change size of data points
    fig.update_traces(marker={'size': 1})

    # Fix the size of figure legend
    fig.update_layout(legend={'itemsizing': 'constant'},
                    plot_bgcolor='white',
                    coloraxis_colorbar_title_text = '', #! DOESN'T WORK
                    xaxis_range=[math.floor(min_x)-0.5, math.ceil(max_x)+0.5],
                    yaxis_range=[math.floor(min_y)-0.5, math.ceil(max_y)+0.5],
                    title=dict(text=f'Clustering (res={res:.6f})', x=0.5, xanchor='center', yanchor='top'))

    fig.update_xaxes(ticks='outside',
                    showline=True,
                    linecolor='black',
                    mirror=True,
                    showgrid=False,
                    zeroline=False)

    fig.update_yaxes(ticks='outside',
                    showline=True,
                    linecolor='black',
                    mirror=True,
                    showgrid=False,
                    zeroline=False)

    fig.write_html(f'{fig_dir}/{exp}_clust_leiden.html')
    fig.write_image(f'{fig_dir}/{exp}_clust_leiden.jpeg')
    fig.write_image(f'{fig_dir}/{exp}_clust_leiden.svg')


if __name__ == '__main__':

    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    res = snakemake.config['res']
    out_dir = exp_dir/'matrix'
    fig_dir = exp_dir/'figures'

    # Load the graph
    print(f'STEP 1 - Loading graph and building igraph object...')
    Gm = get_graph(exp, out_dir)
    print('-----DONE!-----')

    # Find df_peaks with embedding
    df_path = sorted(out_dir.glob(f'*_embed.csv'))[-1]
    df_peaks = pd.read_csv(df_path, header=0, index_col=0)

    # Cluster the data and append info to df peaks
    print(f'STEP 2 - Clustering the data...')
    df_peaks, res = cluster_leiden(Gm, df_peaks, res, fig_dir, exp)
    df_peaks.to_csv(f'{out_dir}/{df_path.stem}_clust.csv')
    print('-----DONE!-----')

    # Plot the clusters
    print(f'STEP 3 - Plot clusters...')
    plot_clusters(df_peaks, res, fig_dir, exp)
    print('-----DONE!-----')
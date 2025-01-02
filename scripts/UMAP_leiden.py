import pandas as pd
import numpy as np
from pathlib import Path, PurePath
import plotly.express as px
import umap
import math

def get_scores(matrix_path):
    """Extracts only scores of each sample across peaks.

    Args:
        matrix_path (Path): Pathlib object of path to annotated matrix.

    Returns:
        dataframe: Annotated matrix.
        dataframe: Matrix of scores only.
    """
    # Read matrix file
    matrix = pd.read_csv(matrix_path,
                        header=0,
                        index_col=0)

    # Drop everything except the scores
    df_scores = matrix.drop(labels=['name',
                                    'chrom',
                                    'start',
                                    'end',
                                    'element',
                                    'family',
                                    'TE'],
                            axis=1)
    # These columns will drop only if human data since kzfp is absent for mouse
    df_scores = df_scores.drop(list(df_scores.filter(regex='kzfp|TSS|transcripts')),
                                axis=1)
    return matrix, df_scores

def get_embed(matrix, df_scores, param_dict, df_metrics, exp, out_dir, matrix_stem):
    """Generates UMAP embedding of data.

    Args:
        matrix (dataframe): Annotated matrix.
        df_scores (dataframe): Matrix of scores only.
        param_dict (dictionary): Dictionary of UMAP parameters.
        df_metrics (dataframe): Saves embeddings of dataset.
        exp (string): Experiment name.
        out_dir (string): Output directory.
        matrix_stem (string): Stem of matrix file.

    Returns:
        dataframe: Updated dataframe of embeddings.
        dataframe: Annotated matrix with embedding for dataset.
    """

    print(f'---Generating the embedding...')
    # Apply the UMAP transformation to all embeddings
    reducer = umap.UMAP(n_neighbors=int(param_dict['n_neighbors']),
                        repulsion_strength=float(param_dict['repulsion_strength']),
                        metric=param_dict['metric'],
                        n_epochs=int(param_dict['n_epochs']),
                        spread=float(param_dict['spread']),
                        min_dist=float(param_dict['min_dist']),
                        random_state=22
                        )
    embedding = reducer.fit_transform(df_scores)
    df_metrics['all'] = [embedding]
    df_embed = pd.DataFrame(embedding, columns=['x','y'])
    print('DONE!')

    print('---Adding embedding to data...')
    df_merged = pd.concat([df_embed, matrix], axis=1)
    df_merged.to_csv(f'{out_dir}/{matrix_stem}_embed.csv')
    print('DONE!')

    print('---Obtaining graph...')
    coo_graph = reducer.graph_.tocoo()
    graph = pd.DataFrame(np.vstack([coo_graph.row, coo_graph.col, coo_graph.data]).T,
                        columns=('source', 'target', 'weight'))
    graph.to_pickle(f'{out_dir}/{exp}_graph_leiden.pkl')
    print('DONE!')

    return df_metrics, df_merged

def plot_umap(df_merged, fig_dir, exp):
    """Plots UMAP for each sample.

    Args:
        df_merged (dataframe): Annotated matrix with embedding for dataset.
        fig_dir (Path): Pathlib object of path to figures directory.
        exp (string): Experiment name.
    """
    # Identify ranges for plotting the axes
    min_x = min(df_merged.x)
    min_y = min(df_merged.y)
    max_x = max(df_merged.x)
    max_y = max(df_merged.y)

    cols = df_merged.columns
    samples = df_merged.columns[cols.get_loc('end')+1 : cols.get_loc('element')]

    for sample in samples:

        fig = px.scatter(
                        df_merged, x='x', y='y',
                        color=sample,
                        color_continuous_scale=[(0, 'rgb(43, 57, 147)'),
                                                (0.25, 'rgb(0, 136, 208)'),
                                                (0.5, 'rgb(0, 185, 155)'),
                                                (0.75, 'rgb(247, 182, 86)'),
                                                (1, 'rgb(253, 246, 63)')],
                        range_color=[0,100],
                        hover_data=['name','element','family','kzfps'],
                        width=750,
                        height=750,
                        labels={'x': 'UMAP 1',
                                'y': 'UMAP 2',},
                        )
        # Change size of data points
        fig.update_traces(marker={'size': 3},
                        showlegend=True)

        fig.update_layout(plot_bgcolor='rgb(176, 193, 209)',
                        coloraxis_colorbar_title_text = '',
                        xaxis_range=[math.floor(min_x)-0.5, math.ceil(max_x)+0.5],
                        yaxis_range=[math.floor(min_y)-0.5, math.ceil(max_y)+0.5],
                        title=dict(text=sample, x=0.5, xanchor='center', yanchor='top'))

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

        fig.write_html(f'{fig_dir}/{exp}_UMAP_{sample}.html')
        fig.write_image(f'{fig_dir}/{exp}_UMAP_{sample}.jpeg')


if __name__ == '__main__':

    genome = snakemake.config['genome']
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    tss_dist = int(snakemake.config['tss_dist'])
    utils_dir = Path.cwd()/'utils'
    bed_dir = exp_dir/'final_peaks'
    out_dir = exp_dir/'matrix'
    fig_dir = exp_dir/'figures'

    matrix_path = sorted(out_dir.glob('*class_kzfps.csv'))[0]
    matrix_stem = matrix_path.stem

    # Create empty dict to save parameters for UMAP
    params = ['metric', 'n_neighbors', 'repulsion_strength', 'min_dist', 'n_epochs', 'spread']
    param_dict = {key: snakemake.config[key] for key in params}

    print('STEP 1 - Retrieving necessary files...')
    matrix, df_scores = get_scores(matrix_path)
    print('-----DONE!-----')

    # Create df in which to save embeddings
    df_metrics = pd.DataFrame(np.nan,
                                index=[0],
                                columns=['all'])
    df_metrics['metric'] = param_dict['metric']

    # Create UMAP object and save embedding in pickle file
    print(f'STEP 2 - Generating embedding...')
    df_metrics, df_merged = get_embed(matrix,
                                    df_scores,
                                    param_dict,
                                    tss_dist,
                                    df_metrics,
                                    exp,
                                    out_dir,
                                    matrix_stem)

    file = f"{out_dir}/{exp}_leiden_{param_dict['metric']}_nneigh{param_dict['n_neighbors']}_repul{param_dict['repulsion_strength']}_mindist{param_dict['min_dist']}_neps{param_dict['n_epochs']}_spr{param_dict['spread']}.pkl"
    df_metrics.to_pickle(file)
    print('-----DONE!-----')

    print(f'STEP 3 - Plotting UMAP...')
    plot_umap(df_merged, fig_dir, exp)
    print('-----DONE!-----')
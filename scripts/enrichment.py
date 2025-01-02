import pybedtools as bt
import pandas as pd
import numpy as np
from pathlib import Path, PurePath
import seaborn as sns
import scipy.stats as stats
import plotly.express as px
import matplotlib
matplotlib.use('Agg')

def avg_per_cluster(clusters, df_clust):
    """Calculates average score of samples per cluster of peaks.

    Args:
        clusters (list): List of cluster ids.
        df_clust (dataframe): Dataframe of scores and cluster ids.
    Returns:
        dataframe: Dataframe of average scores per cluster.
        dataframe: Dataframe of enriched sample names per cluster.
    """
    avg = {}

    for cluster in clusters:

        peaks = df_clust[df_clust['cluster_id']==cluster]

        # Drop the cluster column before calculating the mean of the df columns
        peaks = peaks.drop(labels=['cluster_id'], axis=1)
        values = list(peaks.describe().loc['mean'].values)

        # Add info to dictionary
        avg[f'{cluster}']=values

    # Create df of cluster averages
    df_avg = pd.DataFrame(data=avg).T
    df_avg.columns = peaks.columns

    #* Create doc with assigned enriched samples in cluster
    # Iterate through rows of df to iterate through cluster_ids
    enrich = {}
    for cluster_id in df_avg.index:

        # Add to dict column names which satisfy thresh of avg score
        enrich[cluster_id] = ', '.join(df_avg.iloc[int(cluster_id)-1,:][df_avg.iloc[int(cluster_id)-1,:] > 30].index)

    df_enrich = pd.DataFrame(data=enrich, index=['samples']).T

    return df_avg, df_enrich

def plot_avg(df_avg):
    """Plots clustermap of average scores of samples per cluster.

    Args:
        df_avg (dataframe): Dataframe of average scores per cluster.

    Returns:
        ClusterGrid: Clustermap of average scores per cluster.
    """
    sns.set_theme(font_scale=7)
    heatmap = sns.clustermap(df_avg.T,
                            yticklabels=True,
                            figsize=(90,90),
                            cmap='viridis',
                            vmin=0,
                            vmax=100,
                            metric='correlation', #! check differences in how it plots
                            cbar_pos=(0.02, 0.4, 0.01, 0.1))
    heatmap.ax_row_dendrogram.set_visible(False)
    heatmap.ax_col_dendrogram.set_visible(False)

    return heatmap

def peaks_in_cluster(clusters, matrix, bed_dir, exp, genome_file):
    """Creates BED file of peaks found in each cluster and dict of peak.

    Args:
        clusters (list): List of cluster ids.
        matrix (dataframe): Clustered matrix.
        bed_dir (Path): Path to folder with final peaks.
        exp (string): Name of experiment.
        genome_file (Path): Path to genome domain file.

    Returns:
        dictionary: Dictionary of clusters and BedTool of peaks found in them.
    """
    # Create dict to store peaks for each cluster
    cluster_peaks = {}

    for cluster in ['all'] + clusters:

        if cluster == 'all':
            merged_peaks = bt.BedTool.from_dataframe(
                            matrix[['chrom', 'start', 'end', 'name']]
                                                    ).sort(g=genome_file)

        else:
            merged_peaks = bt.BedTool.from_dataframe(
                matrix[matrix['cluster_id'] == cluster][["chrom", "start", "end", "name"]]
                ).sort(g=genome_file)

        # Save BED file of peaks for each cluster
        merged_peaks.saveas(f'{bed_dir}/{exp}_merged_peaks_cluster_{cluster}.bed')

        # Store peaks for cluster
        #merged_peaks = merged_peaks.each(midpoint).saveas()
        cluster_peaks[cluster] = merged_peaks

    return cluster_peaks

def tes_in_cluster(clusters, matrix, te_subfamilies):
    """Creates dictionary of TE subfamilies that come up in each cluster.

    Args:
        clusters (list): List of cluster ids.
        matrix (dataframe): Clustered matrix.
        te_subfamilies (set): Set of TE subfamilies present in dataset.

    Returns:
        dictionary: Dictionary of TE subfamilies present in each cluster.
    """
    # Create empty dict to store cluster nr and subfamilies
    te_per_cluster = {}

    for cluster in ['all'] + clusters:
        if cluster == 'all':
            te_per_cluster[cluster] = te_subfamilies
        else:
            # Get 'element' column and exclude those not-overlapping TEs which lead to None
            tes_in_df = list(matrix.element[(matrix.cluster_id == cluster) & (~matrix.element.isnull()) & (matrix.element != 'None')])

            # Deal with cases where peak overlaps multiple elements
            tes_extend_in_df = []
            tes_extend_in_df.extend(subfam for multi_subfam in tes_in_df for subfam in multi_subfam.split(','))

            # Make set from list of all subfamilies
            te_per_cluster[cluster] = set(tes_extend_in_df)

    return te_per_cluster

def label_te(te, metadata, column):
    """Extracts the TE metadata for a specific TE subfamily and adds it to a dataframe.

    Args:
        te (string): TE subfamily name.
        metadata (datafram): Dataframe of TE information.
        column (string): Name of column of information wanted.

    Returns:
        list: List of content to be filled into dataframe.
    """
    return metadata.loc[metadata["name"] == te, column].tolist()

def te_fisher(cluster, te_df, merged_peaks, te_per_cluster, te_metadata, genome_file):
    """Calculates Fisher enrichment of TEs per cluster.

    Args:
        cluster (string): Cluster id.
        te_df (dataframe): Dataframe of all TEs in domain.
        merged_peaks (BedTool): BedTool object of peaks in cluster.
        te_per_cluster (dictionary): Dictionary of TE subfamilies present in each cluster.
        te_metadata (dataframe): Dataframe of TE information.
        genome_file (Path): Path to genome domain file.

    Returns:
        dataframe: Dataframe of enrichment scores for TEs in cluster.
    """
    fisher_list_te = []

    # Calculate enrichment of TE subfamilies present in cluster only
    for family in te_per_cluster[cluster]:

        target_te_df = te_df.loc[te_df['name'] == family][["chrom", "start", "end"]]
        target_te_bed = bt.BedTool.from_dataframe(target_te_df).sort(g=genome_file).saveas()
        n_genome = len(target_te_bed)

        fisher = merged_peaks.fisher(b=target_te_bed, g=genome_file, u=True)
        overlaps = int(fisher.text.split("\n")[2].split(" ")[-1])

        # Extract the counts from the bt Fisher
        list_1 = fisher.text.split('\n')[8].split('|')
        row_1 = [s.strip() for s in list_1][1:3]

        list_2 = fisher.text.split('\n')[9].split('|')
        row_2 = [s.strip() for s in list_2][1:3]

        # Use fisher from stats as bt has overflow errors
        _ , p_val = stats.fisher_exact([row_1, row_2], alternative='greater')

        if p_val == 0:
            p_score = np.inf
        else:
            p_score = -np.log10(p_val)

        fisher_list_te.append([family, p_val, p_score, overlaps, n_genome])

    fisher_df_te = pd.DataFrame.from_records(fisher_list_te, columns=["family", "p_val", "p_score", "overlaps", "n_genome"])
    fisher_df_te["p_score"] = fisher_df_te["p_score"].apply(lambda x: 320 if x == np.inf else x)
    fisher_df_sorted_te = fisher_df_te.sort_values(by='p_score', ascending=False)

    # Extracts the family of the TE based on subfamily
    fisher_df_sorted_te['type'] = fisher_df_sorted_te['family'].apply(lambda x: label_te(x, te_metadata, 'type'))

    return fisher_df_sorted_te

def plot_te(cluster, top_te, enrich_thresh_te, fisher_df_sorted_te):
    """Plots top enriched TEs in cluster.

    Args:
        cluster (string): Cluster id.
        top_te (int): Number of top enriched TEs to plot.
        enrich_thresh_te (float): Threshold of -log10(p) for enrichment of TEs.
        fisher_df_sorted_te (dataframe): Dataframe of enrichment scores for TEs in cluster.

    Returns:
        figure: Barplot of top enriched TEs.
    """
    cluster_top_tes = fisher_df_sorted_te.head(top_te).copy()
    fig = px.bar(cluster_top_tes,
                x=list(cluster_top_tes['family']),
                y=list(cluster_top_tes['p_score']),
                color_discrete_sequence=['darkcyan'],
                hover_data=['overlaps', 'n_genome'],
                title=f'<b>CLUSTER #{cluster} | Top {top_te} TEs</b>')
    fig.add_hline(y=enrich_thresh_te, line_dash='dash', line_color='black')
    fig.update_layout(font_family='Tahoma',
                    font_size=16,
                    title_font_size=24,
                    xaxis_title='TE subfamily',
                    yaxis_title='-log10(p)',
                    showlegend=False)
    fig.update_xaxes(tickangle=-45)

    return fig

def plot_clustmap(cluster_rank, index, out_dir, exp, enrich_thresh, types):
    """Plots clustermap of TEs/KZFPs enriched above threshold.

    Args:
        cluster_rank (dictionary): Dictionary of TEs/KZFPs above threshold [clust_id: peaks/TE/KZFP].
        index (set): Set of TE subfamilies/KZFPs present in dataset.
        out_dir (Path): Path to output directory.
        exp (string): Name of experiment.
        enrich_thresh (float): Threshold of -log10(p) for enrichment.
        types (string): 'te' or 'kzfps' depending on target.

    Returns:
        dataframe: Dataframe of full enrichment scores.
        dataframe: Dataframe of enrichment scores without clusters of no enrichment.
        ClusterGrid: Clustermap of enriched TEs/KZFPs.
    """
    list_index = list(index)
    enrichment_scores = pd.DataFrame(columns=cluster_rank.keys(), index=list_index)

    # Iterate through dictionary contents to extract values
    for cluster, rank in cluster_rank.items():

        # Eg: for i in te_subfamilies
        for i in index:
            series = rank.loc[rank.family == i, 'p_score']

            if len(series) != 0:
                value = series.item()
            else:
                fisher_sorted_df = pd.read_csv(f'{out_dir}/{exp}_{types}_rank_cluster_{cluster}.csv',
                                                header=0, index_col=0)
                not_sig_series = fisher_sorted_df.loc[fisher_sorted_df.family == i, 'p_score']
                value = not_sig_series.item() if len(not_sig_series) != 0 else 0

            enrichment_scores.loc[i, cluster] = float(value)

    enrichment_scores = enrichment_scores.apply(pd.to_numeric)

    # Remove 'all' calculations and TEs only enriched in 'all']
    enrichment_scores_not_all = enrichment_scores.drop(columns=['all'])
    enrichment_scores_not_all['max'] = enrichment_scores_not_all.max(axis=1)
    enrichment_scores_not_all = enrichment_scores_not_all[enrichment_scores_not_all['max'] > enrich_thresh]
    enrichment_scores_not_all = enrichment_scores_not_all.drop(columns=['max'])

    # Remove all cluster columns where nothing is enriched
    enrichment_scores_not_all = enrichment_scores_not_all[enrichment_scores_not_all.columns[enrichment_scores_not_all.max() > enrich_thresh]]

    enrichment_scores_not_all = enrichment_scores_not_all.apply(pd.to_numeric)
    print(f"Total elements: {len(enrichment_scores_not_all.index)}")
    sns.set_theme(font_scale=4.0)
    plot_clustmap = sns.clustermap(enrichment_scores_not_all,
                            yticklabels=True,
                            figsize=(40, 60),
                            cmap='viridis',
                            vmin=0,
                            vmax=100,
                            cbar_pos=(0.02, 0.4, 0.01, 0.1))
    plot_clustmap.ax_row_dendrogram.set_visible(False)
    plot_clustmap.ax_col_dendrogram.set_visible(False)

    return enrichment_scores, enrichment_scores_not_all, plot_clustmap

def kzfps_in_cluster(clusters, matrix, kzfp_present):
    """Creates a dictionary of KZFPs which target TE subfamilies in each cluster.

    Args:
        clusters (list): List of cluster ids.
        matrix (dataframe): Clustered matrix.
        kzfp_present (set): Set of KZFPs binding in dataset.

    Returns:
        dictionary: Dictionary of KZFPs binding in each cluster.
    """
    # Create empty dict to store cluster nr and KZFPs
    kzfp_per_cluster = {}

    for cluster in ['all'] + clusters:
        if cluster == 'all':
            kzfp_per_cluster[cluster] = kzfp_present
        else:
            # Get 'kzfps' column and exclude those not-bound which lead to None
            kzfps_in_df = list(matrix.kzfps[(matrix.cluster_id == cluster) & (~matrix.kzfps.isnull()) & (matrix.kzfps != 'None')])

            # Deal with cases where peak is bound by multiple KZFPs
            kzfps_extend_in_df = []
            kzfps_extend_in_df.extend(gene for multi_gene in kzfps_in_df for gene in multi_gene.split(','))

            # Make set from list of all kzfps
            kzfp_per_cluster[cluster] = set(kzfps_extend_in_df)

    return kzfp_per_cluster

def label_kzfp(kzfp, metadata, column):
    """Extracts the KZFP metadata for a specific KZFP and adds it to a dataframe.

    Args:
        kzfp (string): KZFP name.
        metadata (dataframe): Dataframe of KZFP information.
        column (string): Name of column of information wanted.

    Returns:
        list: List of content to be filled into dataframe.
    """    
    return metadata.loc[metadata["Gene symbol"] == kzfp, column].item()

def kzfp_fisher(cluster, kzfp_df, merged_peaks, kzfp_per_cluster, kzfp_metadata, genome_file):
    """Calculates Fisher enrichment of KZFPs per cluster.

    Args:
        cluster (string): Cluster id.
        kzfp_df (dataframe): Dataframe of KZFP binding sites in domain.
        merged_peaks (BedTool): BedTool object of peaks in cluster.
        kzfp_per_cluster (dictionary): Dictionary of KZFPs present in each cluster.
        kzfp_metadata (dataframe): Dataframe of KZFP information.
        genome_file (Path): Path to genome domain file.

    Returns:
        dataframe: Dataframe of enrichment scores for KZFP binding sites in cluster.
    """    
    fisher_list_kzfp = []

    # Calculate enrichment of KZFP present in cluster only
    for kzfp in kzfp_per_cluster[cluster]:

        target_kzfp_df = kzfp_df.loc[kzfp_df['name'] == kzfp][["chrom", "start", "end"]]
        target_kzfp_bed = bt.BedTool.from_dataframe(target_kzfp_df).sort(g=genome_file).saveas()
        n_genome = len(target_kzfp_bed)

        fisher = merged_peaks.fisher(b=target_kzfp_bed, g=genome_file, u=True)
        overlaps = int(fisher.text.split("\n")[2].split(" ")[-1])

        # Extract the counts from the bt Fisher
        list_1 = fisher.text.split('\n')[8].split('|')
        row_1 = [s.strip() for s in list_1][1:3]

        list_2 = fisher.text.split('\n')[9].split('|')
        row_2 = [s.strip() for s in list_2][1:3]

        # Use fisher from stats as bt has overflow errors
        _ , p_val = stats.fisher_exact([row_1, row_2], alternative='greater')

        if p_val == 0:
            p_score = np.inf
        else:
            p_score = -np.log10(p_val)

        fisher_list_kzfp.append([kzfp, p_val, p_score, overlaps, n_genome])

    fisher_df_kzfp = pd.DataFrame.from_records(fisher_list_kzfp, columns=["family", "p_val", "p_score", "overlaps", "n_genome"])
    fisher_df_kzfp["p_score"] = fisher_df_kzfp["p_score"].apply(lambda x: 320 if x == np.inf else x)
    fisher_df_sorted_kzfp = fisher_df_kzfp.sort_values(by='p_score', ascending=False)

    fisher_df_sorted_kzfp['domain'] = fisher_df_sorted_kzfp['family'].apply(lambda x: label_kzfp(x, kzfp_metadata, 'KRAB domain configuration'))
    fisher_df_sorted_kzfp['age'] = fisher_df_sorted_kzfp['family'].apply(lambda x: label_kzfp(x, kzfp_metadata, 'DNA binding homolog found at the latest'))
    fisher_df_sorted_kzfp['gene_id'] = fisher_df_sorted_kzfp['family'].apply(lambda x: label_kzfp(x, kzfp_metadata, 'Ensembl ID'))

    return fisher_df_sorted_kzfp

def plot_kzfp(cluster, top_kzfp, enrich_thresh_kzfp, fisher_df_sorted_kzfp):
    """Plots top enriched KZFPs in cluster.

    Args:
        cluster (string): Cluster id.
        top_kzfp (int): Number of top enriched KZFPs to plot.
        enrich_thresh_kzfp (float): Threshold of -log10(p) for enrichment of KZFPs.
        fisher_df_sorted_kzfp (dataframe): Dataframe of enrichment scores for KZFPs in cluster.

    Returns:
        figure: Barplot of top enriched KZFPs.
    """    
    cluster_top_kzfps = fisher_df_sorted_kzfp.head(top_kzfp).copy()
    fig = px.bar(cluster_top_kzfps,
                x=list(cluster_top_kzfps['family']),
                y=list(cluster_top_kzfps['p_score']),
                color_discrete_sequence=['cornflowerblue'],
                hover_data=['overlaps', 'n_genome'],
                title=f'<b>CLUSTER #{cluster} | Top {top_kzfp} KZFPs</b>')
    fig.add_hline(y=enrich_thresh_kzfp, line_dash='dash', line_color='black')
    fig.update_layout(font_family='Tahoma',
                    font_size=16,
                    title_font_size=24,
                    xaxis_title='KZFP',
                    yaxis_title='-log10(p)',
                    showlegend=False)
    fig.update_xaxes(tickangle=-45)

    return fig


if __name__ == '__main__':

    genome = snakemake.config['genome']
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    tss_dist = int(snakemake.config['tss_dist'])
    top_te = int(snakemake.config['top_te'])
    top_kzfp = int(snakemake.config['top_kzfp'])
    enrich_thresh_te = float(snakemake.config['enrich_thresh_te'])
    enrich_thresh_kzfp = float(snakemake.config['enrich_thresh_kzfp'])

    out_dir = exp_dir/'matrix'
    fig_dir = exp_dir/'figures'
    bed_dir = exp_dir/'final_peaks'
    utils_dir = Path.cwd()/'utils'

    kzfp_bed = bt.BedTool(f'{utils_dir}/{genome}_kzfps_combined.bed')
    kzfp_metadata = pd.read_csv(f'{utils_dir}/kzfp_info.csv')
    te_metadata = pd.read_csv(f'{utils_dir}/te_info.csv')
    genome_file = f'{utils_dir}/{genome}_domain.txt'
    domain = bt.BedTool(f'{utils_dir}/{genome}_domain.bed')
    te_bed = bt.BedTool(f'{utils_dir}/{genome}_repeats_sorted_filtered.bed')

    #* Check if TE intersection with domain exists
    te_in_domain = utils_dir/f'{genome}_te_in_domain.bed'

    if te_in_domain.is_file():
        # Read existing file
        te_in_domain = bt.BedTool(f'{utils_dir}/{genome}_te_in_domain.bed')

    else:
        # Creates bed object with all the TE intersections in the genome
        te_in_domain = te_bed.intersect(domain, u=True).saveas()

        # Sorts the file based on the genome file (chr are in descending length order)
        te_in_domain = te_in_domain.sort(g=genome_file).saveas()

        # Save bed for future uses
        te_in_domain.saveas(f'{utils_dir}/{genome}_te_in_domain.bed')

    # Turn into df for easier use
    te_df = te_in_domain.to_dataframe()

    #* Check if kzfp binding sites intersection with domain exists
    kzfp_in_domain = utils_dir/f'{genome}_kzfp_in_domain.bed'

    if kzfp_in_domain.is_file():
        # Read existing file
        kzfp_in_domain = bt.BedTool(f'{utils_dir}/{genome}_kzfp_in_domain.bed')

    else:
        # Creates bed object with all the kzfp binding site intersections in the genome
        kzfp_in_domain = kzfp_bed.intersect(domain, u=True).saveas()

        # Sorts the file based on the genome file (chr are in descending length order)
        kzfp_in_domain = kzfp_in_domain.sort(g=genome_file).saveas()

        # Save bed for future uses
        kzfp_in_domain.saveas(f'{utils_dir}/{genome}_kzfp_in_domain.bed')

    # Turn into df for easier use
    kzfp_df = kzfp_in_domain.to_dataframe()

    print('STEP 1 - Retrieving necessary files...')

    matrix_path = sorted(out_dir.glob(f'*_clust.csv'))[0]
    matrix_stem = matrix_path.stem

    matrix = pd.read_csv(matrix_path, header=0, index_col=0)
    df_clust = matrix.drop(labels=['x',
                                    'y',
                                    'name',
                                    'chrom',
                                    'start',
                                    'end',
                                    'element',
                                    'family',
                                    'TE',
                                    f'transcripts_{tss_dist}',
                                    f'TSS_{tss_dist}',
                                    'kzfps',
                                    'countkzfps',
                                    'kzfp_geneid',
                                    'kzfp_domain',
                                    'kzfp_age'],
                            axis=1)

    clusters = sorted(list(df_clust['cluster_id'].unique()))
    print('-----DONE!-----')

    print('STEP 2 - Calculating average scores of samples per cluster...')
    df_avg, df_enrich = avg_per_cluster(clusters, df_clust)
    df_avg.to_csv(f'{out_dir}/{matrix_stem}_avg.csv')
    df_enrich.to_csv(f'{out_dir}/{matrix_stem}_enriched_samples.csv')
    print('-----DONE!-----')

    print('STEP 3 - Plotting average scores of samples per cluster...')
    score_clustmap = plot_avg(df_avg)
    score_clustmap.savefig(f'{fig_dir}/avg_score_heatmap.jpeg')
    print('-----DONE!-----')

    print('STEP 4 - Creating BED files of peaks for each cluster...')
    cluster_peaks = peaks_in_cluster(clusters, matrix, bed_dir, exp, genome_file)
    print('-----DONE!-----')

    print('STEP 5 - Extracting all TE subfamilies which appear in dataset...')
    # Create a list of all entries in the element column of df
    te_from_matrix = list(matrix.element[(~matrix.element.isnull()) & (matrix.element != 'None')])

    # Deal with cases where peak overlaps multiple elements
    te_from_matrix_full = []
    te_from_matrix_full.extend(subfam for multi_subfam in te_from_matrix for subfam in multi_subfam.split(','))
    te_subfamilies = set(te_from_matrix_full)
    print('-----DONE!-----')

    print('STEP 6 - Extracting TE subfamilies for each cluster...')
    te_per_cluster = tes_in_cluster(clusters, matrix, te_subfamilies)
    print('-----DONE!-----')

    # Create dictionary to store TEs enriched above a threshold
    cluster_ranked_tes = {}
    # Expand the dictionary of cluster and peaks
    for cluster, merged_peaks in cluster_peaks.items():

        print(f'STEP 7.{cluster} - Calculating the TE Fisher enrichment for CLUSTER #{cluster}...')
        fisher_df_sorted_te = te_fisher(cluster, te_df, merged_peaks, te_per_cluster, te_metadata, genome_file)
        fisher_df_sorted_te.to_csv(f'{out_dir}/{exp}_te_rank_cluster_{cluster}.csv')
        print('-----DONE!-----')

        print(f'STEP 8.{cluster} - Plotting top {top_te} TEs in CLUSTER #{cluster}...')
        fig = plot_te(cluster, top_te, enrich_thresh_te, fisher_df_sorted_te)
        fig.write_html(f'{fig_dir}/{exp}_te_enrich_cluster_{cluster}_thresh_{enrich_thresh_te}.html')
        print('-----DONE!-----')

        # Saving the TEs above enrichment threshold for every cluster
        passed = fisher_df_sorted_te[fisher_df_sorted_te['p_score'] > enrich_thresh_te]
        cluster_ranked_tes[cluster] = passed

    print('STEP 9 - Plotting the heatmap of enriched TEs...')
    te_enrichment_scores, te_enrichment_scores_not_all, te_clustmap = plot_clustmap(cluster_ranked_tes, te_subfamilies, out_dir, exp, enrich_thresh_te, types='te')
    te_enrichment_scores.to_csv(f'{out_dir}/{exp}_te_enrich_clustmap_scores_thresh_{enrich_thresh_te}_full.csv')
    te_enrichment_scores_not_all.to_csv(f'{out_dir}/{exp}_te_enrich_clustmap_scores_thresh_{enrich_thresh_te}_small.csv')
    te_clustmap.savefig(f'{fig_dir}/{exp}_te_clustmap_thresh_{enrich_thresh_te}.jpeg')
    print('-----DONE!-----')

    print('STEP 10 - Extracting all KZFP which appear in dataset...')
    # Create a list of all entries in the element column of df
    kzfp_from_matrix = list(matrix.kzfps[(~matrix.kzfps.isnull()) & (matrix.kzfps != 'None')])

    # Deal with cases where peak is bound by multiple KZFPs
    kzfp_from_matrix_full = []
    kzfp_from_matrix_full.extend(gene for multi_gene in kzfp_from_matrix for gene in multi_gene.split(','))
    kzfp_present = set(kzfp_from_matrix_full)
    print('-----DONE!-----')

    print('STEP 11 - Extracting KZFPs for each cluster...')
    kzfp_per_cluster = kzfps_in_cluster(clusters, matrix, kzfp_present)
    print('-----DONE!-----')

    # Create dictionary to store KZFPs enriched above a threshold
    cluster_ranked_kzfps = {}
    # Expand the dictionary of cluster and peaks
    for cluster, merged_peaks in cluster_peaks.items():

        print(f'STEP 12.{cluster} - Calculating the KZFP Fisher enrichment for CLUSTER #{cluster}...')
        fisher_df_sorted_kzfp = kzfp_fisher(cluster, kzfp_df, merged_peaks, kzfp_per_cluster, kzfp_metadata, genome_file)
        fisher_df_sorted_kzfp.to_csv(f'{out_dir}/{exp}_kzfp_rank_cluster_{cluster}.csv')
        print('-----DONE!-----')

        print(f'STEP 13.{cluster} - Plotting top {top_kzfp} KZFPs in CLUSTER #{cluster}...')
        fig = plot_kzfp(cluster, top_kzfp, enrich_thresh_kzfp, fisher_df_sorted_kzfp)
        fig.write_html(f'{fig_dir}/{exp}_kzfp_enrich_cluster_{cluster}_thresh_{enrich_thresh_kzfp}.html')
        print('-----DONE!-----')

        # Saving the KZFPs above enrichment threshold for every cluster
        passed = fisher_df_sorted_kzfp[fisher_df_sorted_kzfp['p_score'] > enrich_thresh_kzfp]
        cluster_ranked_kzfps[cluster] = passed

    print('STEP 14 - Plotting the clustmap of enriched KZFPs...')
    kzfp_enrichment_scores, kzfp_enrichment_scores_not_all, kzfp_heatmap = plot_clustmap(cluster_ranked_kzfps, kzfp_present, out_dir, exp, enrich_thresh_kzfp, types='kzfp')
    kzfp_enrichment_scores.to_csv(f'{out_dir}/{exp}_kzfp_enrich_clustmap_scores_thresh_{enrich_thresh_kzfp}_full.csv')
    kzfp_enrichment_scores_not_all.to_csv(f'{out_dir}/{exp}_kzfp_enrich_clustmap_scores_thresh_{enrich_thresh_kzfp}_small.csv')
    kzfp_heatmap.savefig(f'{fig_dir}/{exp}_kzfp_clustmap_thresh_{enrich_thresh_kzfp}.jpeg')
    print('-----DONE!-----')
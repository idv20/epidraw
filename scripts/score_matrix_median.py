import numpy as np
from pathlib import Path, PurePath
import pybedtools
import pandas as pd
from scipy.stats.mstats import gmean

def combine_counts(out_dir, merged_peaks):
    """Combines read counts of all samples into one dataframe.

    Args:
        out_dir (path): Path to matrix output directory.
        merged_peaks (bedtool): BedTool object containing merged peaks.

    Returns:
        dataframe: Read count matrix across all samples.
    """
    counts_files = [counts_file for counts_file in sorted(out_dir.glob('*_read_counts.csv'))]

    count_matrix = merged_peaks.to_dataframe()
    count_matrix = count_matrix.drop(labels=['score', 'strand'], axis=1)

    # Rearrange column order by moving last as first
    cols = count_matrix.columns.tolist()
    cols = cols[-1:] + cols[:-1]

    count_matrix = count_matrix[cols]

    for df in counts_files:

        df_sample = pd.read_csv(df, header=0, index_col=0)
        sample_name = str(df_sample.columns[-1])
        count_matrix = pd.concat([count_matrix,
                                    df_sample[sample_name]], axis=1)

    return count_matrix

def normalise_matrix(count_matrix):
    """Normalise read count matrix using median of ratios method (see DESeq2).

    Args:
        count_matrix (dataframe): Read count matrix across all samples.

    Returns:
        dataframe: Normalised read count matrix across all samples.
    """
    print(count_matrix.head())
    # Transform count matrix to array and add 1 to fix count=0 issues
    counts = np.array(count_matrix.iloc[:,4:]) + 1

    # Calculate geometric mean for each peak
    geom_means = gmean(counts, axis=1)

    # Calculate ratio between counts and geometric mean for each peak
    ratio_sample_ref = np.empty(counts.shape)

    for i in range(len(counts)):
        ratio_sample_ref[i] = np.divide(counts[i],geom_means[i])

    # Calculate median of ratios per sample and divide each raw count by median
    norm_factors = np.empty(counts.shape[1])
    counts_T = np.transpose(counts)
    norm_counts_T = np.empty(counts_T.shape)
    ratio_sample_ref_T = np.transpose(ratio_sample_ref)

    for j in range(len(norm_factors)):
        norm_factors[j] = np.median(ratio_sample_ref_T[j])
        norm_counts_T[j] = np.divide(counts_T[j], norm_factors[j])

    norm_counts = np.transpose(norm_counts_T)

    norm_counts = pd.DataFrame(norm_counts, columns=count_matrix.columns[4:])
    norm_matrix = pd.concat([count_matrix.reset_index(drop=True).iloc[:,0:4],
                            norm_counts.reset_index(drop=True)], axis=1)

    return norm_matrix

if __name__ == '__main__':

    # Extract from snakemake config
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    bam_dir_name = snakemake.config['bam_dir_name']
    bam_dir = exp_dir/bam_dir_name
    bed_dir_name = 'final_peaks'
    bed_dir = exp_dir/bed_dir_name
    out_dir = exp_dir/'matrix'
    merged_peaks = pybedtools.BedTool(sorted(out_dir.glob(f'{exp}_merged_peaks.bed'))[0])

    print('STEP 1 - Combining score matrices...')
    count_matrix = combine_counts(out_dir, merged_peaks)
    count_matrix.to_csv(f'{out_dir}/{exp}_count_matrix.csv')
    print('-----DONE!-----')

    print('STEP 2 - Normalising matrix...')
    norm_matrix = normalise_matrix(count_matrix)
    norm_matrix.to_csv(f'{out_dir}/{exp}_score_matrix_norm.csv')
    print('-----DONE!-----')
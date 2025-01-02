import pandas as pd
import pybedtools as bt
from pathlib import Path, PurePath

def store_info(entry, with_text=False):
    """Helper function which stores info from intersections.

    Args:
        entry (pd.array): Dataframe array.
        with_text (bool, optional): Keep or not entry text. Defaults to False.

    Returns:
        int OR string: 0,1 OR 'None',entry.
    """
    # If no overlap with TEs
    if entry == '.':
        # To mark 0 in TE column for absence
        if not with_text:
            return 0
        # To replace . with None in element/family columns for absence
        else:
            return 'None'

    # If overlap with TEs
    else:
        # To mark 1 in TE column for presence
        if not with_text:
            return 1
        # To keep element/family values where presence
        else:
            return entry

def intersect_te(matrix, bt_repeats):
    """Intersects peaks with repeats.

    Args:
        matrix (dataframe): Dataframe of normalised peaks.
        bt_repeats (BedTool): BED file of annotated repeats.

    Returns:
        dataframe: Dataframe with TE intersection info.
        bedTool: BedTool of peak info [chrom, start, end].
    """
    # Extract peaks from matrix
    peaks = matrix[['chrom', 'start', 'end']]
    bt_peaks = bt.BedTool.from_dataframe(peaks).sort()
    print(bt_peaks.to_dataframe())

    # Intersect peaks with TEs
    bt_int_rep = bt_peaks.intersect(bt_repeats, wao=True, f=0.5, F=0.5, e=True)
    bt_int_rep_merged = bt_int_rep.merge(c='7,10', o='collapse,collapse')

    df_int_rep = bt_int_rep_merged.to_dataframe()
    df_int_rep.columns = ['chrom', 'start', 'end', 'element', 'family']

    df_int_rep['TE'] = df_int_rep['element'].apply(store_info)
    df_int_rep['element'] = df_int_rep['element'].apply(store_info, with_text=True)
    df_int_rep['family'] = df_int_rep['family'].apply(store_info, with_text=True)

    return df_int_rep, bt_peaks

def intersect_tss(bt_tss, genome_path, bt_peaks, tss_dist):
    """Intersects peaks with TSS found within user-defined distance of eachother.

    Args:
        bt_tss (BedTool): BedTool of TSSs.
        genome_path (string): Full path of genome .txt file containing chr domains.
        bt_peaks (BedTool): BedTool of peaks in dataset [chrom, start, end].
        tss_dist (int): Maximum distance from TSS.

    Returns:
        dataframe: Dataframe with peak info and TSS annotations.
    """
    # Extend the tss by ± tss_dist
    bt_tss_dist = bt_tss.slop(g=genome_path, b=tss_dist).saveas()

    # Intersect peaks with tss ± tss_dist
    bt_int_tss_dist = bt_peaks.intersect(bt_tss_dist, wao=True, f=0.5, F=0.5, e=True)
    bt_int_tss_dist_merged = bt_int_tss_dist.merge(c='7', o='distinct')

    # Convert bt to df
    df_int_tss_dist = bt_int_tss_dist_merged.to_dataframe()

    # Check if any TSS or not and add to peak df
    df_int_tss_dist.columns = ['chrom', 'start', 'end', f'transcripts_{tss_dist}']

    df_int_tss_dist[f'TSS_{tss_dist}'] = df_int_tss_dist[f'transcripts_{tss_dist}'].apply(store_info)
    df_int_tss_dist[f'transcripts_{tss_dist}'] = df_int_tss_dist[f'transcripts_{tss_dist}'].apply(
                                            store_info, with_text=True)

    return df_int_tss_dist

def te_matrix(matrix, df_int_rep, tss_dist=None, df_int_tss_dist=None):
    """Combine all information into one matrix.

    Args:
        matrix (dataframe): Original normalised matrix.
        df_int_rep (dataframe): Dataframe with TE information.
        tss_dist (int): Maximum distance from TSS. Defaults to None.
        df_int_tss_dist (dataframe, optional): Dataframe with peak info and TSS annotations. Defaults to None.

    Returns:
        dataframe: Concatenated dataframe with peaks and TE info (and TSS if human data analysed).
    """
    # Check if TSS info present ⟹ human data
    if df_int_tss_dist is not None:

        df_te = pd.concat([matrix,
                            df_int_rep[['element', 'family', 'TE']],
                            df_int_tss_dist[[f'transcripts_{tss_dist}', f'TSS_{tss_dist}']]], axis=1)

    # If no TSS info present ⟹ mouse data
    else:
        df_te = pd.concat([matrix,
                            df_int_rep[['element', 'family', 'TE']]], axis=1)

    return df_te

def label_kzfps(kzfps, metadata, column):
    """Extracts KZFP information from metadata for each ZNF whose binding site overlaps a peak.

    Args:
        kzfps (array): Array of KZFPs whose BS overlaps with a peak.
        metadata (dataframe): Dataframe of KZFP metadata.
        column (string): Name of column from metadata dataframe.

    Returns:
        list: List of metadata information separated by commas.
    """
    # Returns 'None' if no KZFP overlap
    if kzfps == 'None':
        return 'None'

    else:
        # Splits the string of multiple KZFPs (if multi-hits)
        kzfps_list = kzfps.split(',')
        label_list = []

        for kzfp in kzfps_list:
            label_list.append(metadata.loc[metadata['Gene symbol'] == kzfp, column].item())

        return ','.join(label_list)

def intersect_kzfp(bt_peaks, bt_kzfp, kzfp_meta):
    """Intersects KZFP binding site with peaks.

    Args:
        bt_peaks (BedTool): BedTool of peak info [chrom, start, end]
        bt_kzfp (BedTool): BedTool of KZFP binding sites.
        kzfp_meta (dataframe): Dataframe of KZFP metadata.

    Returns:
        dataframe: Peak info and KZFP annotation with metadata.
    """
    # Intersect peaks with KZFP
    bt_int = bt_peaks.intersect(bt_kzfp, wao=True)
    bt_int_merged = bt_int.merge(c='7', o='collapse')

    df_int = bt_int_merged.to_dataframe()

    # Store name of KZFPs which bind and rename col
    df_int['name'] = df_int['name'].apply(store_info, with_text=True)
    df_int = df_int.rename(columns={'name':'kzfps'})

    # Count KZFPs and add gene_id, domain and age info
    df_int['countkzfps'] = df_int.apply(lambda row:
                                        len(str(row.kzfps).split(','))
                                        if not row.kzfps == 'None'
                                        else 0,
                                        axis=1)

    df_int['kzfp_geneid'] = df_int.apply(lambda row:
                                        label_kzfps(row.kzfps,
                                                    kzfp_meta,
                                                    'Ensembl ID'),
                                        axis=1)

    df_int['kzfp_domain'] = df_int.apply(lambda row:
                                        label_kzfps(row.kzfps,
                                                    kzfp_meta,
                                                    'KRAB domain configuration'),
                                        axis=1)

    df_int['kzfp_age'] = df_int.apply(lambda row:
                                    label_kzfps(row.kzfps,
                                                kzfp_meta,
                                                'DNA binding homolog found at the latest'),
                                    axis=1)

    return df_int

#   function: to ouput final files
def final_matrix(matrix, df_int):
    """Concatenates everything

    Args:
        matrix (dataframe): Peak info and TE annotation with/without TSS.
        df_int (dataframe): Intersection with KZFPs.

    Returns:
        dataframe: Final concatenated dataframe.
        dataframe: Peak info and KZFP info only.
    """
    df_final = pd.concat([matrix,
                        df_int[['kzfps', 'countkzfps', 'kzfp_geneid',
                                'kzfp_domain', 'kzfp_age']]],
                        axis=1)

    df_znflist = df_int.drop(['countkzfps'], axis=1)

    return df_final, df_znflist


if __name__ == '__main__':

    genome = snakemake.config['genome']
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    tss_dist = int(snakemake.config['tss_dist'])
    utils_dir = Path.cwd()/'utils'
    out_dir = exp_dir/'matrix'
    matrix_name = f'{exp}_score_matrix_norm.csv'

    # Read the matrix file
    matrix = pd.read_csv(out_dir/matrix_name,
                        header=0,
                        index_col=0)

    print('STEP 1 - Reading files...')
    # Load the RepeatMasker bedfile
    bt_repeats = bt.BedTool(f'{utils_dir}/{genome}_repeats_sorted_filtered.bed')

    # Load the genome size file
    genome_path = f'{utils_dir}/{genome}_domain.txt'
    print('-----DONE!-----')

    print('STEP 2 - Intersecting with TEs...')
    df_int_rep, bt_peaks = intersect_te(matrix, bt_repeats)
    print('-----DONE!-----')

    print('STEP 3 - Intersecting with TSS...')
    bt_tss = bt.BedTool(f'{utils_dir}/{genome}_tss_full.bed')
    df_int_tss_dist = intersect_tss(bt_tss, genome_path, bt_peaks, tss_dist)
    print('-----DONE!-----')

    print('STEP 4 - Adding info to matrix...')
    df_te = te_matrix(matrix, df_int_rep, tss_dist, df_int_tss_dist)
    print('-----DONE!-----')

    print('STEP 5 - Saving file...')
    matrix_stem = Path(out_dir/matrix_name).stem
    df_te.to_csv(f'{out_dir}/{matrix_stem}_class.csv', header=True)
    print('-----DONE!-----')

    print('STEP 6 - Intersect with KZFP binding sites...')
    bt_kzfp = bt.BedTool(f'{utils_dir}/{genome}_kzfps_combined.bed')
    kzfp_meta = pd.read_csv(f'{utils_dir}/kzfp_info.csv')

    df_int = intersect_kzfp(bt_peaks, bt_kzfp, kzfp_meta)
    print('-----DONE!-----')

    print('STEP 7 - Adding info to matrix...')
    df_final, df_znflist = final_matrix(df_te, df_int)
    print('-----DONE!-----')

    print('STEP 8 - Saving files...')
    df_final.to_csv(f'{out_dir}/{matrix_stem}_class_kzfps.csv', header=True)
    df_znflist.to_csv(f'{out_dir}/{matrix_stem}_PeakZNFList.csv')
    print('-----DONE!-----')
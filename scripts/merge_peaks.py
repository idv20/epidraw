from pathlib import Path, PurePath
import pybedtools

def merge_peaks(bedtool_beds):
    """Concatenates peaks to create a final merged set.

    Args:
        bedtool_beds (list): List of all .bed files in directory.

    Returns:
        bedtool: BedTool object containing merged peaks.
    """
    merged_peaks = bedtool_beds[0].cat(*bedtool_beds[1:]) \
    if len(bedtool_beds) > 1 \
    else bedtool_beds[0]

    merged_peaks_df = merged_peaks.to_dataframe()
    merged_peaks_df['name'] = list(map(lambda x: f'peak_{x}', merged_peaks_df.index))
    merged_peaks_df['score'] = 0
    merged_peaks_df['strand'] = '+'
    labeled_merged_peaks = pybedtools.BedTool.from_dataframe(merged_peaks_df)

    return labeled_merged_peaks

if __name__ == '__main__':

    # Extract from snakemake config
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    bam_dir_name = snakemake.config['bam_dir_name']

    # Set directories
    bam_dir = exp_dir/bam_dir_name
    bed_dir_name = 'final_peaks'
    bed_dir = exp_dir/bed_dir_name
    out_dir = exp_dir/'matrix'

    bedtool_beds = [pybedtools.BedTool(str(bed_dir/f'{file}')) for file in bed_dir.iterdir()]

    print('STEP 1 - Merging peaks...')
    merged_peaks = merge_peaks(bedtool_beds)
    merged_peaks.saveas(f'{out_dir}/{exp}_merged_peaks.bed')
    print('-----DONE!-----')
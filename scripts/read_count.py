import numpy as np
from pathlib import Path, PurePath
import pybedtools
import pysam
import pandas as pd
from collections import deque

def generate_count_matrix(bam_file, sample_name, merged_peaks):
    """Counts number of reads per BED interval.

    Args:
        bam_file (path): Path to BAM file.
        sample_name (string): Stem of BAM file path.
        merged_peaks (bedtool): BedTool object containing merged peaks.

    Returns:
        dataframe: Read counts at each peak.
    """
    score_matrix = deque()
    pysam_bam = pysam.AlignmentFile(str(bam_file), 'rb')

    # Iterate through each peak and extract information
    for interval in merged_peaks:
        interval_info = [interval.name, interval.chrom,
                        interval.start, interval.stop]

        score = pysam_bam.count(
                                contig=interval.chrom,
                                start=interval.start,
                                stop=interval.stop
                                )

        # Create rows with peak info and number of reads
        row = interval_info + [score] #concat list + list
        score_matrix.append(row)

    return pd.DataFrame(np.array(score_matrix), columns=['name', 'chrom', 'start', 'end'] + [sample_name])

if __name__ == '__main__':

    # Extract from snakemake config
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    out_dir = exp_dir/'matrix'
    bam_file = Path(snakemake.input[0])
    sample_name = bam_file.stem
    merged_peaks = pybedtools.BedTool(f'{out_dir}/{exp}_merged_peaks.bed')

    print('STEP 1 - Generating counts...')
    score_matrix = generate_count_matrix(bam_file, sample_name, merged_peaks)
    score_matrix.to_csv(f'{out_dir}/{sample_name}_read_counts.csv')
    print('-----DONE!-----')
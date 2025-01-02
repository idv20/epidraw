import glob
import os
from pathlib import Path
from subprocess import PIPE, Popen
import pandas as pd
import plotly.express as px
import pybedtools
import pysam

class OutputPaths:
    """Defines name of relevant directories."""

    def __init__(self, root):
        """Forms the name of relevant directories.

        Args:
            root (string): Path to parent directory of the subdirectories.
        """
        self.target_dir = root
        self.excluded_dir = self.target_dir + '/excluded'
        self.accepted_dir = self.target_dir + '/accepted'
        self.final_dir = self.target_dir + '/final_peaks'
        self.stats_dir = self.target_dir + '/stats'
        self.samples_rejected_dir = self.target_dir + '/samples_rejected'

def bam_read_count(bamfile):
    """Returns a tuple of the number of mapped and unmapped reads in BAM (chr1-22,X).

    Args:
        bamfile (string): String of the path to BAM.

    Returns:
        tuple: (mapped, unmapped) reads in BAM.
    """
    # Gets counts for chr1-22,X
    chrs = [f'chr{x}' for x in range(1,23)]
    chrs.append('chrX')
    p = Popen(['samtools', 'idxstats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        # Ref seq name (chr nr), Ref seq lenth, nr mapped reads, nr unmapped reads
        rname, rlen, nm, nu = line.rstrip().split()
        # Only count reads from chr1-22,X
        if rname.decode('ascii') in chrs:
            mapped += int(nm)
            unmapped += int(nu)
    return (mapped, unmapped)

def count_reads(bed_file, bam):
    """Counts reads and parameters for peaks [outside blacklist and in domain].

    Args:
        bed_file (bedtool): BedTool object of raw BED file.
        bam (string): String of the path to BAM.

    Returns:
        dataframe: Summary of counts and metrics for each peak.
    """
    sortedbamfile = pysam.AlignmentFile(bam, 'rb')

    bed_results = []
    total_reads_bam = bam_read_count(bam)[0]

    for peak in bed_file:
        total_count = 0
        pos_count = 0
        neg_count = 0
        mapq_count = 0
        nm_count = 0
        pos_neg_ratio = 0
        unique_reads = set()
        window = sortedbamfile.fetch(peak.chrom, peak.start, peak.end)

        for alignment in window:
            # Only does calculations for mapped reads
            if alignment.is_unmapped is False:

                # Counts reads in peak interval
                total_count += 1

                # Adds to positive strand count
                if alignment.is_reverse is False:
                    pos_count += 1
                    unique_reads.add(alignment.reference_start) # is a set - no duplicates

                # Adds to negative strand count
                else:
                    neg_count += 1
                    unique_reads.add(alignment.reference_end)

                # Adds all the mapq scores to divide later by count
                mapq_count += alignment.mapping_quality
                # NM count of mismatches (edit distance from reference)
                nm_count += alignment.get_tag('NM') if alignment.has_tag('NM') else 0

            unique_reads_count = len(unique_reads)

            # Avoid division by zero
            unique_ratio = total_count / max(1, unique_reads_count)

            # Expect equal distribution of reads on both strands
            if neg_count != 0:
                pos_neg_ratio = min(100, float(pos_count) / float(neg_count)) # removes artifacts
            else:
                pos_neg_ratio = 100
            if pos_neg_ratio < 1:
                if pos_neg_ratio == 0:
                    pos_neg_ratio = 100
                else:
                    pos_neg_ratio = min(100, 1 / pos_neg_ratio)

            total_count = max(total_count, 1)

        bed_result = {
            'chr': peak.chrom,
            'start': peak.start,
            'end': peak.end,
            'name': peak.name,
            'score': peak.score,
            'total_reads_peak': total_count,
            'pos_reads': pos_count,
            'neg_reads': neg_count,
            'pos_neg_ratio': pos_neg_ratio,
            'mapq_mean': float(mapq_count) / float(total_count),
            'nm_mean': float(nm_count) / float(total_count),
            'total_reads_bam': total_reads_bam,
            'unique_ratio': unique_ratio,
            'unique_reads': unique_reads_count,}
        bed_results.append(bed_result)

    sortedbamfile.close()

    if bed_results:
        return pd.DataFrame(bed_results, columns=bed_results[0].keys())

def df_filter_to_csv(filt, df, output_file):
    """Helper function which allows easier filtering of dataframe.

    Args:
        filt (bool array): Conditions to subset dataframe by.
        df (dataframe): Dataframe which will be subsetted.
        output_file (string): Path and name of file in which to save filtered dataframe as BED.

    Returns:
        dataframe: Filtered dataframe.
    """
    filtered = df.loc[filt]
    with open(output_file, 'w') as output:
        filtered.to_csv(output, sep='\t', index=False, header=False, columns=['chr', 'start', 'end', 'name', 'score'])

    return filtered


if __name__ == '__main__':

    # Extract from snakemake config
    genome = snakemake.config['genome']
    exp_dir = Path(snakemake.config['exp_dir'])
    bam_dir_name = snakemake.config['bam_dir_name']
    bed_dir_name = snakemake.config['bed_dir_name']
    bed = snakemake.input[0] # path to bed file
    bed_file = pybedtools.BedTool(bed).saveas()

    # Extract filtering parameters from snakemake config
    reset = snakemake.config['reset']
    pos_neg_ratio_threshold = float(snakemake.config['pos_neg_ratio_threshold'])
    score_threshold = float(snakemake.config['score_threshold'])
    unique_ratio_threshold = float(snakemake.config['unique_ratio_threshold'])
    mapq_threshold = float(snakemake.config['mapq_threshold'])
    nm_threshold = float(snakemake.config['nm_threshold'])
    total_reads_threshold = float(snakemake.config['total_reads_threshold'])

    # Load blacklist and domain
    utils_dir = Path.cwd()/'utils'
    blacklist = sorted(utils_dir.glob(f'{genome}*blacklist*.bed'))[0]
    blacklist = pybedtools.BedTool(blacklist)
    domain = sorted(utils_dir.glob(f'{genome}*domain.bed'))[0]
    domain = pybedtools.BedTool(domain)

    # Define the output paths
    output_paths = OutputPaths(root=f'{exp_dir}')
    bam_dir = exp_dir/bam_dir_name
    bed_dir = exp_dir/bed_dir_name
    diag_plots_dir = exp_dir/'diagnostics'

    os.chdir(bed_dir)

    # Set format and path to writing bed files
    bed_stem = Path(bed).stem
    bam = glob.glob(f"{bam_dir}/{bed_stem.replace('_peaks', '')}.bam")[0]
    stats_fname = f'{output_paths.stats_dir}/{bed_stem}_stats.csv'

    excluded_path = f'{output_paths.excluded_dir}/{bed_stem}'
    accepted_path = f'{output_paths.accepted_dir}/{bed_stem}'
    final_path = f'{output_paths.final_dir}/{bed_stem}'
    samples_rejected_path = f'{output_paths.samples_rejected_dir}/{bed_stem}'

    # Remove peaks in blacklist or out of domain (exclude chr outside chr1-22,X)
    init_bed_file = bed_file
    init_bed_file_len = len(init_bed_file)
    print(f'----ALL PEAKS: {init_bed_file_len}')
    bed_file = bed_file.intersect(blacklist, v=True).saveas()
    bed_file = bed_file.intersect(domain, u=True).saveas()
    temp_bed_file_len = len(bed_file)
    print(f'----PEAKS outside BLACKLIST AND in DOMAIN: {temp_bed_file_len}')
    removed = init_bed_file.intersect(bed_file, v=True).saveas(f'{excluded_path}_on_blacklist_out_domain.bed')

    # If peaks in domain and not on blacklist
    if len(bed_file) > 0:

        # If stats file exists and reset=False ==> grabs the existing stats
        if os.path.isfile(stats_fname) and reset is not True:

            # Read pre-existing file with peaks that pass signal quality thresholds
            bed_signal_filtered = pybedtools.BedTool(f'{accepted_path}_processed_signal.bed')

        # If stats file does not exists and/or reset=True ==> it recalculates the stats
        else:
            combined = count_reads(bed_file, bam)
            # Save calculations into stats csv
            output = open(stats_fname, 'w')
            combined.to_csv(output, index=False)
            output.close()

            output_files_and_filters = {
                f'{excluded_path}_bad_ratio.bed': combined['pos_neg_ratio'] > pos_neg_ratio_threshold,
                f'{excluded_path}_bad_unique_ratio.bed': combined['unique_ratio'] > unique_ratio_threshold,
                f'{excluded_path}_low_mapq.bed': combined['mapq_mean'] <= mapq_threshold,
                f'{excluded_path}_many_snps.bed': combined['nm_mean'] >= nm_threshold
                }

            # Create BED of excluded peaks based on above criteria
            for output_file, filt in output_files_and_filters.items():
                df_filter_to_csv(filt=filt,
                                df=combined,
                                output_file=output_file)

            # Creates BED of peaks which fit all signal criteria (MAPQ KEPT EVEN IF EXCLUDED CREATED)
            combined_filtered = df_filter_to_csv(filt=((combined['pos_neg_ratio'] <= pos_neg_ratio_threshold) &
                                                (combined['unique_ratio'] <= unique_ratio_threshold) &
                                                (combined['nm_mean'] < nm_threshold)),
                                                df=combined,
                                                output_file=f'{accepted_path}_processed_signal.bed')

            # Creates separate files for accepted peaks based on MAPQ
            output_files_and_filters = {
                f'{accepted_path}_processed_signal_low_mapq.bed': combined_filtered['mapq_mean'] <= mapq_threshold,
                f'{accepted_path}_processed_signal_high_mapq.bed': combined_filtered['mapq_mean'] > mapq_threshold
                }

            for output_file, filt in output_files_and_filters.items():
                df_filter_to_csv(filt=filt,
                                df=combined_filtered,
                                output_file=output_file)

            # Read the BED of accepted peaks that pass signal quality thresh created above
            bed_signal_filtered = pybedtools.BedTool(f'{accepted_path}_processed_signal.bed')

            # Plot histogram of scores to help adjust the MACS score threshold when repeating if needed
            df_signal_filtered = bed_signal_filtered.to_dataframe()
            fig = px.histogram(df_signal_filtered, x='score',
                                nbins=int(df_signal_filtered['score'].max()+1),
                                title='Counts of signal-filtered peaks based on MACS score')
            fig.update_layout(xaxis = dict(tickmode='linear')) # adds all x-axis ticks
            fig.write_html(f'{diag_plots_dir}/{bed_stem.replace("_peaks","")}_peak_counts_post_qual_filt.html')

        # Remove peaks from original .bed with low MACS score
        bed_score_filtered = bed_file.filter(lambda x: float(x.score) >= score_threshold).saveas()

        # Keep peaks with high MACS score that overlap with quality signal peaks
        bed_final_filtered = bed_score_filtered.intersect(bed_signal_filtered, u=True).saveas()

        # Save final set of peaks that pass BOTH signal AND MACS checks
        bed_final_filtered.saveas(f'{accepted_path}_processed_signal_score.bed')
        bed_final_filtered.saveas(f'{final_path}_processed_signal_score.bed')

        # Create a file for rejected peaks (that are in original bed but not in final)
        if len(bed_final_filtered) > 0:
            excluded_total = bed_file.intersect(bed_final_filtered, v=True).saveas()
            excluded_total.saveas(f'{excluded_path}_excluded_total.bed')

        else:
            print(f'!!! No peaks passed the quality thresholds and MACS score for {bed_stem} !!!')

    else:
        print(f'!!! {bed_stem} contains only peaks on ENCODE blacklist and/or outside of domain !!!')
        with open(f'{samples_rejected_path}/{bed_stem}_only_blacklist_out_domain.txt', 'w') as f:
                f.write('This sample had no peaks in domain and outside of blacklist')

pybedtools.cleanup()
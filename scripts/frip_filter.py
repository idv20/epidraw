from pathlib import Path, PurePath
import deeptools.countReadsPerBin as crpb
import pysam
import pybedtools

def frip_calc(bed_file, bam_dir, frip_dir, acc_dir, rej_dir, alpha):
    """Calculates fraction of reads in peaks score and splits into 'accepted' and 'rejected'.

    Args:
        bed_file (path): Path of BED file.
        bam_dir (path): Path to BAM directory.
        frip_dir (path): Path to FRiP analysis output directory.
        acc_dir (path): Path to directory of samples that pass alpha threshold.
        rej_dir (path): Path to directory of samples that fail alpha threshold.
        alpha (float): Threshold of pass/fail FRiP analysis.
    """
    with open(bed_file, 'r') as read_file:
        # Read the first character
        first_char = read_file.read(1)
        bed_name = bed_file.stem
        sample_name = bed_name.replace('_peaks', '')

        # If variable exists, calculate frip
        if first_char:

            bam_name = bed_name.replace('_peaks', '.bam')
            bam_file = Path(bam_dir/bam_name)

            # Run deeptools to count reads founds in each peak and add them
            cr = crpb.CountReadsPerBin([bam_file],
                                    bedFile=str(bed_file),
                                    numberOfProcessors=8)
            reads_at_peaks = cr.run()
            total_peaks_reads = int(reads_at_peaks.sum(axis=0).item())

            # Find total number of reads mapped to the genome
            bam_pysam = pysam.AlignmentFile(bam_file)
            total_reads = bam_pysam.mapped

            # Calculate FRiP score and store sample and score in separate folders
            frip = total_peaks_reads/total_reads
            print(frip)

            if frip >= alpha:
                frip = str(round(frip, 3))
                status = 'ACCEPTED'
                with open(acc_dir/f'{sample_name}.txt', 'w') as f:
                    f.write('\t'.join([sample_name, frip]))
            else:
                frip = str(round(frip, 3))
                status = 'REJECTED'
                with open(rej_dir/f'{sample_name}.txt', 'w') as f:
                    f.write('\t'.join([sample_name, frip]))

            # Make summary of numbers
            with open(frip_dir/f'{sample_name}_summary.txt', 'w') as f:
                l1 = f'Summary of {sample_name} stats\n'
                l2 = '------------------------------\n'
                l3 = f'nr_peaks = {len(pybedtools.BedTool(bed_file))}\n'
                l4 = f'nr_reads_peaks = {total_peaks_reads}\n'
                l5 = f'nr_reads_total = {total_reads}\n'
                l6 = f'frip_score = {frip}\n'
                l7 = '------------------------------\n'
                l8 = f'{status}'
                f.writelines([l1, l2, l3, l4, l5, l6, l7, l8])

        # If variable doesn't exist, reject sample
        else:
            status = 'REJECTED'

            with open(rej_dir/f'{sample_name}.txt', 'w') as f:
                f.write('No peak calls for sample')

            with open(frip_dir/f'{sample_name}_summary.txt', 'w') as f:
                l1 = f'Summary of {sample_name} stats\n'
                l2 = '------------------------------\n'
                l3 = 'No peak calls for sample\n'
                l4 = '------------------------------\n'
                l5 = f'{status}'
                f.writelines([l1, l2, l3, l4, l5])

if __name__ == '__main__':

    # Extract from snakemake config
    exp_dir = Path(snakemake.config['exp_dir'])
    exp = PurePath(snakemake.config['exp_dir']).name
    bam_dir_name = snakemake.config['bam_dir_name']
    bam_dir = exp_dir/bam_dir_name
    bed_dir_name = snakemake.config['bed_dir_name']
    bed_dir = exp_dir/bed_dir_name
    alpha = float(snakemake.config['frip_alpha'])
    bed_file = Path(snakemake.input[0])

    # Set directories
    frip_dir = exp_dir/'frip'
    acc_dir = frip_dir/'accepted'
    rej_dir = frip_dir/'rejected'

    print('STEP 1 - Calculating FRiP score...')
    frip_calc(bed_file, bam_dir, frip_dir, acc_dir, rej_dir, alpha)
    print('-----DONE!-----')
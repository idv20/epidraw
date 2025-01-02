import os
from pathlib import Path
import pysam
import subprocess as sub
import pandas as pd

def check_reads(bam_file):
    """Detects if any contigs have no reads.

    Args:
        bam_file (string): Name of BAM file.
    """
    try:
        contigs = pysam.AlignmentFile(bam_file, 'rb').get_index_statistics()
        errors = []

        for contig in contigs:
            if contig.total == 0:
                errors.append(contig.contig)

        if len(errors) == 0:
            print('No empty contigs found.')
        else:
            print(f"WARNING: The following contigs have 0 reads: {', '.join(errors)} The BAM file may be corrupted.")
    except ValueError:
        print('No mapping information available in index, or index not recorded.')

def macs2(bam_file, exp_dir, bam_dir_name, bed_dir_name, pe_mode, inp):
    """Calls MACS2 for peak calling.

    Args:
        bam_file (string): Name of BAM file.
        exp_dir (path): Path to the experiment directory.
        bam_dir_name (string): Name of directory containing BAM files.
        bed_dir_name (string): Name of directory in which to save BED files.
        pe_mode (string): Paired-end or single-end mode (MACS2: 'BAMPE' or 'BAM').
        inp (bool): If sample has an input.
    """
    bam_name = Path(bam_file).stem
    bed_dir = Path(exp_dir/bed_dir_name)

    if inp == True:

        samples = pd.read_csv(snakemake.config['samples'])
        sample_sub = samples[samples.name == bam_name]
        input_name = sample_sub.input.item()

        command = f'macs2 callpeak -t {bam_file} -c {exp_dir}/{bam_dir_name}/{input_name}.bam -f {pe_mode} -g hs -n {bam_name} --outdir {bed_dir}'

    else:
        command = f'macs2 callpeak -t {bam_file} -f {pe_mode} -g hs -n {bam_name} --outdir {bed_dir}'

    sub.run(command.split(), text=True)

def bedtools(bam_file, exp_dir, bed_dir_name, genome, slop, out_file):
    """Uses bedtools to extend around summit of peaks.

    Args:
        bam_file (string): Name of BAM file.
        exp_dir (path): Path to the experiment directory.
        bed_dir_name (string): Name of directory in which to save BED files.
        genome (string): 'hg19' or 'hg38' or 'mm10'.
        slop (int): Number of bases to extend on each side of peak summit.
        out_file (string): String of the path and name of output file inferred from wildcards.
    """
    bam_name = Path(bam_file).stem
    bed_dir = Path(exp_dir/bed_dir_name)

    summits = Path(bed_dir)/(bam_name+'_summits.bed')
    genome_file = Path.cwd()/'utils'/(genome+'.chrom.sizes')

    command = f'bedtools slop -i {summits} -g {genome_file} -b {slop} > {out_file}'
    print(command)
    os.system(command)

if __name__ == "__main__":

    # Extract from snakemake config
    genome = snakemake.config['genome']
    pe_mode = snakemake.config['paired-end']
    inp = snakemake.config['input']
    slop = int(snakemake.config['slop'])
    exp_dir = Path(snakemake.config['exp_dir'])
    bam_dir_name = snakemake.config['bam_dir_name']
    bed_dir_name = snakemake.config['bed_dir_name']
    samples = pd.read_csv(snakemake.config['samples'])
    sample_name = Path(snakemake.input[0]).stem

    print('STEP 1 - Checking total and per contigs reads...')
    check_reads(snakemake.input[0])
    print('-----DONE!-----')

    print('STEP 2 - Peak calling with MACS2...')
    # Check if paired-end or single-end
    if pe_mode == True:
        pe_mode = 'BAMPE'

    elif pe_mode == 'Mixed':
        pe_mode = samples[samples['name'] == sample_name]['type_end'].item()
        if pe_mode == 'paired':
            pe_mode = 'BAMPE'
        else:
            pe_mode = 'BAM'

    else:
        pe_mode = 'BAM'

    # Check if input is present for sample
    if inp == 'Mixed':
        input_file = samples[samples['name'] == sample_name]['input'].item()
        if input_file == 'None':
            inp = False
        else:
            inp = True

    macs2(snakemake.input[0], exp_dir, bam_dir_name, bed_dir_name, pe_mode, inp)
    print('-----DONE!-----')

    print('STEP 3 - Extending bases around peak summits...')
    out_file = snakemake.output[0]
    bedtools(snakemake.input[0], exp_dir, bed_dir_name, genome, slop, out_file)
    print('-----DONE!-----')
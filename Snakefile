import pandas as pd
from pathlib import Path, PurePath

configfile: "config.yaml"

samples = pd.read_csv(f"{config['samples']}")
SAMPLES = list(samples['name'])

exp_dir = config['exp_dir']
exp = PurePath(exp_dir).name
fig_dir = Path(exp_dir)/'figures'
bed_dir_name = config['bed_dir_name']
bam_dir_name = config['bam_dir_name']
enrich_thresh_te = config['enrich_thresh_te']

# MAKE ALL NECESSARY DIRECTORIES
bed_dir = Path(f'{exp_dir}/{bed_dir_name}')
frip_dir = Path(f'{exp_dir}/frip')
acc_dir = Path(f'{frip_dir}/accepted')
rej_dir = Path(f'{frip_dir}/rejected')
acc_bed_dir = Path(f'{exp_dir}/accepted')
excl_dir = Path(f'{exp_dir}/excluded')
final_dir = Path(f'{exp_dir}/final_peaks')
stats_dir = Path(f'{exp_dir}/stats')
sample_rej_dir = Path(f'{exp_dir}/samples_rejected')
diag_plots_dir = Path(f'{exp_dir}/diagnostics')
out_dir = Path(f'{exp_dir}/matrix')
fig_dir = Path(f'{exp_dir}/figures')

directories = [bed_dir, frip_dir, acc_dir, rej_dir, acc_bed_dir, excl_dir, final_dir,
                stats_dir, sample_rej_dir, diag_plots_dir, out_dir, fig_dir]
for direc in directories:
    if not direc.is_dir():
        direc.mkdir()

rule all:
    input:
        expand(f"{bed_dir}/{{sample}}_peaks.bed", sample=SAMPLES),
        expand(f"{bed_dir}/{{sample}}_summits.bed", sample=SAMPLES),
        expand(f"{frip_dir}/{{sample}}_summary.txt", sample=SAMPLES),
        expand(f"{final_dir}/{{sample}}_peaks_processed_signal_score.bed",sample=SAMPLES),
        f"{out_dir}/{exp}_merged_peaks.bed",
        expand(f"{out_dir}/{{sample}}_read_counts.csv", sample=SAMPLES),
        f"{out_dir}/{exp}_count_matrix.csv",
        f"{out_dir}/{exp}_score_matrix_norm.csv",
        f"{out_dir}/{exp}_score_matrix_norm_class.csv",
        expand(f"{fig_dir}/{exp}_UMAP_{{sample}}.jpeg", sample=SAMPLES),
        f"{out_dir}/{exp}_graph_leiden.pkl",
        f"{fig_dir}/{exp}_clust_leiden.html",
        f"{fig_dir}/{exp}_te_clustmap_thresh_{enrich_thresh_te}.jpeg"

rule peak_call:
    input:
        f"{exp_dir}/{bam_dir_name}/{{sample}}.bam"
    output:
        f"{exp_dir}/{bed_dir_name}/{{sample}}_peaks.bed",
        f"{exp_dir}/{bed_dir_name}/{{sample}}_summits.bed"
    threads: 4
    script:
        "scripts/peak_call.py"

rule frip_filter:
    input:
        f"{exp_dir}/{bed_dir_name}/{{sample}}_peaks.bed",
        f"{exp_dir}/{bed_dir_name}/{{sample}}_summits.bed"
    output:
        f"{exp_dir}/frip/{{sample}}_summary.txt"
    script:
        "scripts/frip_filter.py"

rule bed_filter:
    input:
        f"{exp_dir}/{bed_dir_name}/{{sample}}_peaks.bed",
        f"{exp_dir}/frip/{{sample}}_summary.txt"
    output:
        f"{exp_dir}/final_peaks/{{sample}}_peaks_processed_signal_score.bed"
    script:
        "scripts/bed_filtering.py"

rule merge_peaks:
    input:
        expand(f"{exp_dir}/final_peaks/{{sample}}_peaks_processed_signal_score.bed",
                sample=SAMPLES)
    output:
        f"{exp_dir}/matrix/{exp}_merged_peaks.bed"
    script:
        "scripts/merge_peaks.py"

rule count_reads:
    input:
        f"{exp_dir}/{bam_dir_name}/{{sample}}.bam",
        f"{out_dir}/{exp}_merged_peaks.bed"
    output:
        f"{out_dir}/{{sample}}_read_counts.csv"
    script:
        "scripts/read_count.py"

rule score_matrix:
    input:
        expand(f"{out_dir}/{{sample}}_read_counts.csv",
                sample=SAMPLES)
    output:
        f"{out_dir}/{exp}_count_matrix.csv",
        f"{out_dir}/{exp}_score_matrix_norm.csv"
    script:
        "scripts/score_matrix_median.py"

rule classify_peaks:
    input:
        f"{out_dir}/{exp}_score_matrix_norm.csv"
    output:
        f"{out_dir}/{exp}_score_matrix_norm_class.csv"
    script:
        "scripts/peak_classifier.py"

rule umap:
    input:
        f"{out_dir}/{exp}_score_matrix_norm_class.csv"
    output:
        expand(f"{fig_dir}/{exp}_UMAP_{{sample}}.jpeg", sample=SAMPLES),
        f"{out_dir}/{exp}_graph_leiden.pkl",
    script:
        "scripts/UMAP_leiden.py"

rule cluster:
    input:
        f"{out_dir}/{exp}_graph_leiden.pkl"
    output:
        f"{fig_dir}/{exp}_clust_leiden.html"
    script:
        "scripts/clustering_leiden.py"

rule enrich:
    input:
        f"{fig_dir}/{exp}_clust_leiden.html",
    output:
        f"{fig_dir}/{exp}_te_clustmap_thresh_{enrich_thresh_te}.jpeg"
    script:
        "scripts/enrichment.py"

#* REMOVE samples_rejected folder if empty after all jobs
samples_rej = sorted(sample_rej_dir.glob('*'))
if (len(samples_rej) == 0) & (sample_rej_dir.is_dir()):
    sample_rej_dir.rmdir()
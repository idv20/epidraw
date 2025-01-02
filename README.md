# EpiDRAW: Epigenetic Dimensionality Reduction Analysis Workflow

## Method overview
The EpiDRAW input consists in genomic alignment files in BAM format.

(A) Peaks are called in each sample analysed. Instead of relying on a binary peak quantification based on an arbitrarily set threshold of significance, EpiDRAW uses a ‘peak union’ method. If a peak was called in one sample at a genomic region, all samples are queried at the same location and the reads are counted.

(B) The resulting read count matrix is normalised and UMAP is applied for dimensionality reduction. The 2D projection of the matrix displays every peak as a data point with a coloured overlay based on the magnitude of the normalised read count in each sample.

(C) The Leiden clustering algorithm is used on the underlying UMAP graph.

(D) Each cluster of peaks can then be used for downstream enrichment analyses of features of interest.


## Installing the Linux environment
You can install the environment on a Linux machine from the epidraw.yml file. To install the environment using Conda use:

```
conda create env --name <env> --file epidraw.yml
```

## EpiDRAW input requirements
### 1. CSV file of sample names
The CSV file contains the names of samples to be analysed. Names need to be provided without the __.bam__ file extension. All BAM files need to be indexed prior to starting the analysis. The CSV file should follow this template:

| __name__ | __input__ |
| :--------: | :---------: |
| _sample_1_ | _sample_1_IN_ |
| _sample_2_ | _sample_2_IN_ |
| _..._ | _..._ |
| _sample_n_ | _sample_n_IN_ |

If the samples do not have inputs, you can leave the __input__ column header in the CSV file without filling any cells below it.

### 2. Config file
The config file specifies the parameters required for the EpiDRAW analysis.

### 3. Snakefile
The Snakefile specifies the order in which the Python scripts should be run as part of the Snakemake EpiDRAW workflow.

## Running the EpiDRAW workflow.
Provided that all 3 requirements mentioned above have been fulfilled, the EpiDRAW workflow can be run fully.

Activate your environment (e.g. _epidraw_) using
```
conda activate epidraw
```
Then make sure you are in the same directory as your Snakemake file is
```
cd path/to/directory/containing/Snakefile
```
Then run the workflow by specifying the number of cores you want to utilise
```
snakemake --cores <number of cores>
```

## More information
You can find more information on how to run a Snakemake workflow by following [this link](https://snakemake.github.io/).# epidraw

# Analyse Breaks

- Table of contents
  - [Workflow](#workflow)
  - [Installation and usage](#installation-and-usage)
  - [Parameters](#parameters)
  - [Outputs](#outputs)

A Nextflow pipeline designed to process INDUCE-Seq reads and provide sample-level information on break counts.
Primary outputs include normalised break counts and assignment of samples to control/treatment groups.

The pipeline assumes data has been through:
- Mapping to reference genome
- Conversion to genomic intervals stored in BED files, including read metadata and MapQ scores
- Processed such that coordinates refer to only the break sites

A reference file (BED format) of canonical cut sites must be provided, with one for AsiSI sites on GRCh38 chromosome 21 included by default.

## Workflow
The pipeline carries out the following steps:
1. **Filtering** of each input BED based on MapQ value (default: 30)
2. Further **filtering** of each BED to remove sites which do not appear in the reference file (default: `chr21_AsiSI_sites.t2t.be`)
3. Per-sample **normalisation** of canonical breaks against total break count
4. **Aggregation** of sample statistics
5. **Analysis** of aggregate data to :
   - Produce QC plots
   - Assign to control/treated/ambiguous groups (Gaussian Mixture model)
   - Statistical testing of assigned groups (Mann-Whitney U test)

## Installation and usage
This pipeline runs via Nextflow. It was written against v25.1 and has not been tested with other versions.

By default it is executed with SLURM, it may be changed to local execution as instructed in `nextflow.config`.

Software dependencies are managed with either Mamba (default) or Conda.

Clone the repository and execute a test run as follows:

``` bash
git clone https://github.com/sholt6/analyse-breaks.git

nextflow run main.nf
```

The default `nextflow.config` included will run an example analysis using files included under the `data/` and `ref/` directories.

## Parameters
Several parameters are exposed, with defaults provided via `nextflow.config`:
```
# Location of a directory containing appropriately formatted BED files:
params.inputDir = "${projectDir}/data/"

# Minimum MapQ value, reads below this quality will be filtered out:
params.mapQThreshold = 30 

# Path to reference BED file containing canonical cut sites:
params.referenceBED = "${projectDir}/ref/chr21_AsiSI_sites.t2t.bed"

# Normalised break sites metric is calculated as:
# $$ totalBreaks / ( canonicalBreaks / normalisationValue ) $$
params.normalisationValue = 1000

# A top-level name to be applied to output files:
params.experimentName = "example"
```

The default settings will run analysis of the example dataset.

Alter parameters with a copy of the config file or via the command line:
``` bash
nextflow run main.nf \
  --inputDir=/path/to/bed_dir/
  --referenceBED=/path/to/reference.bed
```

The assumed filename format for input bed files is `{Sample_Name}.breakends.bed`.


## Outputs
Three output files are provided in the `results` directory:

- `{output_prefix}.summary.txt`: provides top-level statistics for amount of samples assigned to each group and median values for control and treatment categories.
Includes results of a Mann-Whitney U test which evaluates the statistical significance of the difference between the control/treated clusters. 
While not an outright validation of the validity of individual assignments, this value helps with evaluation of the clustering process.

- `{output_prefix}.assignments.tsv`: contains sample-level information breakend counts, as well as assignment to control, treated, or ambiguous categories and probability values for this assignment.
*Note that for small batches of samples, clustering may become less meaningful*.

- `{output_prefix}.qc_plots.png`: histogram, boxplot, and violin plot showing distribution of normalised break count values. This should be considered carefully for QC purposes, as small samples batches or those without a bimodal distribution may cluster poorly.

In addition to these, various processes produce log files for tracking of inputs and the durations of individual functions.
These include:
- `mapQFilterReads`
- `countNormaliseBreaks`
- `plotStats`

Logs are not published as results but may be found in the relevant process directoies.
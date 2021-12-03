# MASV: A **M**is-**A**ssembly detection and **S**tructural **V**ariant calling pipeline.
A **M**is-**A**ssembly detection and **S**tructural **V**ariant calling pipeline.

![logo.png](https://github.com/Dfupa/MASV-pipeline/blob/master/masv-logo-small.png)

## Getting Started
First it is required to pre-install conda with either anaconda or miniconda3. 

**anaconda** - please follow these [instructions](https://docs.anaconda.com/anaconda/install/)

Alternatively, you can install:

**miniconda3** - please follow these [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Pipeline Installation
After installing conda, download and install the pipeline

```Shell
# First clone the MASV repository
$ git clone https://github.com/Dfupa/MASV-pipeline.git
# Change to MASV directory
$ cd MASV_pipeline/
# Set the conda environment with all the necessary dependencies
$ conda env create -f MASV_pipeline.yml
# Activate the conda environment
$ conda activate MASV_pipeline
```

# Input
This pipline uses as an input a config file located in the MASV_pipeline/lib/config/ directory.  The user can find in that directoy a script called MASV_get_config.py which will build a config file with the parameters provided by the user.

More information is provided by using the help parameter of MASV_get_config.py

```Shell
# Always assuming the MASV_pipeline environment is active
$ python3 lib/config/MASV_get_config.py -h
```

# Target rules
The target rules currently available to use are:

- rule **all**: Outputs both [Svim](https://github.com/eldariont/svim) and [Sniffles](https://github.com/fritzsedlazeck/Sniffles) filtered calls, as well alignment QC and evaluation output (only if a truth dataset has been provided to the config file).
- rule **mapping_only**: Outputs only the BAM alignment between the reads and the reference.
- rule **sniffles**: Outputs Sniffles filtered calls.
- rule **svim**: Outputs Svim filtered calls.
- rule **eval_sniffles**: Outputs evaluation output for Sniffles only.
- rule **eval_svim**: Outputs evaluation output for Sniffles only.
- rule **sanity_check**: Outputs alignment QC based on [mosdepth](https://github.com/brentp/mosdepth) and [Nanoplot.](https://github.com/wdecoster/NanoPlot)

# How to run snakemake
In order to run the MASV pipeline just use the following command

```Shell
# Always assuming the MASV_pipeline environment is active
$ snakemake -s MASV_pipeline.smk -r all -j 16 -n
```
The -n option is for a dry run, avoding the execution of any job yet. -j is the number of threads provided to the pipeline.  In this specific example, we selected the rule all using the parameter -r. 

More information regarding Snakemake and its commands can be found through Snakemake [documentation](https://snakemake.readthedocs.io/en/stable/index.html).

## Important considerations

- Right now the pipeline is able to handle .fastq and .fastq.gz. Please keep this in mind when providing your own reads.
- Right now the pipeline uses either [Minimap2](https://github.com/lh3/minimap2) or [Ngmlr](https://github.com/philres/ngmlr)n aligner, but not both.
- The evaluation script requires one specific SV feature ( default is DEL) and a truth dataset for tha feature types. The user must manually edit the bedtools_eval.smk file in lib/rules/ and account for the feature change. 
- The pipeline was developed and test in Unix like environments. As such, functionality is not assured when used from Windows or Mac.

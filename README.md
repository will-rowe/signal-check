# signal check

[![Build Status](https://travis-ci.org/will-rowe/signal-check.svg?branch=master)](https://travis-ci.org/will-rowe/signal-check)

> WIP

This is a pipeline and analysis notebook for checking the usefulness of signal-data (for viral metagenomics), helping decide if/where to keep it.

## running the analysis

* create and activate the conda environment for the notebook:
  
```
conda env create -f environments/notebook-analysis.yaml 
conda activate notebook-analysis
```

* open the notebook:

```
jupyter notebook analysis.ipynb
```

## standalone running of the pipeline

If you have nextflow and conda installed, or are using the *notebook-analysis* environment from above, you just need:

```
nextflow run main.nf --reads path/to/<reads> --fast5 full/path/to/<fast5 directory> --output <output directory> --assembler miniasm -profile conda
```

## issues

* full path is needed to the fast5 directory - nanopolish and nextflow aren't playing nicely at the moment
* travis does not appear to support the conda environment feature of nextflow, so I've created a single environment for CI testing (environments/travis-testing.yaml) 
* to update an existing notebook-analysis environment: `conda env update -f environments/notebook-analysis.yaml --prune`
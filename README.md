# signal check

[![Build Status](https://travis-ci.org/will-rowe/signal-check.svg?branch=master)](https://travis-ci.org/will-rowe/signal-check)

> WIP

This is a pipeline and analysis notebook for checking the usefulness of signal-data (for viral metagenomics), helping decide if/where to keep it.

## workflow

* demux reads and trim adapters
    * qcat
* assemble
    * miniasm / redbean
* subsample reads
    * pomoxis
* correct assembly
    * racon
* polish assembly
    * medaka
    * nanopolish
    * medaka + nanopolish
* evaluate assemblies
    * quast
    * fastani
    * nucdiff

## running the analysis

* create and activate the conda environment for the notebook:
  
```
conda env create -f pipeline/environments/notebook-analysis.yaml 
conda activate notebook-analysis
```

* open the notebook:

```
jupyter notebook analysis.ipynb
```

## standalone running of the pipeline

If you have nextflow and conda installed, you just need:

```
nextflow run pipeline/long-read-assembly-pipeline.nf --inputDir <full/path/to/directory> --barcodes 09,10,11 --output <output directory> -profile conda --cpus 6 --mem 12GB
```

> to run using Docker instead, swap the `-profile` over to docker

## todo

* add in help message and full param descript
* get more info from the qcat process (using the parsing script)
* add in pre-run checks for reads etc.
* output visualisation of MSA / pileup
* add pycoqc

## issues

* full path is needed to the input directory - nanopolish and nextflow aren't playing nicely at the moment when it comes to collecting the fast5
* travis does not appear to support the conda environment feature of nextflow, so I've created a single environment for CI testing (environments/travis-testing.yaml) 
* to update an existing notebook-analysis environment: `conda env update -f environments/notebook-analysis.yaml --prune`
* `redbean` is currently only supported for linux platforms - I'll hopefully get around to making a conda recipe for OSX

## pipeline dag

![dag](long-read-assembly-pipeline/flowchart.png)
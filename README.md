<div align="center">
    <h3>SIGNAL CHECK</h3>
    <hr>
    <a href="https://travis-ci.org/will-rowe/signal-check"><img src="https://travis-ci.org/will-rowe/signal-check.svg?branch=master" alt="travis"></a>
    <a href="https://www.nextflow.io"><img src="https://img.shields.io/badge/nextflow-%E2%89%A519.07.0-brightgreen.svg" alt="Nextflow"></a>
    <a href=""><img src="https://img.shields.io/badge/status-WIP-orange" alt="Status"></a>
</div>

***

## Overview

This is a set of pipelines and workbooks for checking the usefulness of signal-data (for viral metagenomics); helping to decide if/where to keep it.

### pipelines

There are two nextflow pipelines for long read genome assembly:

* [long-read-assembly-dn](pipelines/long-read-assembly-dn.nf) for `de novo`
* [long-read-assembly-rg](pipelines/long-read-assembly-rg.nf) for `reference guided`

The basic steps of the pipelines are:

* demux and trim basecalled reads
* assemble
* correct the assemblies
* polish, either:
  * without signal
  * with signal
  * with signal first, then without
  * without signal first, then with
* basic assessment of the assemblies

There are a few choices of software for each step (e.g. miniasm or redbean for de-novo assembly)

### workbooks

Each workbook will take you from `data download` -> `de novo genome assembly` -> `assessment` -> `reference guided assembly` -> `assessment`.

There are two works books:

* [analysis-r941_min_fast](analysis-r941_min_fast.ipynb) for analysing data produced from GUPPY `fast` basecalling
* [analysis-r941_min_high](analysis-r941_min_high.ipynb) for analysing data produced from GUPPY `high accuracy` basecalling

> note: only reads generated using the fast model are available for download, you will have to basecall using HAC yourself.

### data

We use data from the latest [artic data release](http://artic.network/protocol_validation_2019.html). In particular, the Ebola virus (EBOV) minion run that amplicon sequenced 3 strains of the virus (Mayinga, Kikwit, Makona) using the rapid PCR kit.


## Running the analysis

* create and activate the conda environment for the notebook:
  
```
conda env create -f pipelines/environments/notebook-analysis.yaml 
conda activate notebook-analysis
```

* open the notebook:

```
jupyter notebook analysis-r941_min_fast.ipynb
```

## Standalone running of the pipeline

If you have nextflow and conda installed, you just need:

```
nextflow run pipelines/long-read-assembly-pipeline-dn.nf --inputDir <full/path/to/directory> --barcodes 09,10,11 --output <output directory> -profile conda --cpus 6 --mem 12GB
```

> to run using Docker instead, swap the `-profile` over to docker

## Pipeline dags

### de novo long read assembly pipeline

![dag](pipelines/long-read-assembly-dn.png)

## Todo

* add in help message and full param descript
* get more info from the qcat process (using the parsing script)
* add in pre-run checks for reads etc.
* output visualisation of MSA / pileup
* add pycoqc
* Redbean assemblies aren't great, I need to try parameterising this better

## Notes

* medaka renames the contigs to include range data, which then breaks Nanopolish - so contigs are renamed sequentially after medaka

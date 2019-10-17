#!/usr/bin/env bash
pipelineDir=${PWD}"/pipelines"
fastqDir=${pipelineDir}"/data/ebov-test-data/fastq_pass"
fast5Dir=${pipelineDir}"/data/ebov-test-data/fast5_pass"
refGenome=${pipelineDir}"/data/ebov-reference-genomes/NC_002549.fasta"
resultsDir=${pipelineDir}"/test-results"
resultCheck1=${resultsDir}"/de-novo-assembly/testing-barcode-09.dn-assembly.raw.fasta"
resultCheck2=${resultsDir}"/de-novo-assembly/testing-barcode-09.dn-assembly.racon.medaka.fasta"
resultCheck3=${resultsDir}"/de-novo-assembly/testing-barcode-09.dn-assembly.racon.nanopolish.fasta"
resultCheck4=${resultsDir}"/reference-guided-assembly/testing-barcode-09.rg-assembly.racon.fasta"

cd ${pipelineDir}

pipeline="nextflow run long-read-assembly.nf \
    -profile docker \
    -with-dag flow.png \
    --cpus 2 \
    --mem 4GB \
    --fastqDir ${fastqDir} \
    --fast5Dir ${fast5Dir} \
    --barcodes 09, \
    --output ${resultsDir} \
    --subSamplingDepth 1 \
    --label testing \
    --refGenome ${refGenome} \
    -resume"

echo "starting assembly pipeline..."
echo $pipeline
echo "-------------------------------------------------------"
eval $pipeline
echo "-------------------------------------------------------"

echo "checking output..."
if test -f "$resultCheck1"; then
    echo "check for initial dn-assembly passed"
else
    echo "check for initial dn-assembly failed"
    exit 1
fi
if test -f "$resultCheck2"; then
    echo "check for medaka polished dn-assembly passed"
else
    echo "check for medaka polished dn-assembly failed"
    exit 1
fi
if test -f "$resultCheck3"; then
    echo "check for nanopolish polished dn-assembly passed"
else
    echo "check for nanopolish polished dn-assembly failed"
    exit 1
fi
if test -f "$resultCheck4"; then
    echo "check for initial rg-assembly passed"
else
    echo "check for initial rg-assembly failed"
    exit 1
fi

echo "cleaning up..."
rm -rf work
rm -rf ${resultsDir}
rm -rf .nextflo*

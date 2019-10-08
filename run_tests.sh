#!/usr/bin/env bash
pipelineDir=${PWD}"/pipelines"
fastqDir=${pipelineDir}"/data/test-ebov-subset/fastq_pass"
fast5Dir=${pipelineDir}"/data/test-ebov-subset/fast5_pass"
resultsDir=${pipelineDir}"/test-results"
resultCheck1=${resultsDir}"/barcode-09.assembly-unpolished.fasta"
resultCheck2=${resultsDir}"/barcode-09.assembly-corrected.medaka-polished.fasta"
resultCheck3=${resultsDir}"/barcode-09.sub-sampled.assembly-corrected.nanopolish-polished.fasta"

cd ${pipelineDir}

testCmd="nextflow run long-read-assembly-dn.nf \
    -profile docker \
    -with-dag flowchart.png \
    --cpus 2 \
    --mem 4GB \
    --fastqDir ${fastqDir} \
    --fast5Dir ${fast5Dir} \
    --barcodes 09, \
    --output ${resultsDir} \
    --subSamplingDepth 1 \
    -resume"

echo "starting nextflow pipeline..."
echo $testCmd
echo "-------------------------------------------------------"
eval $testCmd
echo "-------------------------------------------------------"

echo "checking output..."
if test -f "$resultCheck1"; then
    echo "check for initial assembly passed"
else
    echo "check for initial assembly failed"
    exit 1
fi
if test -f "$resultCheck2"; then
    echo "check for medaka polished assembly passed"
else
    echo "check for medaka polished assembly failed"
    exit 1
fi
if test -f "$resultCheck3"; then
    echo "check for nanopolish polished assembly passed"
else
    echo "check for nanopolish polished assembly failed"
    exit 1
fi

echo "cleaning up..."
rm -rf work
rm -rf ${resultsDir}
rm -rf .nextflo*

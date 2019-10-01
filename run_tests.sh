#!/usr/bin/env bash
pipelineDir=${PWD}"/pipeline"
inputData=${pipelineDir}"/data/ebov-subset"
resultsDir=${pipelineDir}"/test-results"
resultCheck1=${resultsDir}"/barcode-09.assembly-unpolished.fasta"
resultCheck2=${resultsDir}"/barcode-09.assembly-corrected.medaka-polished.fasta"


testCmd="nextflow run long-read-assembly.nf \
    -profile docker \
    --cpus 2 \
    --mem 4GB \
    --inputDir ${inputData} \
    --barcodes 09, \
    --output ${resultsDir} \
    --subSamplingDepth 1 \
    -resume"

cd ${pipelineDir}
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
    exit -1
fi
if test -f "$resultCheck2"; then
    echo "check for medaka polished assembly passed"
else
    echo "check for medaka polished assembly failed"
    exit -1
fi

echo "cleaning up..."
rm -rf work/
rm -rf .nextflow.log*
rm -rf ${resultsDir}

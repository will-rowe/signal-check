#!/usr/bin/env bash
resultsDir=${PWD}"/test-results"
resultCheck1=${resultsDir}"/barcode-09.assembly-unpolished.fasta"

testCmd="nextflow run long-read-assembly.nf -profile conda --cpus 2 --mem 4GB --inputDir $PWD/data/ebov-subset --barcodes 09, --output ${resultsDir} --subSamplingDepth 1 -resume"

echo "Starting nextflow..."
echo $testCmd
echo "-------------------------------------------------------"
eval $testCmd
echo "-------------------------------------------------------"

if test -f "$resultCheck1"; then
    echo "check for initial assembly passed"
else
    echo "check for initial assembly failed"
    exit -1
fi
echo "cleaning up..."
eval "make"
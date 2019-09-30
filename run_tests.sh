#!/usr/bin/env bash
resultsDir=${PWD}"./test-results"
resultCheck1=${resultsDir}"/barcode-09.assembly-unpolished.fasta"

cmd="nextflow run main.nf -profile conda --cpus 2 --mem 4GB --inputDir $PWD/data/ebov-subset --barcodes 09, --output ${resultsDir} --assembler miniasm --subSamplingDepth 1"
echo "Starting nextflow..."
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
echo "-------------------------------------------------------"

if test -f "$resultCheck1"; then
    echo "check for initial assembly passed"
else
    echo "check for initial assembly failed"
    exit -1
fi
echo "cleaning up..."
eval "make"
#!/usr/bin/env bash
resultCheck="results/assembly-unpolished.fasta"
cmd="nextflow run main.nf --cpus 2 --mem 4GB --reads ./data/ecoli_small.fasta --output results -profile conda --assembler miniasm"
echo "Starting nextflow..."
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
echo "-------------------------------------------------------"

if test -f "$resultCheck"; then
    echo "check for initial assembly passed"
else
    echo "check for initial assembly failed"
    exit -1
fi
echo "cleaning up..."
eval "make"
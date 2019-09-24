#!/usr/bin/env bash
cmd="nextflow run main.nf --cpus 2 --mem 4GB --reads ./data/ecoli_small.fasta --output results -profile conda --assembler miniasm"
echo "Starting nextflow..."
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
echo "-------------------------------------------------------"
echo "cleaning up..."
eval "make"
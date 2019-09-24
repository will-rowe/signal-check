#!/usr/bin/env bash
git clone https://github.com/will-rowe/signal-check && cd signal-check

# run tests
run_name="Test run: "$(date +%s)
cmd="nextflow run main.nf --cpus 2 --mem 4GB --reads ./data/ecoli_small.fasta --output results -profile conda --assembler miniasm"
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
#!/usr/bin/env python3

"""
sample: the file of demuxed reads (e.g. barcode-09.fastq)
refMultiFasta: the multifasta of reference genomes
outfile: the name to give the reference fasta

this script will find the reference genome in a multifasta for the provided sample
"""

import argparse
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser()
parser.add_argument("sample")
parser.add_argument("refMultiFasta")
parser.add_argument("outfile")
args = parser.parse_args()

found = False
seqIter = SeqIO.parse(args.refMultiFasta, 'fasta')
for seq in seqIter:
    barcode = seq.id.split("---")[0]
    if barcode in args.sample:
        found = True
        SeqIO.write(seq, args.outfile, "fasta")
        break
if found == False:
    sys.exit("error: could not find reference in multifasta")
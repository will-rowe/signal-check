#!/usr/bin/env python3

"""
inputDir: the directory containing the output from qcat
barcodes: the list of barcodes to keep
label: prepends a label to the output filenames

this script will remove all barcodes we don't want from the qcat output
it will also remain those we do want
"""

import argparse
from os import listdir, rename, remove
from os.path import isfile, join


parser = argparse.ArgumentParser()
parser.add_argument("inputDir")
parser.add_argument("barcodes")
parser.add_argument("label")
args = parser.parse_args()

for fname in listdir(args.inputDir):
    if (isfile(join(args.inputDir, fname))):
        base = fname.split(".")[0]
        if "barcode" in base:
            bc = base.replace("barcode", "")
            if bc in args.barcodes:
                # for now just rename to barcode
                # TODO: use a lookup table to add in a meaningful name here
                if (args.label != ""):
                    rename(args.inputDir + "/" + fname, args.inputDir + "/" + args.label + "-barcode-" + bc + ".fastq")
                else:
                    rename(args.inputDir + "/" + fname, args.inputDir + "/" + "barcode-" + bc + ".fastq")
            else:
                remove(args.inputDir + "/" + fname)
        else:
            remove(args.inputDir + "/" + fname)


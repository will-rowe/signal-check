#!/usr/bin/env python3

"""
Plot the read depth of across an alignment

depthFile: the output of `samtools depth -a sample.bam > sample.depth`
outPrefix: the prefix for the output file
smoothing: subsample the depth values

Note: there is no checking to see if the right file type was provided
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="ticks")

# set up the CLI
parser = argparse.ArgumentParser()
parser.add_argument("depthFile")
parser.add_argument("outPrefix")
parser.add_argument("--smoothing", action="store_true", help="subsample the depth values")
args = parser.parse_args()

# customise some things
outSuffix = "-depthPlot.png"
titlePrefix = "Depth plot: "
xLab = "genomic position"
yLab = "number of reads"

# collect the output from samtools depth
df = pd.read_csv(args.depthFile, sep='\t', header=None)
df.columns = ["refName", "position", "depth"]

# check all ref positions were depth checked
startPos = df["position"][0]
endPos = df["position"][df.index[-1]]
assert len(df.index) == ((endPos - startPos) + 1), "depth needs to be reported at all positions (use `samtools depth -a in.bam`)"

# subsample the depth readings if requested
if (args.smoothing == True):
    #SAMPLE_FREQ = int(0.1 * len(df.index))
    SAMPLE_FREQ = 10
    df = df.iloc[::SAMPLE_FREQ, :]

# check that there aren't too many positions in the depth data
if len(df.index) > 20000:
    sys.stderr.write("\nerror: too many data points to plot, try using the `--smoothing` option\n")
    sys.exit(1)

# set up the plot
plt.clf()
plt.style.use("seaborn")
title = titlePrefix + df["refName"][0] + ":" + str(startPos) + "-"+  str(endPos)
plt.vlines(df["position"], 0, df["depth"], colors=['blue'])

#p1=sns.lineplot(data=df, x="positions", y="depth", color="r")
#p1=sns.lineplot(data=df, x="positions", y="d2", color="b")

# add some bits
plt.title(title, fontsize=18)
plt.gca().set_xlabel(xLab)
plt.gca().set_ylabel(yLab)
plt.tight_layout()

# write the plot
outfile = args.outPrefix + outSuffix
plt.savefig(outfile, dpi=500)

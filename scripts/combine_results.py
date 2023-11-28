#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--pools", required=True, nargs="+", type=str, help="")
parser.add_argument("--counts", required=True, nargs="+", type=str, help="")
parser.add_argument("--barcodes", required=True, nargs="+", type=str, help="")
parser.add_argument("--bam", required=True, nargs="+", type=str, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

import pandas as pd

# Construct the file_directories.tsv file that points to the files we will use for WG1.
file_directories = pd.DataFrame({
    'Pool': args.pools,
    'Counts': args.counts,
    'Barcodes': args.barcodes,
    'Bam': args.bam,
})
print("Saving file directories for WG1.")
file_directories.to_csv(os.path.join(args.out, "wg0_file_directories.tsv"), sep="\t", index=False)

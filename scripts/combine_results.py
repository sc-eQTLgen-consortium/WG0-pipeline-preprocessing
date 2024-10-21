#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--samplesheet", required=True, type=str, help="")
parser.add_argument("--cellbender_selection", required=False, type=str, default=None, help="")
parser.add_argument("--basedir", required=True, type=str, help="")
parser.add_argument("--outfile", required=True, type=str, help="")
args = parser.parse_args()

dir = os.path.dirname(args.outfile)
if dir != "":
    os.makedirs(dir, exist_ok=True)

print("Options in effect:")
arguments = {}
for arg in vars(args):
    print("  --{} {}".format(arg, getattr(args, arg)))
    arguments[arg] = getattr(args, arg)
print("")

import pandas as pd

print("Loading the sample sheet ...")
sample_df = pd.read_csv(args.samplesheet, sep="\t", dtype=str)
SAMPLES = set(sample_df["Sample"].unique())
print("  Loaded dataframe: {} with shape: {}".format(os.path.basename(args.samplesheet), sample_df.shape))

print("Constructing output dataframe ...")
df = sample_df[["Sample"]].rename(columns={"Sample": "Pool"})
del sample_df

if args.cellbender_selection is None:
    print("\tUsing 'Counts' from: CellRanger")
    df["Counts"] = args.basedir + "CellRanger/" + df["Pool"] + "/outs/filtered_feature_bc_matrix.h5"

    print("\tUsing 'Barcodes' from: CellRanger")
    df["Barcodes"] = args.basedir + "CellRanger/" + df["Pool"] + "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
else:
    print("Loading the CellBender selection")
    cb_selection_df = pd.read_csv(args.cellbender_selection, sep="\t", header=0, index_col=None)
    print("  Loaded dataframe: {} with shape: {}".format(os.path.basename(args.cellbender_selection), cb_selection_df.shape))

    # Validate CellBender selection file.
    passed_select_df = cb_selection_df.loc[(cb_selection_df["FINISHED"]) & (cb_selection_df["PASSED"]), ["Sample", "Run"]].copy()
    del cb_selection_df
    if passed_select_df.shape[0] < len(SAMPLES) or len(SAMPLES.difference(set(passed_select_df["Sample"].values))) != 0:
        print("Error, --cellbender_selection is invalid.")
        exit()

    # Merge the matrices together.
    passed_select_df.rename(columns={"Sample": "Pool"}, inplace=True)
    df = df.merge(passed_select_df, on="Pool", how="inner")
    del passed_select_df
    if df.shape[0] < len(SAMPLES) or len(SAMPLES.difference(set(df["Pool"].values))) != 0:
        print("Error, --cellbender_selection is invalid.")
        exit()

    df["PoolRun"] = df["Pool"] + "Run" + df["Run"].astype(str)

    print("\tUsing 'Counts' from: CellBender")
    df["Counts"] = args.basedir +  "CellBender/" + df["PoolRun"] + "/cellbender_feature_bc_matrix_filtered.h5"

    print("\tUsing 'Barcodes' from: CellBender")
    df["Barcodes"] = args.basedir +  "CellBender/" + df["PoolRun"] + "/cellbender_feature_bc_matrix_cell_barcodes.csv"

# Add CellRanger BAM file.
print("\tUsing 'Bam' from: CellRanger")
df["Bam"] = args.basedir +  "CellRanger/" + df["Pool"] + "/outs/possorted_genome_bam.bam"
df = df[["Pool", "Counts", "Barcodes", "Bam"]]

# Check if all files exist.
for _, row in df.iterrows():
    for colname, fpath in row.iteritems():
        if colname == "Pool":
            continue
        if not os.path.exists(fpath):
            print("Error, output file '{}' does not exist.".format(fpath))
            exit()

print("Saving output ...")
df.to_csv(args.outfile, sep="\t", index=False)

print("Done.")

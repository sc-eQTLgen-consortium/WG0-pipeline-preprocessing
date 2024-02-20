#!/usr/bin/env python
# Adapted from https://github.com/broadinstitute/CellBender/blob/master/cellbender/remove_background/report.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--cellranger_dir", required=True, type=str, help="")
parser.add_argument("--cellranger_infolders", required=True, nargs="+", type=str, help="")
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

from scipy.cluster.hierarchy import ward, dendrogram, leaves_list
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


def metrics_heatmap(df, annot_df, outfile, title="", scale=0.75):
    data = df.to_numpy()
    colnames = df.columns.tolist()
    rownames = df.index.tolist()

    fig, ax = plt.subplots(figsize=(df.shape[1] * scale, df.shape[0] * scale))
    im = ax.imshow(data, cmap='bwr')
    fig.colorbar(im, ax=ax)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(colnames)), labels=colnames)
    ax.set_yticks(np.arange(len(rownames)), labels=rownames)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(rownames)):
        for j in range(len(colnames)):
            text = ax.text(j, i, annot_df.iloc[i, j],
                           ha="center", va="center", color="black")

    ax.set_title(title)
    fig.tight_layout()
    plt.savefig(outfile, bbox_inches="tight")


################################################################################

infolders = list(args.cellranger_infolders)
infolders.sort()

# Loading metrics.
metrics = []
for infolder in infolders:
    metrics_summary_path = os.path.join(args.cellranger_dir, infolder, "outs", "metrics_summary.csv")
    metrics_summary = pd.read_csv(metrics_summary_path, sep=',', header=0, index_col=None)
    metrics_summary.index = [infolder]
    metrics.append(metrics_summary)

# Removing unwanted characters.
metrics_df = pd.concat(metrics, axis=0)
for column in metrics_df.columns:
    metrics_df[column] = [float(re.sub("[^0-9^.]", "", str(value))) for value in metrics_df[column]]

# Save output.
metrics_df.to_csv(os.path.join(args.cellranger_dir, 'metrics_summary.tsv.gz'), sep="\t", header=True, index=True, compression="gzip")

# Hierarchical clustering.
Z = ward(pdist(metrics_df))
order = leaves_list(Z)
metrics_df = metrics_df.iloc[order, :]

# Z-transform for color of the heatmap.
z_metrics_df = (metrics_df - metrics_df.mean(axis=0)) / metrics_df.std(axis=0)

# Summarise labels for annotation of the heatmap.
for column in metrics_df.columns:
    if metrics_df[column].mean() > 1e6:
        metrics_df[column] = [str(round((value / 1e6), 1)) + "M" for value in metrics_df[column]]
    elif metrics_df[column].mean() > 1e3:
        metrics_df[column] = [str(round((value / 1e3), 1)) + "K" for value in metrics_df[column]]
    else:
        metrics_df[column] = metrics_df[column].round(1)

# Plot.
metrics_heatmap(df=z_metrics_df,
                annot_df=metrics_df,
                title="CellRanger Metrics",
                outfile=os.path.join(args.out, 'CellRanger_metrics.png'))


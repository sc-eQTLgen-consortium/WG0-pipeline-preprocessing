#!/usr/bin/env python3

"""
File:         visualise_cellbender_results.py
Created:      2023/07/07
Last Changed: 2023/09/13
Author:       M.Vochteloo
"""

# Standard imports.
from __future__ import print_function
from datetime import datetime
import argparse
import pathlib
import math
import os
import re

# Third party imports.
import numpy as np
import pandas as pd
import h5py
import torch
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Local application imports.

# Metadata
__program__ = "Visualise CellBender Results"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

parser = argparse.ArgumentParser(prog=__program__,
                                 description=__description__)

# Add other arguments.
parser.add_argument("--pools", nargs="*", type=str, required=True, help="")
parser.add_argument("--raw_h5", nargs="*", type=str, required=True, help="")
parser.add_argument("--posterior_h5", nargs="*", type=str, required=True, help="")
parser.add_argument("--metrics", nargs="*", type=str, required=True, help="")
parser.add_argument("--log", nargs="*", type=str, required=True, help="")
parser.add_argument("--output", type=str, required=True, help="")
args = parser.parse_args()


def parse_log_file(logfile):
    data = {}

    start_datetime = None
    end_datetime = None
    completed = False
    with open(logfile, 'r') as f:
        for line in f:
            if line.startswith("cellbender remove-background"):
                for part in line.split("--")[1:]:
                    if "=" in part:
                        key, value = part.split("=")
                        try:
                            data["parameter:" + key] = float(value.strip("\n"))
                        except ValueError:
                            pass
                    else:
                        data[part] = 1
            elif "Completed remove-background." in line:
                data["passed"] = True
            elif re.search("(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})", line):
                tmp_datetime = datetime.strptime(re.search("(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})", line).group(0), "%Y-%m-%d %H:%M:%S")
                if start_datetime is None:
                    start_datetime = tmp_datetime
                elif completed:
                    end_datetime = tmp_datetime
            elif re.search("Features in dataset: ([0-9]+) Gene Expression", line):
                data["features"] = int(re.search("Features in dataset: ([0-9]+) Gene Expression", line).group(1))
            elif re.search("([0-9]+) features have nonzero counts.", line):
                data["nonzero_genes"] = int(re.search("([0-9]+) features have nonzero counts.", line).group(1))
            elif re.search("Prior on counts for cells is ([0-9]+)", line):
                data["prior_counts_empty"] = int(re.search("Prior on counts for cells is ([0-9]+)", line).group(1))
            elif re.search("Prior on counts for empty droplets is ([0-9]+)", line):
                data["prior_counts_cell"] = int(re.search("Prior on counts for empty droplets is ([0-9]+)", line).group(1))
            elif re.search("Excluding ([0-9]+) features that are estimated to have <= 0.1 background counts in cells.", line):
                data["excluded_features"] = int(re.search("Excluding ([0-9]+) features that are estimated to have <= 0.1 background counts in cells.", line).group(1))
            elif re.search("Including ([0-9]+) features in the analysis.", line):
                data["included_features"] = int(re.search("Including ([0-9]+) features in the analysis.", line).group(1))
            elif re.search("Excluding barcodes with counts below ([0-9]+)", line):
                data["barcodes_threshold"] = int(re.search("Excluding barcodes with counts below ([0-9]+)", line).group(1))
            elif re.search("Using ([0-9]+) probable cell barcodes, plus an additional ([0-9]+) barcodes, and ([0-9]+) empty droplets.", line):
                match = re.search("Using ([0-9]+) probable cell barcodes, plus an additional ([0-9]+) barcodes, and ([0-9]+) empty droplets.", line)
                data["probable_cell_barcodes"] = int(match.group(1))
                data["additional_barcodes"] = int(match.group(2))
                data["empty_droplets"] = int(match.group(3))
            elif re.search("Largest surely-empty droplet has ([0-9]+) UMI counts.", line):
                data["max_empty_droplet_count"] = int(re.search("Largest surely-empty droplet has ([0-9]+) UMI counts.", line).group(1))
            elif "Learning failed.  Retrying with learning-rate " in line:
                data = {}
                start_datetime = end_datetime
                end_datetime = None
            elif "Completed remove-background." in line:
                completed = True
            else:
                pass
    f.close()

    if start_datetime is not None and end_datetime is not None:
        data["time"] = (end_datetime - start_datetime).total_seconds() / 60.0

    return data


def barplot(df, panels=None, x="x", y="y", panel_column=None,
            palette=None, order=None, title="", filename="plot"):
    if panels is None:
        panels = df[panel_column].unique()
    if order is None:
        order = df[x].unique()
        order.sort()

    nplots = len(panels)
    ncols = math.ceil(np.sqrt(nplots))
    nrows = math.ceil(nplots / ncols)

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex='all',
                             sharey='none',
                             figsize=(6 * ncols, 6 * nrows))
    sns.set(color_codes=True)

    row_index = 0
    col_index = 0
    for i in range(ncols * nrows):
        if nrows == 1 and ncols == 1:
            ax = axes
        elif nrows == 1 and ncols > 1:
            ax = axes[col_index]
        elif nrows > 1 and ncols == 1:
            ax = axes[row_index]
        else:
            ax = axes[row_index, col_index]

        if i < nplots:
            data = df.loc[df[panel_column] == panels[i], :]
            if data.shape[0] == 0:
                continue

            sns.barplot(
                data=data,
                x=x,
                y=y,
                color="black",
                palette=palette,
                order=order,
                ax=ax
            )

            ax.set_xticklabels(ax.get_xmajorticklabels(), rotation=90)
            ax.set_ylim(int(data[y].min() - (data[y].max() * 0.05)), ax.get_ylim()[1])

            ax.set_xlabel("",
                          fontsize=10,
                          fontweight='bold')
            ax.set_ylabel("",
                          fontsize=10,
                          fontweight='bold')
            ax.set_title(panels[i],
                         fontsize=14,
                         fontweight='bold')

        else:
            ax.set_xticklabels(ax.get_xmajorticklabels(), rotation=90)

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    fig.suptitle(title,
                 fontsize=40,
                 fontweight='bold')

    outpath = os.path.join(args.output, filename)
    fig.savefig(outpath)
    print("\tSaved figure: {}".format(os.path.basename(outpath)))
    plt.close()


def plot_per_sample(df, plottype, samples=None, subtitles=None, x="x", y="y", y2=None,
                    alpha=1.0, sample_column="sample", hue=None, color="black", palette=None,
                    xlabel="", ylabel="", ylabel2="", title="", filename="plot"):
    if samples is None:
        samples = df[sample_column].values.tolist()
        samples.sort()
    if subtitles is not None and len(subtitles) != len(samples):
        print("Error, subtitles are not the same length as the samples")
        return

    nplots = len(samples)
    ncols = math.ceil(np.sqrt(nplots))
    nrows = math.ceil(nplots / ncols)

    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex='none',
                             sharey='none',
                             figsize=(6 * ncols, 6 * nrows))
    sns.set(color_codes=True)

    row_index = 0
    col_index = 0
    for i in range(ncols * nrows):
        if nrows == 1 and ncols == 1:
            ax = axes
        elif nrows == 1 and ncols > 1:
            ax = axes[col_index]
        elif nrows > 1 and ncols == 1:
            ax = axes[row_index]
        else:
            ax = axes[row_index, col_index]

        if i < nplots:
            data = df.loc[df[sample_column] == samples[i], :].dropna()
            if data.shape[0] == 0:
                continue

            subtitle = ""
            if subtitles is not None:
                subtitle = subtitles[i]

            sns.despine(fig=fig, ax=ax)

            if plottype == "line":
                sns.lineplot(data=data,
                             x=x,
                             y=y,
                             hue=hue,
                             style=hue,
                             color=color,
                             alpha=alpha,
                             markers=["o"] * len(data[hue].unique()) if hue in data else "o",
                             markeredgewidth=0.0,
                             palette=palette,
                             ax=ax)
            elif plottype == "scatter":
                sns.scatterplot(data=data,
                                x=x,
                                y=y,
                                hue=hue,
                                color=color,
                                alpha=alpha,
                                s=10,
                                palette=palette,
                                ax=ax)
            elif plottype == "multi":
                g = sns.lineplot(data=data,
                                 x=x,
                                 y=y,
                                 hue=hue,
                                 color="black",
                                 alpha=1,
                                 palette=palette,
                                 ax=ax)
                ax.set(yscale="log")
                ax2 = ax.twinx()

                sns.scatterplot(data=data,
                                x=x,
                                y=y2,
                                hue=hue,
                                color=color,
                                alpha=alpha,
                                s=10,
                                palette=palette,
                                ax=ax2)

                ax2.set_ylabel(ylabel2,
                               color=color,
                               fontsize=10,
                               fontweight='bold')

            ax.set_xlabel(xlabel,
                          fontsize=10,
                          fontweight='bold')
            ax.set_ylabel(ylabel,
                          fontsize=10,
                          fontweight='bold')
            ax.set_title("{}{}".format(samples[i], subtitle),
                         fontsize=14,
                         fontweight='bold')


        else:
            ax.set_axis_off()

        col_index += 1
        if col_index > (ncols - 1):
            col_index = 0
            row_index += 1

    fig.suptitle(title,
                 fontsize=16,
                 color="#000000",
                 weight='bold')

    outpath = os.path.join(args.output, filename)
    fig.savefig(outpath)
    print("\tSaved figure: {}".format(os.path.basename(outpath)))
    plt.close()


#############################
############ STATS ##########
#############################
stats_data = []
for pool, metrics_path, log_path in zip(args.pools, args.metrics, args.log):
    # Load the metrics.
    metrics_df = pd.read_csv(metrics_path, sep=",")
    metrics_df.columns = ["key", "value"]
    metrics_df["value"] = metrics_df["value"].astype(float)
    metrics_df["id"] = pool + "_" + pathlib.PurePath(metrics_path).parent.name

    # Add the metrics from the log file.
    log_df = pd.DataFrame(parse_log_file(log_path), index=[0]).T.reset_index(drop=False)
    log_df.columns = ["key", "value"]
    metrics_df["id"] = pool + "_" + pathlib.PurePath(log_path).parent.name

    # Save the data.
    stats_data.append(metrics_df)
    stats_data.append(log_df)

stats_df = pd.concat(stats_data, axis=1)
barplot(
    df=stats_df,
    x="id",
    y="value",
    panel_column="key",
    filename="stats.png"
)
del stats_data, stats_df

##########################################
############ TRAINING PROCEDURE ##########
##########################################
training_procedure_data = []
for pool, raw_h5_path in zip(args.pools, args.raw_h5):
    hf = h5py.File(raw_h5_path, 'r')
    training_procedure_data.append(pd.DataFrame({
        "epoch": np.array(hf.get('metadata/learning_curve_train_epoch')),
        "loss": np.array(hf.get('metadata/learning_curve_train_elbo')),
        "id": pool + "_" + pathlib.PurePath(raw_h5_path).parent.name,
        "group": "Train"
    }))
    training_procedure_data.append(pd.DataFrame({
        "epoch": np.array(hf.get('metadata/learning_curve_test_elbo')),
        "loss": np.array(hf.get('metadata/learning_curve_test_epoch')),
        "id": pool + "_" + pathlib.PurePath(raw_h5_path).parent.name,
        "group": "Test"
    }))
    hf.close()

training_procedure_df = pd.concat(training_procedure_data, axis=1)
plot_per_sample(
    df=training_procedure_df,
    plottype="line",
    x="epoch",
    y="loss",
    hue="group",
    xlabel="Epoch",
    ylabel="ELBO",
    title="Progress of the training procedure",
    filename="training_procedure.png"
)
del training_procedure_data, training_procedure_df

########################################
############ CELL PROBABILITY ##########
########################################
cell_probability_data = []
for pool, posterior_h5_path in zip(args.pools, args.posterior_h5):
    hf = h5py.File(posterior_h5_path, 'r')
    cell_probability_data.append(pd.DataFrame({
        "d": np.array(hf.get('droplet_latents_map/d')),
        "p": np.array(hf.get('droplet_latents_map/p')),
        "id": pool + "_" + pathlib.PurePath(posterior_h5_path).parent.name
    }))
    hf.close()

cell_probability_df = pd.concat(cell_probability_data, axis=1)
# TODO: reorder based on cell size. I tried this but then the cell procentage scatterplot does not match the
#  CellBender output figure. If I don't reorder, then the UMI line does not match. I think I need to
#  order the UMI line but not the counts but that feels wrong.
plot_per_sample(
    df=cell_probability_df,
    plottype="multi",
    x="index",
    y="d",
    y2="p",
    color="red",
    alpha=0.3,
    xlabel="Barcode index",
    ylabel="UMI counts",
    ylabel2="Cell probability",
    title="Determination of which barcodes contain cells",
    filename="cell_probability.png"
)
del cell_probability_data, cell_probability_df

############################################
############ LATENT GENE ENCODING ##########
############################################
latent_gene_encoding_data = []
for pool, posterior_h5_path in zip(args.pools, args.posterior_h5):
    hf = h5py.File(posterior_h5_path, 'r')

    # Extract the data.
    p = np.array(hf.get('droplet_latents_map/p'))
    z = np.array(hf.get('droplet_latents_map/z'))

    # Calculate PCA.
    A = torch.as_tensor(z[p >= 0.5]).float()
    U, S, V = torch.pca_lowrank(A)
    z_pca = torch.matmul(A, V[:, :2]).to_numpy()

    # Save data.
    latent_gene_encoding_data.append(pd.DataFrame({
        "PC0": z_pca[:, 0],
        "PC1": z_pca[:, 1],
        "id": pool + "_" + pathlib.PurePath(posterior_h5_path).parent.name
    }))
    hf.close()

latent_gene_encoding_df = pd.concat(latent_gene_encoding_data, axis=1)
plot_per_sample(
    df=latent_gene_encoding_df,
    plottype="scatter",
    x="PC0",
    y="PC1",
    alpha=0.3,
    xlabel="PC 0",
    ylabel="PC 1",
    title="PCA of latent encoding of cell gene expression",
    filename="latent_gene_encoding"
)
del latent_gene_encoding_data, latent_gene_encoding_df

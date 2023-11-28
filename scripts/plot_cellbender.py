#!/usr/bin/env python
# Adapted from https://github.com/broadinstitute/CellBender/blob/master/cellbender/remove_background/report.py and
# https://github.com/broadinstitute/CellBender/blob/master/cellbender/remove_background/data/io.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--samples", required=True, nargs="+", type=str, help="")
parser.add_argument("--input_h5", required=True, nargs="+", type=str, help="")
parser.add_argument("--settings", required=True, nargs="+", type=str, help="")
parser.add_argument("--raw_h5", required=True, nargs="+", type=str, help="")
parser.add_argument("--posterior_h5", required=True, nargs="+", type=str, help="")
parser.add_argument("--max_plots_per_page", required=False, default=5, type=int, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import tables
import scipy.sparse as sp
import h5py
import torch


setting_abbr = {
    "use_cuda": "cuda",
    "expected_cells": "expected_cell_count",
    "total_droplets_included": "total_droplets",
    "force_cell_umi_prior": "force_cell_umi_prior",
    "force_empty_umi_prior": "force_empty_umi_prior",
    "model": "model",
    "epochs": "epochs",
    "low_count_threshold": "low_count_threshold",
    "z_dim": "z_hidden_dims",
    "z_layers": "z_layers",
    "training_fraction": "training_fraction",
    "empty_drop_training_fraction": "fraction_empties",
    "ignore_features": "blacklisted_genes",
    "fpr": "fpr",
    "exclude_feature_types": "exclude_features",
    "projected_ambient_count_threshold": "ambient_counts_in_cells_low_limit",
    "learning_rate": "learning_rate",
    "checkpoint_mins": "checkpoint_min",
    "final_elbo_fail_fraction": "final_elbo_fail_fraction",
    "epoch_elbo_fail_fraction": "epoch_elbo_fail_fraction",
    "num_training_tries": "num_training_tries",
    "learning_rate_retry_mult": "learning_rate_retry_mult",
    "posterior_batch_size": "posterior_batch_size",
    # "posterior_regularization": "posterior_regularization",
    "alpha": "prq_alpha",
    "q": "cdf_threshold_q",
    # "estimator": "estimator",
    "estimator_multiple_cpu": "use_multiprocessing_estimation",
    "constant_learning_rate": "constant_learning_rate",
    "debug": "debug"
}


def get_matrix_from_cellranger_h5(filename):
    # Detect CellRanger version.
    with tables.open_file(filename, 'r') as f:
        cellranger_version = 2
        try:
            # This works for version 3 but not for version 2.
            getattr(f.root.matrix, 'features')
            cellranger_version = 3
        except tables.NoSuchNodeError:
            pass

    with tables.open_file(filename, 'r') as f:
        # Initialize empty lists.
        csc_list = []

        # CellRanger v2:
        # Each group in the table (other than root) contains a genome,
        # so walk through the groups to get data for each genome.
        if cellranger_version == 2:
            for group in f.walk_groups():
                try:
                    # Read in data for this genome, and put it into a
                    # scipy.sparse.csc.csc_matrix
                    data = getattr(group, 'data').read()
                    indices = getattr(group, 'indices').read()
                    indptr = getattr(group, 'indptr').read()
                    shape = getattr(group, 'shape').read()
                    csc_list.append(sp.csc_matrix((data, indices, indptr), shape=shape))
                except tables.NoSuchNodeError:
                    # This exists to bypass the root node, which has no data.
                    pass

        # CellRanger v3:
        # There is only the 'matrix' group.
        elif cellranger_version == 3:
            data = getattr(f.root.matrix, 'data').read()
            indices = getattr(f.root.matrix, 'indices').read()
            indptr = getattr(f.root.matrix, 'indptr').read()
            shape = getattr(f.root.matrix, 'shape').read()
            csc_list.append(sp.csc_matrix((data, indices, indptr), shape=shape))

    # Put the data together (possibly from several genomes for v2 datasets).
    count_matrix = sp.vstack(csc_list, format='csc')
    count_matrix = count_matrix.transpose().tocsr()

    return count_matrix


def plot_settings(ax, settings, sample):
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame({setting_abbr[key]: str(value) for key, value in settings.items() if key in setting_abbr}, index=[0]).T
    max_key_length = max([len(str(key)) for key in df.index])
    max_value_length = max([len(str(value)) for value in df[0]])
    total_length = max_key_length + max_value_length
    table = ax.table(cellText=df.values, colWidths=[total_length / max_key_length, total_length / max_value_length], rowLabels=df.index, loc='center', edges='open')
    table.auto_set_column_width(col=list(range(len(df.columns))))
    table.scale(1, 0.8)
    ax.set_title(sample)


def plot_train_error(ax, loss):
    try:
        ax.plot(loss['train']['elbo'], '.--', label='Train')

        # Plot the test error, if there was held-out test data.
        if 'test' in loss.keys():
            if len(loss['test']['epoch']) > 0:
                ax.plot(loss['test']['epoch'],
                         loss['test']['elbo'], 'o:', label='Test')
                ax.legend()

        ylim_low = max(loss['train']['elbo'][0], loss['train']['elbo'][-1] - 2000)
        try:
            ylim_high = max(max(loss['train']['elbo']), max(loss['test']['elbo']))
        except ValueError:
            ylim_high = max(loss['train']['elbo'])
        ylim_high = ylim_high + (ylim_high - ylim_low) / 20
        ax.set_ylim([ylim_low, ylim_high])
    except:
        print('Error plotting the learning curve. Skipping.')
        pass

    ax.set_xlabel('Epoch')
    ax.set_ylabel('ELBO')
    ax.set_title('Progress of the training procedure')


def plot_barcodes_and_inferred_cell_prop(ax, umi_counts, p):
    count_order = np.argsort(umi_counts)[::-1]
    ax.semilogy(umi_counts[count_order], color='black')
    ax.set_ylabel('UMI counts')
    ax.set_xlabel('Barcode index, sorted by UMI count')
    if p is not None:  # The case of a simple model.
        right_ax = ax.twinx()
        right_ax.plot(p[count_order], '.:', color='red', alpha=0.3, rasterized=True)
        right_ax.set_ylabel('Cell probability', color='red')
        right_ax.set_ylim([-0.05, 1.05])
        ax.set_title('Determination of which barcodes contain cells')
    else:
        ax.set_title('The subset of barcodes used for training')


def plot_latent_encoding_z(ax, p, z):
    if p is None:
        p = np.ones(z.shape[0])

    # Do PCA on the latent encoding z.
    A = torch.as_tensor(z[p >= 0.5]).float()
    U, S, V = torch.pca_lowrank(A)
    z_pca = torch.matmul(A, V[:, :2])

    # Plot the latent encoding via PCA.
    ax.plot(z_pca[:, 0], z_pca[:, 1], '.', ms=3, color='black', alpha=0.3, rasterized=True)
    ax.set_ylabel('PC 1')
    ax.set_xlabel('PC 0')
    ax.set_title('PCA of latent encoding of gene expression in cells')


def plot_report(input_files, suffix="1"):
    nrows = len(input_files)

    fig, axs = plt.subplots(nrows, 4, figsize=(24, 6 * nrows), gridspec_kw={"width_ratios": [0.1, 0.3, 0.3, 0.3]})
    if nrows == 1:
        axs = axs[np.newaxis, ...]

    for row_index, (sample, input_h5_path, settings_path, raw_h5_path, posterior_h5_path) in enumerate(input_files):
        print("\tRow {} - sample: {}:".format(row_index, sample))
        print("\t  --input_h5 {}".format(input_h5_path))
        print("\t  --settings {}".format(settings_path))
        print("\t  --raw_h5 {}".format(raw_h5_path))
        print("\t  --posterior_h5 {}".format(posterior_h5_path))
        print("")

        fh = open(settings_path)
        settings = json.load(fh)
        fh.close()
        plot_settings(
            ax=axs[row_index, 0],
            settings=settings,
            sample=sample
        )

        hf = h5py.File(raw_h5_path, 'r')
        loss = {'train': {'epoch': np.array(hf.get('metadata/learning_curve_train_epoch')),
                          'elbo': np.array(hf.get('metadata/learning_curve_train_elbo'))},
                'test': {'epoch': np.array(hf.get('metadata/learning_curve_test_epoch')),
                         'elbo': np.array(hf.get('metadata/learning_curve_test_elbo'))},
                'learning_rate': {'epoch': np.array(hf.get('metadata/learning_rate_epoch')),
                                  'value': np.array(hf.get('metadata/learning_rate_value'))}}
        analyzed_barcode_inds = np.array(hf.get('metadata/barcodes_analyzed_inds'))
        analyzed_gene_inds = np.array(hf.get('metadata/features_analyzed_inds'))
        hf.close()

        plot_train_error(
            ax=axs[row_index, 1],
            loss=loss
        )

        hf = h5py.File(posterior_h5_path, 'r')
        p = np.array(hf.get('droplet_latents_map/p'))
        z = np.array(hf.get('droplet_latents_map/z'))
        hf.close()

        count_matrix = get_matrix_from_cellranger_h5(input_h5_path)
        trimmed_bc_matrix = count_matrix[analyzed_barcode_inds, :].tocsc()
        trimmed_matrix = trimmed_bc_matrix[:, analyzed_gene_inds].tocsr()
        counts = np.array(trimmed_matrix.sum(axis=1)).squeeze()

        plot_barcodes_and_inferred_cell_prop(
            ax=axs[row_index, 2],
            umi_counts=counts,
            p=p
        )

        plot_latent_encoding_z(
            ax=axs[row_index, 3],
            p=p,
            z=z
        )

    fig.tight_layout()
    outfile = os.path.join(args.out, 'CellBender_report.{}.png'.format(suffix))
    plt.savefig(outfile, bbox_inches="tight")
    return outfile


################################################################################

input_files = list(zip(args.samples, args.input_h5, args.settings, args.raw_h5, args.posterior_h5))
input_files.sort()

nplots = len(input_files)

start_indices = [start_index for start_index in range(0, nplots, args.max_plots_per_page)]
n_figures = len(start_indices)
start_indices.append(start_indices[-1] + args.max_plots_per_page)

print("Plotting")
for i in range(n_figures):
    reportfile = plot_report(
        input_files=input_files[start_indices[i]:start_indices[i + 1]],
        suffix=str(i)
    )
    print("\tSaved {}".format(reportfile))
    print("")

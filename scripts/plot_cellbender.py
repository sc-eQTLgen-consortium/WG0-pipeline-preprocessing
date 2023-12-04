#!/usr/bin/env python
# Adapted from https://github.com/broadinstitute/CellBender/blob/master/cellbender/remove_background/report.py
import argparse
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("--samples", required=True, nargs="+", type=str, help="")
parser.add_argument("--input_h5", required=True, nargs="+", type=str, help="")
parser.add_argument("--settings", required=True, nargs="+", type=str, help="")
parser.add_argument("--raw_h5", required=True, nargs="+", type=str, help="")
parser.add_argument("--max_plots_per_page", required=False, default=5, type=int, help="")
parser.add_argument("--out", required=True, type=str, help="The output directory where results will be saved.")
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out, exist_ok=True)

from cellbender.remove_background.downstream import load_anndata_from_input, load_anndata_from_input_and_output
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import scipy.stats
import torch


settings_info = {
    "use_cuda": ("cuda", None),
    "expected_cells": ("expected_cell_count", None),
    "total_droplets_included": ("total_droplets", None),
    "force_cell_umi_prior": ("force_cell_umi_prior", None),
    "force_empty_umi_prior": ("force_empty_umi_prior", None),
    "model": ("model", "full"),
    "epochs": ("epochs", 150),
    "low_count_threshold": ("low_count_threshold", 5),
    "z_dim":( "z_hidden_dims", 64),
    "z_layers": ("z_layers", [512]),
    "training_fraction": ("training_fraction", 0.9),
    "empty_drop_training_fraction": ("fraction_empties", 0.2),
    "ignore_features": ("blacklisted_genes", []),
    "fpr": ("fpr", [0.01]),
    "exclude_feature_types": ("exclude_features", []),
    "projected_ambient_count_threshold": ("ambient_counts_in_cells_low_limit", 0.1),
    "learning_rate": ("learning_rate", 1e-4),
    "checkpoint_mins": ("checkpoint_min", 7.),
    "final_elbo_fail_fraction": ("final_elbo_fail_fraction", None),
    "epoch_elbo_fail_fraction": ("epoch_elbo_fail_fraction", None),
    "num_training_tries": ("num_training_tries", 1),
    "learning_rate_retry_mult": ("learning_rate_retry_mult", 0.2),
    "posterior_batch_size": ("posterior_batch_size", 128),
    # "posterior_regularization": ("posterior_regularization", None),
    "alpha": ("prq_alpha", None),
    "q": ("cdf_threshold_q", None),
    # "estimator": ("estimator", None),
    "estimator_multiple_cpu": ("use_multiprocessing_estimation", False),
    "constant_learning_rate": ("constant_learning_rate", False),
    "debug": ("debug", False)
}

def assess_overall_count_removal(adata):
    cells = (adata.obs['cell_probability'] > 0.5)
    initial_counts = adata.layers['raw'][cells].sum()
    removed_counts = initial_counts - adata.layers['cellbender'][cells].sum()
    removed_percentage = removed_counts / initial_counts * 100

    estimated_ambient_per_droplet = np.exp(adata.uns['empty_droplet_size_lognormal_loc']).item()
    expected_fraction_removed_from_cells = estimated_ambient_per_droplet * cells.sum() / initial_counts

    fpr = adata.uns['target_false_positive_rate'].item()  # this is an np.ndarray with one element
    cohort_mode = False
    if type(fpr) != float:
        cohort_mode = True

    if not cohort_mode:
        expected_percentage = (expected_fraction_removed_from_cells + fpr) * 100
    else:
        expected_percentage = expected_fraction_removed_from_cells * 100

    if np.abs(expected_percentage - removed_percentage) <= 0.5:
        return {"overall_count_removal": "great"}
    elif np.abs(expected_percentage - removed_percentage) <= 1:
        return {"overall_count_removal": "decent"}
    elif np.abs(expected_percentage - removed_percentage) <= 5:
        if removed_percentage < expected_percentage:
            return {"overall_count_removal": "bit less"}
        else:
            return {"overall_count_removal": "bit more"}
    elif removed_percentage - expected_percentage > 5:
        return {"overall_count_removal": "more"}
    elif expected_percentage - removed_percentage > 5:
        return {"overall_count_removal": "fewer"}

def assess_count_removal_per_gene(adata,
                                  raw_full_adata,
                                  input_layer_key='raw',
                                  r_squared_cutoff=0.5):
    # how well does it correlate with our expectation about the ambient RNA profile?
    cells = (adata.obs['cell_probability'] > 0.5)
    counts = np.array(raw_full_adata.X.sum(axis=1)).squeeze()
    clims = [adata.obs[f'n_{input_layer_key}'][~cells].mean() / 2,
             np.percentile(adata.obs[f'n_{input_layer_key}'][cells].values, q=2)]
    # if all are called "cells" then clims[0] will be a nan
    if np.isnan(clims[0]):
        clims[0] = counts.min()
    if 'approximate_ambient_profile' in adata.uns.keys():
        approximate_ambient_profile = adata.uns['approximate_ambient_profile']
    else:
        empty_count_matrix = raw_full_adata[(counts > clims[0]) & (counts < clims[1])].X
        if empty_count_matrix.shape[0] > 100:
            approximate_ambient_profile = np.array(raw_full_adata[(counts > clims[0])
                                                                  & (counts < clims[1])].X.mean(axis=0)).squeeze()
        else:
            # a very rare edge case I've seen once
            approximate_ambient_profile = np.array(raw_full_adata[counts < clims[1]].X.mean(axis=0)).squeeze()
        approximate_ambient_profile = approximate_ambient_profile / approximate_ambient_profile.sum()
    y = adata.var['n_removed'] / adata.var['n_removed'].sum()

    cutoff = 1e-6
    logic = np.logical_not((approximate_ambient_profile < cutoff) | (y < cutoff))
    r_squared_result = scipy.stats.pearsonr(np.log(approximate_ambient_profile[logic]),
                                            np.log(y[logic]))
    if hasattr(r_squared_result, 'statistic'):
        # scipy version 1.9.0+
        r_squared = r_squared_result.statistic
    else:
        r_squared = r_squared_result[0]

    if r_squared > r_squared_cutoff:
        return {"per_gene_removal": "correct"}
    else:
        return {"per_gene_removal": "incorrect"}

def assess_learning_curve(adata,
                          spike_size: float = 0.5,
                          deviation_size: float = 0.25,
                          monotonicity_cutoff: float = 0.1):
    warnings = {
        "short_run": False,
        "large_spikes_in_train": np.nan,
        "large_deviation_in_train": np.nan,
        "low_end_in_train": np.nan,
        "non_monotonic": np.nan,
        "backtracking": np.nan,
        "runaway_test": np.nan,
        "non_convergence": np.nan
    }

    if 'learning_curve_train_elbo' not in adata.uns.keys():
        return warnings

    if adata.uns['learning_curve_train_epoch'][-1] < 50:
        warnings["short_run"] = True
        return warnings

    train_elbo_min_max = np.percentile(adata.uns['learning_curve_train_elbo'], q=[5, 95])
    train_elbo_range = train_elbo_min_max.max() - train_elbo_min_max.min()

    # look only from epoch 45 onward for spikes in train ELBO
    warnings["large_spikes_in_train"] = np.any((adata.uns['learning_curve_train_elbo'][46:]
                                                - adata.uns['learning_curve_train_elbo'][45:-1])
                                               < -train_elbo_range * spike_size)

    second_half_train_elbo = (adata.uns['learning_curve_train_elbo']
                              [(len(adata.uns['learning_curve_train_elbo']) // 2):])
    warnings["large_deviation_in_train"] = np.any(second_half_train_elbo
                                                  < np.median(second_half_train_elbo)
                                                  - train_elbo_range * deviation_size)

    half = len(adata.uns['learning_curve_train_elbo']) // 2
    threequarter = len(adata.uns['learning_curve_train_elbo']) * 3 // 4
    typical_end_variation = np.std(adata.uns['learning_curve_train_elbo'][half:threequarter])
    warnings["low_end_in_train"] = (adata.uns['learning_curve_train_elbo'][-1]
                                    < adata.uns['learning_curve_train_elbo'].max() - 5 * typical_end_variation)

    # look only from epoch 45 onward for spikes in train ELBO
    non_monotonicity = ((adata.uns['learning_curve_train_elbo'][46:]
                         - adata.uns['learning_curve_train_elbo'][45:-1])
                        < -3 * typical_end_variation).sum() / len(adata.uns['learning_curve_train_elbo'])
    warnings["non_monotonic"] = (non_monotonicity > monotonicity_cutoff)

    def windowed_cumsum(x, n=20):
        return np.array([np.cumsum(x[i:(i + n)])[-1] for i in range(len(x) - n)])

    windowsize = 20
    tracking_trace = windowed_cumsum(adata.uns['learning_curve_train_elbo'][1:]
                                     - adata.uns['learning_curve_train_elbo'][:-1],
                                     n=windowsize)
    big_dip = -1 * (adata.uns['learning_curve_train_elbo'][-1]
                    - adata.uns['learning_curve_train_elbo'][5]) / 10
    warnings["backtracking"] = (tracking_trace.min() < big_dip)

    halftest = len(adata.uns['learning_curve_test_elbo']) // 2
    threequartertest = len(adata.uns['learning_curve_test_elbo']) * 3 // 4
    typical_end_variation_test = np.std(adata.uns['learning_curve_test_elbo'][halftest:threequartertest])
    warnings["runaway_test"] = (adata.uns['learning_curve_test_elbo'][-1]
                                < adata.uns['learning_curve_test_elbo'].max() - 4 * typical_end_variation_test)

    warnings["non_convergence"] = (np.mean([adata.uns['learning_curve_train_elbo'][-1]
                                            - adata.uns['learning_curve_train_elbo'][-2],
                                            adata.uns['learning_curve_train_elbo'][-2]
                                            - adata.uns['learning_curve_train_elbo'][-3]])
                                   > 2 * typical_end_variation)

    return warnings


def plot_table(ax, values, colors=None, title=""):
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame(values, index=[0]).T
    max_key_length = max([len(str(key)) for key in df.index])
    max_value_length = max([len(str(value)) for value in df[0]])
    total_length = max_key_length + max_value_length
    table = ax.table(cellText=df.values, colWidths=[total_length / max_key_length, total_length / max_value_length], rowLabels=df.index, loc='center', edges='open')
    table.auto_set_column_width(col=list(range(len(df.columns))))
    table.scale(1, 0.8)
    if colors is not None:
        for row_index in range(len(values)):
            key = table[(row_index, -1)].get_text()
            key.set_color(colors[key.get_text()])
            value = table[(row_index, 0)].get_text()
            value.set_color(colors[key.get_text()])
    ax.set_title(title)


def plot_train_error(adata, ax):
    if 'learning_curve_train_elbo' not in adata.uns.keys():
        print('No learning curve recorded!')
        ax.set_axis_off()

    try:
        ax.plot(adata.uns['learning_curve_train_epoch'],
                adata.uns['learning_curve_train_elbo'], '.--', label='Train')
        try:
            ax.plot(adata.uns['learning_curve_test_epoch'],
                    adata.uns['learning_curve_test_elbo'], 'o:', label='Test')
            ax.legend()
        except Exception:
            pass

        ylim_low = max(adata.uns['learning_curve_train_elbo'][0], adata.uns['learning_curve_train_elbo'][-1] - 2000)
        try:
            ylim_high = max(max(adata.uns['learning_curve_train_elbo']), max(adata.uns['learning_curve_test_elbo']))
        except ValueError:
            ylim_high = max(adata.uns['learning_curve_train_elbo'])
        ylim_high = ylim_high + (ylim_high - ylim_low) / 20
        ax.set_ylim([ylim_low, ylim_high])
    except:
        print('Error plotting the learning curve. Skipping.')
        pass

    ax.set_xlabel('Epoch')
    ax.set_ylabel('ELBO')
    ax.set_title('Progress of the training procedure')


def plot_barcodes_and_inferred_cell_prop(adata, ax):
    in_counts = np.array(adata.layers['raw'][:, adata.var['cellbender_analyzed']].sum(axis=1)).squeeze()
    order = np.argsort(in_counts)[::-1]
    ax.semilogy(in_counts[order], color='black')
    ax.set_ylabel('UMI counts')
    ax.set_xlabel('Barcode index, sorted by UMI count')
    right_ax = ax.twinx()
    right_ax.plot(adata.obs['cell_probability'][order].values, '.:', color='red', alpha=0.3, rasterized=True)
    right_ax.set_ylabel('Cell probability', color='red')
    right_ax.set_ylim([-0.05, 1.05])
    ax.set_title('Determination of which barcodes contain cells')


def plot_latent_encoding_z(adata, ax):
    cells = (adata.obs['cell_probability'] > 0.5)
    adata.obsm['X_pca'] = pca_2d(adata.obsm['cellbender_embedding']).detach().numpy()

    # Plot the latent encoding via PCA.
    ax.plot(adata.obsm['X_pca'][:, 0][cells],
            adata.obsm['X_pca'][:, 1][cells],
            '.',
            ms=3,
            color='black',
            alpha=0.3,
            rasterized=True)
    ax.set_ylabel('PC 1')
    ax.set_xlabel('PC 0')
    ax.set_title('PCA of latent encoding of gene expression in cells')


def pca_2d(mat: np.ndarray) -> torch.Tensor:
    A = torch.as_tensor(mat).float()
    U, S, V = torch.pca_lowrank(A)
    return torch.matmul(A, V[:, :2])

def plot_report(input_files, suffix="1"):
    nrows = len(input_files)

    fig, axs = plt.subplots(nrows, 5, figsize=(30, 6 * nrows), gridspec_kw={"width_ratios": [0.05, 0.05, 0.3, 0.3, 0.3]})
    if nrows == 1:
        axs = axs[np.newaxis, ...]

    for row_index, (sample, input_h5_path, settings_path, raw_h5_path) in enumerate(input_files):
        print("\tRow {} - sample: {}:".format(row_index, sample))
        print("\t  --input_h5 {}".format(input_h5_path))
        print("\t  --settings {}".format(settings_path))
        print("\t  --raw_h5 {}".format(raw_h5_path))
        print("")

        fh = open(settings_path)
        settings = {}
        colors = {}
        for key, value in json.load(fh).items():
            if key in settings_info:
                new_key, default_value = settings_info[key]
                settings[new_key] = str(value)
                if str(value) == str(default_value):
                    colors[new_key] = "black"
                else:
                    colors[new_key] = "red"
        fh.close()
        plot_table(
            ax=axs[row_index, 0],
            values=settings,
            colors=colors,
            title=sample
        )

        # load datasets, before and after CellBender
        adata = load_anndata_from_input_and_output(input_file=input_h5_path,
                                                   output_file=raw_h5_path,
                                                   analyzed_barcodes_only=True,
                                                   input_layer_key='raw',
                                                   truth_file=None)
        raw_full_adata = load_anndata_from_input(input_h5_path)

        # bit of pre-compute
        adata.var['n_removed'] = adata.var[f'n_raw'] - adata.var[f'n_cellbender']

        warnings = {}
        try:
            warnings.update(assess_overall_count_removal(adata))
        except ValueError:
            print('Skipping assessment over overall count removal. Presumably '
                  'this is due to including the whole dataset in '
                  '--total-droplets-included.')

        try:
            warnings.update(assess_learning_curve(adata))
        except Exception:
            pass

        # look at per-gene count removal
        warnings.update(assess_count_removal_per_gene(adata, raw_full_adata=raw_full_adata))
        has_warnings = False
        colors = {}
        for key, value in warnings.items():
            if value is True or value == "incorrect" or value == "more" or value == "fewer":
                has_warnings = True
                colors[key] = "red"
            else:
                colors[key] = "green"
        warnings["--------"] = "--------"
        warnings["Warnings"] = has_warnings
        colors["--------"] = "black"
        colors["Warnings"] = "red" if has_warnings else "green"

        plot_table(
            ax=axs[row_index, 1],
            values=warnings,
            colors=colors,
            title="Warnings"
        )
        plot_train_error(
            adata=adata,
            ax=axs[row_index, 2]
        )
        plot_barcodes_and_inferred_cell_prop(
            adata=adata,
            ax=axs[row_index, 3]
        )
        plot_latent_encoding_z(
            adata=adata,
            ax=axs[row_index, 4]
        )

    fig.tight_layout()
    outfile = os.path.join(args.out, 'CellBender_report.{}.png'.format(suffix))
    plt.savefig(outfile, bbox_inches="tight")
    return outfile


################################################################################

input_files = list(zip(args.samples, args.input_h5, args.settings, args.raw_h5))
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

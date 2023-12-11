#!/usr/bin/env python
import numpy as np
import pandas as pd
import itertools
import json
import math
import os

# Add trailing /.
if not config["inputs"]["scripts_dir"].endswith("/"):
    config["inputs"]["scripts_dir"] += "/"
if not config["inputs"]["fastq_dir"].endswith("/"):
    config["inputs"]["fastq_dir"] += "/"
if not config["refs"]["ref_dir"].endswith("/"):
    config["refs"]["ref_dir"] += "/"
if not config["outputs"]["output_dir"].endswith("/"):
    config["outputs"]["output_dir"] += "/"

# Check if the singularity image exists.
if not os.path.exists(config["inputs"]["singularity_image"]):
    logger.info("Error, the singularity image does not exist.\n\nExiting.")
    exit("MissingSIFFile")

# Check if the samplesheet exists.
if not os.path.exists(config["inputs"]["samplesheet_path"]):
    logger.info("Error, the samplesheet file does not exist.\n\nExiting.")
    exit("MissingSampleSheetFile")

# Loading the input samplesheet.
logger.info("Loading the input samplesheet")
SAMPLE_DF = pd.read_csv(config["inputs"]["samplesheet_path"], sep="\t", dtype=str)
SAMPLE_DF.fillna("NA", inplace=True)
SAMPLE_DF.index = SAMPLE_DF["Sample"]

# Check for missing columns.
missing_columns = [column for column in ["Fastq", "Sample"] if not column in SAMPLE_DF.columns]
if len(missing_columns) > 0:
    logger.info("\tError, missing columns {} in samplesheet file for the selected methods.".format(", ".join(missing_columns)))
    exit()

# Check if the input samplesheet is valid.
samplesheet_is_valid = True
for column in ["Fastq", "Sample"]:
    if not SAMPLE_DF[column].is_unique:
        logger.info("\tYour {} column contains duplicates, please make sure all values are unique.".format(column))
        samplesheet_is_valid = False

if not samplesheet_is_valid:
    logger.info("\n\nExiting.")
    exit("InvalidSampleSheet")

logger.info("\tValid.")
SAMPLES = SAMPLE_DF["Sample"].values.tolist()
FASTQ_DICT = dict(zip(SAMPLE_DF["Sample"], SAMPLE_DF["Fastq"]))

# Check if the reference transcriptome exists.
if not os.path.isdir(config["refs"]["ref_dir"] + config["refs_extra"]["transcriptome_path"]):
    logger.info("Could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["refs"]["ref_dir"] + config["refs_extra"]["transcriptome_path"]))
    exit("MissingReferenceFile")

# Create the aggregation file using for CellRanger_aggr.
# This can be harded coded since CellRanger does not require manual selection.
if len(SAMPLES) > 1:
    aggregation_path = config["outputs"]["output_dir"] + "CellRanger/aggregation.csv"
    if not os.path.exists(aggregation_path):
        if not os.path.isdir(config["outputs"]["output_dir"] + "CellRanger"):
            os.mkdir(config["outputs"]["output_dir"] + "CellRanger")
        aggregation_csv = pd.DataFrame({"sample_id": SAMPLES,
                                        "molecule_h5": [config["outputs"]["output_dir"] + "/CellRanger/{sample}/outs/molecule_info.h5".format(sample=sample) for sample in SAMPLES]})
        aggregation_csv.to_csv(aggregation_path, sep=",", index=False)

def process_manual_selection_method(name, settings=None, extra_settings=None, settings_dtype=None, max_manual_runs=25):
    """
    This function allows for methods to be run with different settings in parallel without overwriting previous results.
    """
    if settings is None:
        settings = []
    if extra_settings is None:
        extra_settings = []
    small_dtype = {"Sample": object}
    full_dtype = {"Sample": object, "Run": int, "FINISHED": bool, "PASSED": bool}
    if settings_dtype is not None:
        small_dtype.update(settings_dtype)
        full_dtype.update(settings_dtype)
    all_settings = settings + extra_settings

    if not os.path.isdir(config["outputs"]["output_dir"] + "manual_selection"):
        os.mkdir(config["outputs"]["output_dir"] + "manual_selection")

    # Step 1. Load the manual selection file and validate the content. Check if all samples have one
    # run that is finished and passed.
    man_select_path = config["outputs"]["output_dir"] + "manual_selection/{name}_manual_selection.tsv".format(name=name)
    if os.path.exists(man_select_path):
        select_df = pd.read_csv(man_select_path, dtype=full_dtype, sep="\t", header=0, index_col=None)
        for index, row in select_df.iterrows():
            settings_path = config["outputs"]["output_dir"] + "{name}/{sample}Run{run}/{name}_settings.json".format(sample=row["Sample"], name=name, run=row["Run"])
            # Load the settings file and check if it matches.
            if not os.path.exists(settings_path):
                # We assume it hasn't started yet.
                continue

            # Load settings. I need to force the type on the output since we can have lists as parameters but those
            # are interpreted as strings since Pandas does not work with lists in a DataFrame.
            fh = open(settings_path)
            used_settings = json.load(fh)
            used_settings_converted  = {}
            for parameter in row.index:
                if parameter not in used_settings:
                    continue
                value = used_settings[parameter]
                dtype = settings_dtype[parameter]
                if value is None and dtype == float:
                    value = np.nan
                used_settings_converted[parameter] = dtype(value)
            used_settings = pd.DataFrame(used_settings_converted, index=[0]).astype(settings_dtype)
            fh.close()

            # Validate settings.
            for setting in all_settings:
                if setting not in used_settings or (row[setting] != used_settings.loc[0, setting]):
                    if math.isnan(row[setting]) and math.isnan(used_settings.loc[0, setting]):
                        # Edge case where nan == nan is False.
                        continue
                    logger.info("\tError, output directory {name}/{sample}Run{run}/{name}_settings.json contains unexpected settings.".format(sample=row["Sample"], name=name, run=row["Run"]))
                    exit("UnexpectedSetting")

            # Check if the run was completed.
            results_path = config["outputs"]["output_dir"] + "{name}/{sample}Run{run}/{name}.done".format(sample=row["Sample"], name=name, run=row["Run"])
            select_df.loc[index, "FINISHED"] = os.path.exists(results_path)

        if select_df["FINISHED"].all():
            logger.info("\tAll expected output files are created.")
        else:
            logger.info("\tWaiting on expected output files.")

        # Check if each sample has exactly one FINISHED run with a PASSED flag, if so: we are done.
        passed_select_df = select_df.loc[(select_df["FINISHED"]) & (select_df["PASSED"]), :].copy()
        if passed_select_df.shape[0] == len(SAMPLES) and len(set(passed_select_df["Pool"].values).symmetric_difference(set(SAMPLES))) == 0:
            return True, [], dict(zip(passed_select_df["Pool"], passed_select_df["Run"])), {}

        if select_df["FINISHED"].all() and (passed_select_df.shape[0] < len(SAMPLES) or len(set(passed_select_df["Pool"].values).symmetric_difference(set(SAMPLES))) != 0):
            logger.info("\tPlease select one accepted run for all samples in the manual_selection file or add additional runs in the manual rerun file.")
        elif select_df["FINISHED"].all() and passed_select_df.shape[0] > len(SAMPLES):
            logger.info("\tPlease select one accepted run for each sample in the manual_selection file or add additional runs in the manual rerun file.")
        else:
            pass

        del passed_select_df
    else:
        # In case the manual_selection.tsv does not exist or got corrupted we can regenerate it based on the files we can find
        # in the output directory.
        select_data = []
        for sample in SAMPLES:
            for run_id in range(1, max_manual_runs):
                results_path = config["outputs"]["output_dir"] + "{name}/{sample}Run{run}/{name}.done".format(sample=sample, name=name, run=run_id)
                settings_path = config["outputs"]["output_dir"] + "{name}/{sample}Run{run}/{name}_settings.json".format(sample=sample, name=name, run=run_id)
                if os.path.exists(settings_path):
                    # Load settings.
                    fh = open(settings_path)
                    used_settings = json.load(fh)
                    fh.close()

                    # Validate settings.
                    valid = True
                    for setting in all_settings:
                        if setting not in used_settings.keys():
                            valid = False
                            break

                    if not valid:
                        logger.info("\tError, output directory {name}/{sample}Run{run}/{name}_settings.json contains unexpected settings.".format(sample=sample, name=name, run=run_id))
                        exit("UnexpectedSetting")

                    result_exists = os.path.exists(results_path)
                    select_data.append([sample, run_id] + [used_settings[setting] for setting in all_settings] + [result_exists, False])
        select_df = pd.DataFrame(select_data, columns=["Sample", "Run"] + all_settings + ["FINISHED", "PASSED"]).astype(full_dtype)

    # Step 2a. Generate a new settings data frame based on default settings.
    # Note that this only includes the sample and the settings (i.e. not Run, FINISHED, PASSED).
    # Also, by using lists we can generate all possible combinations of default settings. Therefore, make sure that
    # standard settings (i.e. settings that might change) are in a list while the extra settings (things the user shouldn't really change)
    # are not lists.
    default_select_data = []
    default_settings = [config[name.lower()][name.lower() + "_" + setting] for setting in settings]
    default_extra_settings = [config[name.lower() + "_extra"][setting] for setting in extra_settings]
    for sample in SAMPLES:
        for ds in list(itertools.product(*default_settings)):
            default_select_data.append([sample] + list(ds) + default_extra_settings)
    default_select_df = pd.DataFrame(default_select_data, columns=["Sample"] + all_settings).astype(small_dtype)
    select_df = select_df.merge(default_select_df, how="outer")
    del default_select_df

    # Step 2b. Add the manual rerun settings data frame.
    # Note that this only includes the sample and the settings (i.e. not Run, FINISHED, PASSED)
    man_rerun_path = config["outputs"]["output_dir"] + "manual_selection/{name}_manual_run.tsv".format(name=name)
    if os.path.exists(man_rerun_path):
        man_rerun_df = pd.read_csv(man_rerun_path, dtype=small_dtype, sep="\t", header=0, index_col=None)
        man_rerun_df = man_rerun_df.loc[man_rerun_df["Sample"] != "Example", :]

        select_df = select_df.merge(man_rerun_df, how="outer")
        del man_rerun_df
    else:
        pd.DataFrame([["Example"] + [config[name.lower() + "_extra"][setting] for setting in all_settings]], columns=["Sample"] + all_settings).to_csv(man_rerun_path, sep="\t", header=True, index=False)

    # Step 3. Assign unique output directories (run IDs) to each run and reformat the data frame.
    select_df = select_df.loc[:, ["Sample", "Run"] + all_settings + ["FINISHED", "PASSED"]]
    select_df[["FINISHED", "PASSED"]] = select_df[["FINISHED", "PASSED"]].fillna(False)
    for index, row in select_df.iterrows():
        # Fill in the missing Run IDs.
        if np.isnan(row["Run"]):
            max_sample_id = select_df.loc[select_df["Sample"] == row["Sample"], "Run"].max()
            if np.isnan(max_sample_id):
                max_sample_id = 0
            select_df.loc[index, "Run"] = max_sample_id + 1
    select_df = select_df.astype(full_dtype)

    # Step 4. Generate the input files.
    output_folders = []
    settings = {row["Sample"]: {} for _, row in select_df.iterrows()}
    for _, row in select_df.iterrows():
        output_folders.append(config["outputs"]["output_dir"] + "{name}/{sample}Run{run}/".format(sample=row["Sample"], name=name, run=row["Run"]))
        settings[row["Sample"]][str(row["Run"])] = row[all_settings].to_dict()

    # Step 5. Safe manual_selection_df.
    select_df.to_csv(man_select_path, sep="\t", header=True, index=False)

    return False, output_folders, {}, settings

#####################
######## ALL ########
#####################
input_files = []

logger.info("Running CellRanger.")
if len(SAMPLES) == 1:
    # One sample so no need to run CellRanger aggr.
    input_files.append(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/web_summary.html".format(sample=SAMPLES[0]))
else:
    # Multiple samples so we do need to run CellRanger aggr.
    input_files.append(config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/web_summary.html")

cellbender_passed = False
CELLBENDER_REPORT_INDICES = []
if config["settings"]["ambient_rna_correction"]:
    logger.info("Running CellBender on {}.".format("GPU" if config["settings"]["use_gpu"] else "CPU"))
    cellbender_passed, cellbender_output_folders, CELLBENDER_SELECTION, CELLBENDER_SETTINGS = process_manual_selection_method(
        name="CellBender",
        settings=[],
        extra_settings=["expected_cells", "total_droplets_included", "force_cell_umi_prior", "force_empty_umi_prior", "model", "epochs",
                        "low_count_threshold", "z_dim", "z_layers", "training_fraction", "empty_drop_training_fraction", "ignore_features",
                        "fpr", "exclude_feature_types", "projected_ambient_count_threshold", "learning_rate", "checkpoint_mins",
                        "final_elbo_fail_fraction", "epoch_elbo_fail_fraction", "num_training_tries", "learning_rate_retry_mult",
                        "posterior_batch_size", "alpha", "q", "estimator_multiple_cpu", "constant_learning_rate", "debug"],
        settings_dtype={
            "expected_cells": float, # int
            "total_droplets_included": float, # int
            "force_cell_umi_prior": float,
            "force_empty_umi_prior": float,
            "model": str,
            "epochs": int,
            "low_count_threshold": int,
            "z_dim": int,
            "z_layers": str, # [int]
            "training_fraction": float,
            "empty_drop_training_fraction": float,
            "ignore_features": str, # [int]
            "fpr": str, # [float]
            "exclude_feature_types": str, # [str]
            "projected_ambient_count_threshold": float,
            "learning_rate": float,
            "checkpoint_mins": float,
            "final_elbo_fail_fraction": float,
            "epoch_elbo_fail_fraction": float,
            "num_training_tries": int,
            "learning_rate_retry_mult": float,
            "posterior_batch_size": int,
            "alpha": float,
            "q": float,
            "estimator_multiple_cpu": bool,
            "constant_learning_rate": bool,
            "debug": bool
        }
    )

    if not cellbender_passed:
        for i in range(0, len(cellbender_output_folders), config["cellbender_extra"]["max_plots_per_page"]):
            CELLBENDER_REPORT_INDICES.append(i)
            input_files.extend([config["outputs"]["output_dir"] + "QC_figures/CellBender_report.{}.png".format(i)])


if not config["settings"]["ambient_rna_correction"] or cellbender_passed:
    # End point of rule combine_results.
    input_files = [config["outputs"]["output_dir"] + "Combine_Results/wg0_file_directories.tsv"]

    include: "includes/combine_results.smk"

rule all:
    input:
        files = input_files

# Import individual rules
include: "includes/cellranger.smk"
include: "includes/cellbender.smk"

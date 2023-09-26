#!/usr/bin/env python
import numpy as np
import pandas as pd
import os

# Get the samples.
SAMPLES = pd.read_csv(config["inputs"]["samplesheet_filepath"], sep="\t")
SAMPLES = SAMPLES.iloc[:2, :] #TODO: remove
SAMPLES.columns = ["Pool", "ID"]
SAMPLES["Pool"] = SAMPLES["Pool"].astype(str)
SAMPLES["ID"] = SAMPLES["ID"].astype(str)
SAMPLE_DICT = dict(zip(SAMPLES["Pool"], SAMPLES["ID"]))

# Determine which methods we will be running.
METHODS = ["CellRanger"]
if config["settings"]["ambient_rna_correction"]:
    logger.info("Running CellBender rules since 'ambient_rna_correction' is set to True.")
    METHODS.append("CellBender")

# Remove trailing / from paths.
if not config["inputs"]["fastqs"].endswith("/"):
    config["inputs"]["fastqs"] += "/"
if not config["outputs"]["output_dir"].endswith("/"):
    config["outputs"]["output_dir"] += "/"

# Create the aggregation file using for CellRanger_aggr.
# This can be harded coded since CellRanger does not require manual selection.
aggregation_csv = pd.DataFrame({"sample_id": SAMPLES["Pool"], "molecule_h5": [config["outputs"]["output_dir"] + pool + "/CellRanger_count/outs/molecule_info.h5" for pool in SAMPLES["Pool"]]})
if not os.path.isdir(config["outputs"]["output_dir"] + "CellRanger_aggr"):
    os.mkdir(config["outputs"]["output_dir"] + "CellRanger_aggr")
aggregation_csv_path = config["outputs"]["output_dir"] + "CellRanger_aggr/aggregation.csv"
aggregation_csv.to_csv(aggregation_csv_path, sep=",", index=False)

# Determine the settings to use for cellbender for each of the pools. By doing it like this I allow
# for default and manual jobs to be run at the same time. By using a unique output directory for each run
# I prevent having to delete previous runs when a run fails. This should allow the user to be able to select
# an old run as the best one if it turns out that the new run does not perform any better.
CELLBENDER_SETTINGS = {pool: {"outdir": "CellBenderDefaultRun1",
                              "epochs": config["cellbender"]["cellbender_epochs"],
                              "learning-rate": config["cellbender"]["cellbender_learning_rate"]}
                       for pool in SAMPLES["Pool"]}
if "CellBender" in METHODS and config["cellbender_manual"]["run_cellbender_manual"]:
    for pool, epochs, learning_rate in zip(config["cellbender_manual"]["cellbender_manual_pools"],
            config["cellbender_manual"]["cellbender_manual_epochs"],
            config["cellbender_manual"]["cellbender_manual_learning_rates"]):

        # Find a CellBender output directory that does not exist (e.g. CellBenderManualRun1, CellBenderManualRun2, etc.)
        found = False
        for i in range(1, 100):
            if not found and not os.path.exists(config["outputs"]["output_dir"] + str(pool) + "/CellBenderManualRun" + str(i)):
                CELLBENDER_SETTINGS[str(pool)] = {"outdir": "CellBenderManualRun" + str(i), "epochs": epochs, "learning-rate": learning_rate}
                found = True
                break

        # This should not happen. Why would anyone run 100 manual runs?!
        if not found:
            logger.info("Error, could not find valid output directory for CellBender manual run.")
            exit()

# Find all previous runs of CellBender. This enables us to plot previous runs in the CellBender_summary rule to compare the outputs.
CELLBENDER_OUTDIRS = []
if "CellBender" in METHODS:
    for pool in SAMPLES["Pool"]:
        # Look for old default runs.
        default_dir = config["outputs"]["output_dir"] + pool + "/CellBenderDefaultRun1"
        if os.path.exists(default_dir):
            CELLBENDER_OUTDIRS.append(default_dir)

        # Look for old manual runs.
        for i in range(1,100):
            manual_dir = config["outputs"]["output_dir"] + pool + "/CellBenderManualRun" + str(i)
            if os.path.exists(manual_dir):
                CELLBENDER_OUTDIRS.append(manual_dir)

        # Add the current runs.
        CELLBENDER_OUTDIRS.append(CELLBENDER_SETTINGS[pool]["outdir"])

#####################
######## ALL ########
#####################
def get_all_input(wildcards):
    """
    Main function of the snakemake. This will determine which output files we need to start on WG1 based on the user settings.
    """
    input_files = expand(config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/possorted_genome_bam.bam", pool=SAMPLES["Pool"])

    if "CellRanger" in METHODS:
        # Files the CellRanger files we need for WG1 or for CellBender.
        input_files.extend(expand(config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/filtered_feature_bc_matrix.h5", pool=SAMPLES["Pool"]))
        input_files.extend(expand(config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", pool=SAMPLES["Pool"]))

        # The aggregated report we want.
        input_files.append(config["outputs"]["output_dir"] + "CellRanger_aggr/outs/web_summary.html")

        # Construct the file_directories.tsv file that points to the files we will use for WG1.
        if "CellBender" not in METHODS:
            samplesheet_data = []
            for pool in SAMPLES["Pool"]:
                samplesheet_data.append([
                    pool,
                    config["outputs"]["output_dir"] + pool + "/CellRanger_count/outs/possorted_genome_bam.bam",
                    config["outputs"]["output_dir"] + pool + "/CellRanger_count/outs/filtered_feature_bc_matrix.h5",
                    config["outputs"]["output_dir"] + pool + "/CellRanger_count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
                ])

            logger.info("Saving file directories for WG1.")
            pd.DataFrame(samplesheet_data, columns=["Pool", "Bam", "Counts", "Barcodes"]).to_csv(config["outputs"]["output_dir"] + "wg0_file_directories.tsv", sep="\t", index=False)

    if "CellBender" in METHODS:
        # Create the manual selection directory.
        if not os.path.isdir(config["outputs"]["output_dir"] + "manual_selections"):
            os.mkdir(config["outputs"]["output_dir"] + "manual_selections")

        # Check if the manual select file exists.
        man_select_path = config["outputs"]["output_dir"] + "manual_selections/CellBender_manual_selection.tsv"
        empty_select_df = pd.DataFrame({"Pool": SAMPLES["Pool"], "Selection": np.nan})
        if os.path.exists(man_select_path):
            # The file exists, here we check if all the input is valid meaning the pipeline was run completely.
            logger.info("Read in the CellBender manual selection file.")
            selection = pd.read_csv(man_select_path, sep="\t", index_col=None)
            selection["Pool"] = selection["Pool"].astype(str)

            # If the selection df equals the empty one snakemake creates than there is no need to check any input (e.g. the user
            # did not change anything)
            if not selection.equals(empty_select_df):
                logger.info("The CellBender_manual_selection.tsv file has been updated, evaluating if selection is valid.")
                selection["Selection"] = selection["Selection"].astype(str)

                # Check if the pools still match.
                if len(set(SAMPLES["Pool"]).intersection(set(selection["Pool"]))) != SAMPLES.shape[0]:
                    logger.info("ERROR: the CellBender_manual_selection.tsv does not contain the same pools as the samplesheet.")
                    exit()

                samplesheet_data = []
                for _, row in selection.iterrows():
                    # Check if the CellRanger bam exists.
                    cell_ranger_bam_path = config["outputs"]["output_dir"] + row["Pool"] + "/CellRanger_count/outs/possorted_genome_bam.bam"
                    if not os.path.exists(cell_ranger_bam_path):
                        logger.info("ERROR: the '{}/CellRanger_count/outs/possorted_genome_bam.bam' does not exist.".format(row["Pool"]))
                        exit()

                    # Check if all the output files for WG1 exist in the selected output directory for this pool.
                    cellbender_filtered_h5_path = config["outputs"]["output_dir"] + row["Pool"] + "/" + row["Selection"] + "/cellbender_feature_bc_matrix_filtered.h5"
                    cellbender_filtered_barcodes_path = config["outputs"]["output_dir"] + row["Pool"] + "/" + row["Selection"] + "/cellbender_feature_bc_matrix_cell_barcodes.csv"
                    if not os.path.exists(cellbender_filtered_h5_path) or not os.path.exists(cellbender_filtered_barcodes_path):
                        logger.info("ERROR: the '{pool}/{selection}/cellbender_feature_bc_matrix_filtered.h5' and/or '{pool}/{selection}/cellbender_feature_bc_matrix_cell_barcodes.csv' does not exist.".format(pool=row["Pool"], selection=row["Selection"]))
                        logger.info("Please check the CellBender outputs and choose the best run (see the docs).")
                        logger.info("Once you are happy with the thresholding, input the correct output directory (e.g. CellBenderDefaultRun1 / CellBenderManualRun1 etc. into the second column of the CellBender_manual_selection.tsv file and restart the snakemake pipeline.")
                        exit()

                    # Save the paths for WG1.
                    samplesheet_data.append([row["Pool"], cell_ranger_bam_path, cellbender_filtered_h5_path, cellbender_filtered_barcodes_path])

                # Save the file_directories.tsv file that points to the files we will use for WG1.
                logger.info("All the CellBender results have PASSED. Saving file directories for WG1.")
                pd.DataFrame(samplesheet_data, columns=["Pool", "Bam", "Counts", "Barcodes"]).to_csv(config["outputs"]["output_dir"] + "wg0_file_directories.tsv", sep="\t", index=False)
        else:
            # Manual selection file did not exist, create it.
            empty_select_df.to_csv(man_select_path, sep="\t", header=True, index=False)

        # Add CellBender jobs. Might be that these were already created but those won't be executed then.
        input_files.extend([config["outputs"]["output_dir"] + pool + "/" + CELLBENDER_SETTINGS[pool]["outdir"]  + "/cellbender_feature_bc_matrix_filtered.h5" for pool in SAMPLES["Pool"]])
        input_files.extend([config["outputs"]["output_dir"] + pool + "/" + CELLBENDER_SETTINGS[pool]["outdir"]  + "/cellbender_feature_bc_matrix_cell_barcodes.csv" for pool in SAMPLES["Pool"]])

        # The summary figures we want for the report. If they exist we should delete them to force snakemake to
        # reproduce them.
        figure_filenames = ["stats", "training_procedure", "cell_probability", "latent_gene_encoding"]
        for fname in figure_filenames:
            if os.path.isfile(fname):
                os.remove(fname)
            input_files.append(config["outputs"]["output_dir"] + "CellBenderCombinedResults/" + fname + ".png")
    return input_files

rule all:
    input:
        files = get_all_input

# Import individual rules
include: "includes/cellranger.smk"
include: "includes/cellbender.smk"
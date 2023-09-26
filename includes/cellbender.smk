#!/usr/bin/env python
import pandas as pd

def calculate_mem_per_thread_gb(wildcards):
    """
    Function to calculate the memory needed for CellBender. Turns out there is a nice linear relationship between CellRanger estimated number
     of cells and the amount of CPU RAM we need for CellBender. Using this function to define the memory reduced the overall memory usage of all jobs.
    """
    # -1 means we want to predict the memory usage.
    if config["cellbender"]["cellbender_memory"] == -1:
        # Find the estimated number of cells from CellRanger.
        metrics_df = pd.read_csv(config["outputs"]["output_dir"] + wildcards.pool + "/CellRanger_count/outs/metrics_summary.csv",sep=",",header=0,index_col=None)
        estimated_number_of_cells = int(metrics_df["Estimated Number of Cells"][0].replace(",",""))

        # Calculate the memory usage.
        return (wildcards.attempt * config["cellbender_extra"]["cellbender_flat_memory"]) + (estimated_number_of_cells * config["cellbender_extra"]["cellbender_scaling_memory"])
    return wildcards.attempt * config["cellbender"]["cellbender_memory"]

rule CellBender:
    input:
        raw_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/raw_feature_bc_matrix.h5",
        metrics_summary = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/metrics_summary.csv"
    output:
        checkpoint = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/ckpt.tar.gz",
        raw_h5 = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix.h5",
        filtered_h5 = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix_filtered.h5",
        filtered_barcodes = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix_cell_barcodes.csv",
        posterior_h5 = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix_posterior.h5",
        metrics = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix_metrics.csv",
        figures = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix.pdf",
        report = report(config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix_report.html", category="CellBender", subcategory="{pool}", caption="../report_captions/CellBender.rst"),
        log = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}/cellbender_feature_bc_matrix.log"
    resources:
        mem_per_thread_gb = calculate_mem_per_thread_gb,
        disk_per_thread_gb = calculate_mem_per_thread_gb,
        threads = config["cellbender"]["cellbender_threads"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + (config["cellbender"]["cellbender_gpu_time"] if config["cellbender"]["cellbender_use_gpu"] else config["cellbender"]["cellbender_cpu_time"])]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "{pool}/CellBender{run_type}Run{run_id}",
        cuda = "--cuda" if config["cellbender"]["cellbender_use_gpu"] else "",
        epochs = lambda wildcards: CELLBENDER_SETTINGS[wildcards.pool]["epochs"],
        learning_rate = lambda wildcards: CELLBENDER_SETTINGS[wildcards.pool]["learning-rate"],
        cpu_threads = "" if config["cellbender"]["cellbender_use_gpu"] else "--cpu-threads " + config["cellbender"]["cellbender_threads"],
    log: config["outputs"]["output_dir"] + "log/CellBender{run_type}Run{run_id}.{pool}.log"
    shell:
        """
        mkdir {params.out} && cd {params.out} || exit
        singularity exec --nv --bind {params.bind} {params.sif} cellbender remove-background \
            --input {input.raw_h5} \
            --output {output.raw_h5} \
            {params.cuda} \
            --checkpoint {output.checkpoint} \
            --epochs {params.epochs} \
            --learning-rate {params.learning_rate} \
            {params.cpu_threads}
        """

# Note that this function actually uses all CellBender runs it can find, not just the one of the current run. This is to allow comparison
# with previous runs.
rule CellBender_report:
    input:
        raw_h5 = expand(config["outputs"]["output_dir"] + "{pool}/{cellbender_outdir}/cellbender_feature_bc_matrix.h5", pool=SAMPLES["Pool"], cellbender_outdir=CELLBENDER_OUTDIRS),
        posterior_h5 = expand(config["outputs"]["output_dir"] + "{pool}/{cellbender_outdir}/cellbender_feature_bc_matrix_posterior.h5", pool=SAMPLES["Pool"], cellbender_outdir=CELLBENDER_OUTDIRS),
        metrics = expand(config["outputs"]["output_dir"] + "{pool}/{cellbender_outdir}/cellbender_feature_bc_matrix_metrics.csv", pool=SAMPLES["Pool"], cellbender_outdir=CELLBENDER_OUTDIRS),
        log = expand(config["outputs"]["output_dir"] + "{pool}/{cellbender_outdir}/cellbender_feature_bc_matrix.log", pool=SAMPLES["Pool"], cellbender_outdir=CELLBENDER_OUTDIRS)
    output:
        stats = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/stats.png", category = "CellBender", subcategory = "CellBender_Summary", caption = "../report_captions/CellBender.rst"),
        training_procedure = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/training_procedure.png", category = "CellBender_Summary", caption = "../report_captions/CellBender.rst"),
        latent_cell_probability = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/cell_probability.png", category = "CellBender_Summary", caption = "../report_captions/CellBender.rst"),
        gene_encoding = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/latent_gene_encoding.png", category = "CellBender", subcategory = "CellBender_Summary", caption = "../report_captions/CellBender.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["cellbender_report_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["cellbender_report_memory"],
        threads = config["cellbender"]["cellbender_report_threads"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["cellbender"]["cellbender_report_time"]]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG0-pipeline-preprocessing/scripts/visualise_cellbender_results.py",
        pools = SAMPLES["Pool"],
        out = config["outputs"]["output_dir"] + "CellBenderCombinedResults"
    log: config["outputs"]["output_dir"] + "log/CellBenderCombinedResults.log"
    shell:
        """
        mkdir {params.out} && cd {params.out} || exit
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
                --pools {params.pools} \
                --raw_h5 {input.raw_h5} \
                --posterior_h5 {input.posterior_h5} \
                --metrics {input.metrics} \
                --log {input.log} \
                --output {params.out}
        """

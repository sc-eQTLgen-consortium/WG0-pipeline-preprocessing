#!/usr/bin/env python
import pandas as pd
import math

def calculate_mem_per_thread_gb(wildcards, attempt):
    """
    Function to calculate the memory needed for CellBender. Turns out there is a nice linear relationship between CellRanger estimated number
     of cells and the amount of RAM we need for CellBender. Using this function to define the memory reduced the overall memory usage of all jobs.
    """
    # -1 means we want to predict the memory usage.
    if config["cellbender"]["cellbender_memory"] == -1:
        # Find the estimated number of cells from CellRanger.
        metrics_df = pd.read_csv(config["outputs"]["output_dir"] + "CellRanger_count/{sample}/outs/metrics_summary.csv".format(sample=wildcards.sample), sep=",", header=0, index_col=None)
        estimated_number_of_cells = int(metrics_df["Estimated Number of Cells"][0].replace(",",""))

        # Calculate the memory usage.
        return math.ceil((attempt * config["cellbender_extra"]["flat_memory"]) + (estimated_number_of_cells * config["cellbender_extra"]["scaling_memory"]))
    return attempt * config["cellbender"]["cellbender_memory"]


# Not implemented options posterior_regularization and estimator, is for "experts" anyway and gave issues.
# Also, need to cd to {params.out} since CellBender likes to make files in the current working directory.
# Note that the python script is there to make the CellBender_settings.json file which we need to allow
# multiple runs without removing previous runs.
rule CellBender:
    input:
        raw_h5 = config["outputs"]["output_dir"] + "CellRanger_count/{sample}/outs/raw_feature_bc_matrix.h5",
        metrics_summary = config["outputs"]["output_dir"] + "CellRanger_count/{sample}/outs/metrics_summary.csv"
    output:
        settings = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/CellBender_settings.json",
        checkpoint = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/ckpt.tar.gz",
        raw_h5 = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix.h5",
        filtered_h5 = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_filtered.h5",
        filtered_barcodes = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_cell_barcodes.csv",
        posterior_h5 = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_posterior.h5",
        metrics = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_metrics.csv",
        figures = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix.pdf",
        report = report(config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_report.html", category="CellBender", subcategory="{sample}", caption="../report_captions/CellBender.rst"),
        log = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix.log",
        done = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/CellBender.done",
    resources:
        mem_per_thread_gb = calculate_mem_per_thread_gb,
        disk_per_thread_gb = calculate_mem_per_thread_gb,
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + (config["cellbender"]["cellbender_gpu_time"] if config["settings"]["use_gpu"] else config["cellbender"]["cellbender_cpu_time"])],
        gpu_sbatch_settings = config["settings"]["gpu_sbatch_settings"] if config["settings"]["use_gpu"] else ""
    threads: config["cellbender"]["cellbender_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG0-pipeline-preprocessing/scripts/cellbender_argparser.py",
        out = config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/",
        cuda = "--cuda" if config["settings"]["use_gpu"] else "",
        expected_cells = lambda wildcards: "--expected-cells " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["expected_cells"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["expected_cells"]) else "",
        total_droplets_included = lambda wildcards: "--total-droplets_-included " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["total_droplets_included"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["total_droplets_included"]) else "",
        force_cell_umi_prior = lambda wildcards: "--force-cell-umi-prior " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["force_cell_umi_prior"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["force_cell_umi_prior"]) else "",
        force_empty_umi_prior = lambda wildcards: "--force-empty-umi-prior " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["force_empty_umi_prior"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["force_empty_umi_prior"]) else "",
        model = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["model"],
        epochs = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["epochs"],
        low_count_threshold = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["low_count_threshold"],
        z_dim = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["z_dim"],
        z_layers = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["z_layers"][1:-1],
        training_fraction = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["training_fraction"],
        empty_drop_training_fraction = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["empty_drop_training_fraction"],
        ignore_features = lambda wildcards: "--ignore-features " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["ignore_features"][1:-1] if CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["ignore_features"] != "[]" else "",
        fpr = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["fpr"][1:-1],
        exclude_feature_types = lambda wildcards: "--exclude-feature-types " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["exclude_feature_types"][1:-1] if CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["exclude_feature_types"] != "[]" else "",
        projected_ambient_count_threshold = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["projected_ambient_count_threshold"],
        learning_rate = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["learning_rate"],
        checkpoint_mins = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["checkpoint_mins"],
        final_elbo_fail_fraction = lambda wildcards: "--final-elbo-fail-fraction " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["final_elbo_fail_fraction"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["final_elbo_fail_fraction"]) else "",
        epoch_elbo_fail_fraction = lambda wildcards: "--epoch-elbo-fail-fraction " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["epoch_elbo_fail_fraction"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["epoch_elbo_fail_fraction"]) else "",
        num_training_tries = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["num_training_tries"],
        learning_rate_retry_mult = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["learning_rate_retry_mult"],
        posterior_batch_size = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["posterior_batch_size"],
        # posterior_regularization = lambda wildcards: "--posterior-regularization " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["posterior_regularization"] if CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["posterior_regularization"] != "null" else "",
        alpha = lambda wildcards: "--alpha " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["alpha"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["alpha"]) else "",
        q = lambda wildcards: "--q " + CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["q"] if not math.isnan(CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["q"]) else "",
        # estimator = lambda wildcards: CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["estimator"],
        estimator_multiple_cpu = lambda wildcards: "--estimator-multiple-cpu" if CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["estimator_multiple_cpu"] else "",
        constant_learning_rate = lambda wildcards: "--constant-learning-rate" if CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["constant_learning_rate"] else "",
        cpu_threads = "--cpu-threads " + str(config["cellbender"]["cellbender_threads"]) if not config["settings"]["use_gpu"] else "",
        debug = lambda wildcards: "--debug" if CELLBENDER_SETTINGS[wildcards.sample][wildcards.run]["debug"] else "",
        nv = "--nv" if config["settings"]["use_gpu"] else ""
    log: config["outputs"]["output_dir"] + "log/CellBender.{sample}.run{run}.log"
    shell:
        """
        cd {params.out} || exit
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --input {input.raw_h5} \
            --output {output.raw_h5} \
            {params.cuda} \
            --checkpoint {output.checkpoint} \
            {params.expected_cells} \
            {params.total_droplets_included} \
            {params.force_cell_umi_prior} \
            {params.force_empty_umi_prior} \
            --model {params.model} \
            --epochs {params.epochs} \
            --low-count-threshold {params.low_count_threshold} \
            --z-dim {params.z_dim} \
            --z-layers {params.z_layers} \
            --training-fraction {params.training_fraction} \
            --empty-drop-training-fraction {params.empty_drop_training_fraction} \
            {params.ignore_features} \
            --fpr {params.fpr} \
            {params.exclude_feature_types} \
            --projected-ambient-count-threshold {params.projected_ambient_count_threshold} \
            --learning-rate {params.learning_rate} \
            --checkpoint-mins {params.checkpoint_mins} \
            {params.final_elbo_fail_fraction} \
            {params.epoch_elbo_fail_fraction} \
            --num-training-tries {params.num_training_tries} \
            --learning-rate-retry-mult {params.learning_rate_retry_mult} \
            --posterior-batch-size {params.posterior_batch_size} \
            {params.alpha} \
            {params.q} \
            {params.estimator_multiple_cpu} \
            {params.constant_learning_rate} \
            {params.cpu_threads} \
            --debug {params.debug}
        
        singularity exec {params.nv} --bind {params.bind} {params.sif} cellbender remove-background \
            --input {input.raw_h5} \
            --output {output.raw_h5} \
            {params.cuda} \
            --checkpoint {output.checkpoint} \
            {params.expected_cells} \
            {params.total_droplets_included} \
            {params.force_cell_umi_prior} \
            {params.force_empty_umi_prior} \
            --model {params.model} \
            --epochs {params.epochs} \
            --low-count-threshold {params.low_count_threshold} \
            --z-dim {params.z_dim} \
            --z-layers {params.z_layers} \
            --training-fraction {params.training_fraction} \
            --empty-drop-training-fraction {params.empty_drop_training_fraction} \
            {params.ignore_features} \
            --fpr {params.fpr} \
            {params.exclude_feature_types} \
            --projected-ambient-count-threshold {params.projected_ambient_count_threshold} \
            --learning-rate {params.learning_rate} \
            --checkpoint-mins {params.checkpoint_mins} \
            {params.final_elbo_fail_fraction} \
            {params.epoch_elbo_fail_fraction} \
            --num-training-tries {params.num_training_tries} \
            --learning-rate-retry-mult {params.learning_rate_retry_mult} \
            --posterior-batch-size {params.posterior_batch_size} \
            {params.alpha} \
            {params.q} \
            {params.estimator_multiple_cpu} \
            {params.constant_learning_rate} \
            {params.cpu_threads} \
            --debug {params.debug}
        
        singularity exec --bind {params.bind} {params.sif} touch {output.done}
        """


# Note that this function actually uses all CellBender runs it can find, not just the one of the current run. This is to allow comparison
# with previous runs.
rule plot_CellBender:
    input:
        input_h5 = [config["outputs"]["output_dir"] + "CellRanger_count/{sample}/outs/raw_feature_bc_matrix.h5".format(sample=sample) for sample in CELLBENDER_SETTINGS.keys()],
        settings = [config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/CellBender_settings.json".format(sample=sample, run=run) for sample in CELLBENDER_SETTINGS.keys() for run in CELLBENDER_SETTINGS[sample].keys()],
        raw_h5 = [config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix.h5".format(sample=sample, run=run) for sample in CELLBENDER_SETTINGS.keys() for run in CELLBENDER_SETTINGS[sample].keys()],
        posterior_h5 = [config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_posterior.h5".format(sample=sample, run=run) for sample in CELLBENDER_SETTINGS.keys() for run in CELLBENDER_SETTINGS[sample].keys()]
    output:
        figure = report(expand(config["outputs"]["output_dir"] + "QC_figures/CellBender_report.{index}.png", index=CELLBENDER_REPORT_INDICES), category="CellBender", caption="../report_captions/CellBender.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["plot_cellbender_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["plot_cellbender_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["cellbender"]["plot_cellbender_time"]],
        gpu_sbatch_settings = ""
    threads: config["cellbender"]["plot_cellbender_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG0-pipeline-preprocessing/scripts/plot_cellbender.py",
        samples = ["{sample}Run{run}".format(sample=sample, run=run) for sample in CELLBENDER_SETTINGS.keys() for run in CELLBENDER_SETTINGS[sample].keys()],
        max_plots_per_page = config["cellbender_extra"]["max_plots_per_page"],
        out = config["outputs"]["output_dir"] + "QC_figures/"
    log: config["outputs"]["output_dir"] + "log/CellBenderCombinedResults.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --samples {params.samples} \
            --input_h5 {input.input_h5} \
            --settings {input.settings} \
            --raw_h5 {input.raw_h5} \
            --posterior_h5 {input.posterior_h5} \
            --max_plots_per_page {params.max_plots_per_page} \
            --out {params.out} 
        """

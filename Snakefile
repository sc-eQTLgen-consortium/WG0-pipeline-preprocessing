#!/usr/bin/env python
import pandas as pd

# Get the samples.
SAMPLES = pd.read_csv(config["inputs"]["samplesheet_filepath"], sep="\t")
SAMPLES = SAMPLES.iloc[:2, :] #TODO: remove
SAMPLES.columns = ["Pool", "ID"]
SAMPLES["Pool"] = SAMPLES["Pool"].astype(str)
SAMPLES["ID"] = SAMPLES["ID"].astype(str)
SAMPLE_DICT = dict(zip(SAMPLES["Pool"], SAMPLES["ID"]))

METHODS = ["CellRanger"]
if config["settings"]["ambient_rna_correction"]:
    logger.info("Running CellBender rules since 'ambient_rna_correction' is set to True.")
    METHODS.append("CellBender")

# Remove trailing /.
if not config["inputs"]["fastqs"].endswith("/"):
    config["inputs"]["fastqs"] += "/"
if not config["outputs"]["output_dir"].endswith("/"):
    config["outputs"]["output_dir"] += "/"

#####################
######## ALL ########
#####################
def get_all_input(wildcards):
    input_files = expand(config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/possorted_genome_bam.bam", pool=SAMPLES["Pool"])

    if "CellRanger" in METHODS:
        # Files we need for WG1.
        input_files.extend(expand(config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/filtered_feature_bc_matrix.h5", pool=SAMPLES["Pool"]))
        input_files.extend(expand(config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", pool=SAMPLES["Pool"]))

    if "CellBender" in METHODS:
        # Files we need for WG1.
        input_files.extend(expand(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_filtered.h5", pool=SAMPLES["Pool"]))
        input_files.extend(expand(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_cell_barcodes.csv", pool=SAMPLES["Pool"]))

        # Files we need for the report.
        input_files.append(config["outputs"]["output_dir"] + "CellBenderCombinedResults/stats.png")
        input_files.append(config["outputs"]["output_dir"] + "CellBenderCombinedResults/training_procedure.png")
        input_files.append(config["outputs"]["output_dir"] + "CellBenderCombinedResults/cell_probability.png")
        input_files.append(config["outputs"]["output_dir"] + "CellBenderCombinedResults/latent_gene_encoding.png")
    return input_files

rule all:
    input:
        files = get_all_input


######################################
############# CellRanger #############
######################################
rule CellRanger:
    input:
        fastqs = config["inputs"]["fastqs"],
        transcriptome = config["refs"]["ref_dir"]
    output:
        raw_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/raw_feature_bc_matrix.h5",
        filtered_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/filtered_feature_bc_matrix.h5",
        filtered_barcodes = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        bam = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/possorted_genome_bam.bam",
        metrics_summary = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/metrics_summary.csv",
        web_summary = report(config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/web_summary.html", category = "CellRanger", subcategory = "{pool}", caption = "../report_captions/CellRanger.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_mem_per_thread_gb"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_disk_per_thread_gb"],
        threads = config["cellranger"]["cellranger_threads"]
    params:
        out = config["outputs"]["output_dir"] + "{pool}/CellRanger",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        id = lambda wildcards: SAMPLE_DICT[wildcards.pool]
    log: config["outputs"]["output_dir"] + "logs/CellRanger.{pool}.log"
    shell:
        """
        cd {params.out} || exit
        singularity exec --bind {params.bind} {params.sif} cellranger count \
            --id {params.id} \
            --fastqs {input.fastqs} \
            --sample {wildcards.pool} \
            --transcriptome {input.transcriptome}
        """

##################################
############ CellBender ##########
##################################
def calculate_mem_per_thread_gb(wildcards):
    metrics_df = pd.read_csv(config["outputs"]["output_dir"] + wildcards.pool + "/CellRanger/outs/metrics_summary.csv", sep=",", header=0, index_col=None)
    estimated_number_of_cells = int(metrics_df["Estimated Number of Cells"][0].replace(",", ""))

    if config["cellbender"]["cellbender_mem_per_thread_gb"] == -1:
        return (wildcards.attempt * config["cellbender_extra"]["cellbender_flat_memory"]) + (estimated_number_of_cells * config["cellbender_extra"]["cellbender_scaling_memory"])
    return wildcards.attempt * config["cellbender_extra"]["cellbender_mem_per_thread_gb"]

rule CellBender:
    input:
        raw_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/raw_feature_bc_matrix.h5",
        metrics_summary = config["outputs"]["output_dir"] + "{pool}/CellRanger/outs/metrics_summary.csv"
    output:
        checkpoint = config["outputs"]["output_dir"] + "{pool}/CellBender/ckpt.tar.gz",
        raw_h5 = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix.h5",
        filtered_h5 = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_filtered.h5",
        filtered_barcodes = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_cell_barcodes.csv",
        posterior_h5 = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_posterior.h5",
        metrics = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_metrics.csv",
        figures = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix.pdf",
        report = report(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_report.html", category="CellBender", subcategory="{pool}", caption="../report_captions/CellBender.rst"),
        log = config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix.log"
    resources:
        mem_per_thread_gb = calculate_mem_per_thread_gb,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["cellbender_disk_per_thread_gb"],
        threads = config["cellbender"]["cellbender_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        out = config["outputs"]["output_dir"] + "{pool}/CellBender",
        cuda = "--cuda" if config["cellbender"]["cellbender_use_gpu"] else "",
        epochs = config["cellbender_extra"]["cellbender_epochs"],
        learning_rate = config["cellbender_extra"]["cellbender_learning_rate"],
        cpu_threads = "" if config["cellbender"]["cellbender_use_gpu"] else "--cpu-threads " + config["cellbender"]["cellbender_threads"],
    log: config["outputs"]["output_dir"] + "logs/CellBender.{pool}.log"
    shell:
        """
        cd {params.out} || exit
        
        singularity exec --nv --bind {params.bind} {params.sif} cellbender remove-background \
            --input {input.raw_h5} \
            --output {output.raw_h5} \
            {params.cuda} \
            --checkpoint {output.checkpoint} \
            --epochs {params.epochs} \
            --learning-rate {params.learning_rate} \
            {params.cpu_threads}
        """

rule CellBender_report:
    input:
        raw_h5 = expand(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix.h5", pool=SAMPLES["Pool"]),
        posterior_h5 = expand(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_posterior.h5", pool=SAMPLES["Pool"]),
        metrics = expand(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix_metrics.csv", pool=SAMPLES["Pool"]),
        log = expand(config["outputs"]["output_dir"] + "{pool}/CellBender/cellbender_feature_bc_matrix.log", pool=SAMPLES["Pool"])
    output:
        stats = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/stats.png", category = "CellBender", subcategory = "CellBender_Summary", caption = "../report_captions/CellBender.rst"),
        training_procedure = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/training_procedure.png", category = "CellBender_Summary", caption = "../report_captions/CellBender.rst"),
        latent_cell_probability = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/cell_probability.png", category = "CellBender_Summary", caption = "../report_captions/CellBender.rst"),
        gene_encoding = report(config["outputs"]["output_dir"] + "CellBenderCombinedResults/latent_gene_encoding.png", category = "CellBender", subcategory = "CellBender_Summary", caption = "../report_captions/CellBender.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["cellbender_report_mem_per_thread_gb"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellbender"]["cellbender_report_disk_per_thread_gb"],
        threads = config["cellbender"]["cellbender_report_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = "/opt/WG0-pipeline-preprocessing/scripts/visualise_cellbender_results.py",
        pools = SAMPLES["Pool"],
        out = config["outputs"]["output_dir"] + "/CellBenderCombinedResults"
    log: config["outputs"]["output_dir"] + "logs/CellBenderCombinedResults.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
                --pools {params.pools} \
                --raw_h5 {input.raw_h5} \
                --posterior_h5 {input.posterior_h5} \
                --metrics {input.metrics} \
                --log {input.log} \
                --output {params.out}
        """

#!/usr/bin/env python


# Create the input file for WG1.
rule combine_results:
    input:
        counts = [config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_filtered.h5".format(sample=sample, run=CELLBENDER_SELECTION[sample]) for sample in SAMPLES] if config["settings"]["ambient_rna_correction"] else expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/filtered_feature_bc_matrix.h5", sample=SAMPLES),
        barcodes = [config["outputs"]["output_dir"] + "CellBender/{sample}Run{run}/cellbender_feature_bc_matrix_cell_barcodes.csv".format(sample=sample, run=CELLBENDER_SELECTION[sample]) for sample in SAMPLES] if config["settings"]["ambient_rna_correction"] else expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample=SAMPLES),
        bam = expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/possorted_genome_bam.bam", sample=SAMPLES),
    output:
        summary = config["outputs"]["output_dir"] + "Combine_Results/wg0_file_directories.tsv"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["combine_results"]["combine_results_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["combine_results"]["combine_results_time"]],
        gpu_sbatch_settings = ""
    threads: config["combine_results"]["combine_results_threads"]
    params:
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        script = config["inputs"]["repo_dir"] + "scripts/combine_results.py",
        samples = " ".join(SAMPLES),
        out = config["outputs"]["output_dir"] + "Combine_Results/",
    log: config["outputs"]["output_dir"] + "log/combine_results.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --pools {params.samples} \
            --counts {input.counts} \
            --barcodes {input.barcodes} \
            --bam {input.bam} \
            --out {params.out}
        """
#!/usr/bin/env python


# Create the input file for WG1.
rule combine_results:
    input:
        samplesheet = config["inputs"]["samplesheet_path"]
    output:
        file_directories = config["outputs"]["output_dir"] + "Combine_Results/wg0_file_directories.tsv"
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
        cellbender_selection = "--cellbender_selection " + config["outputs"]["output_dir"] + "manual_selection/CellBender_manual_selection.tsv" if config["settings"]["ambient_rna_correction"] else "",
        basedir = config["outputs"]["output_dir"]
    log: config["outputs"]["output_dir"] + "log/combine_results.log"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} python {params.script} \
            --samplesheet {input.samplesheet} \
            {params.cellbender_selection} \
            --basedir {params.basedir} \
            --outfile {output.file_directories} > {log} 2>&1
        """
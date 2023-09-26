#!/usr/bin/env python


rule CellRanger_count:
    input:
        fastqs = config["inputs"]["fastqs"],
        transcriptome = config["refs"]["ref_dir"]
    output:
        raw_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/raw_feature_bc_matrix.h5",
        filtered_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/filtered_feature_bc_matrix.h5",
        filtered_barcodes = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        molecule_info_h5 = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/molecule_info.h5",
        bam = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/possorted_genome_bam.bam",
        metrics_summary = config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/metrics_summary.csv",
        web_summary = report(config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/web_summary.html", category = "CellRanger", subcategory = "{pool}", caption = "../report_captions/CellRanger.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_count_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_count_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["cellranger"]["cellranger_count_time"]]
    threads: config["cellranger"]["cellranger_count_threads"]
    params:
        out = config["outputs"]["output_dir"] + "{pool}/CellRanger_count",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        id = lambda wildcards: SAMPLE_DICT[wildcards.pool]
    log: config["outputs"]["output_dir"] + "log/CellRanger_count.{pool}.log"
    shell:
        """
        mkdir {params.out} && cd {params.out} || exit
        singularity exec --bind {params.bind} {params.sif} cellranger count \
            --id {params.id} \
            --fastqs {input.fastqs} \
            --sample {wildcards.pool} \
            --transcriptome {input.transcriptome}
        """

# This rules only function is to create the html file. If it turns out this html is not useful we can also delete this rule.
rule CellRanger_aggr:
    input:
        csv = aggregation_csv_path,
        molecule_info_h5 = expand(config["outputs"]["output_dir"] + "{pool}/CellRanger_count/outs/molecule_info.h5", pool=SAMPLES["Pool"]),
    output:
        summary = config["outputs"]["output_dir"] + "/CellRanger_aggr/outs/summary.json",
        web_summary = report(config["outputs"]["output_dir"] + "CellRanger_aggr/outs/web_summary.html", category = "CellRanger", caption = "../report_captions/CellRanger.rst")
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_aggr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_aggr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][attempt + config["cellranger"]["cellranger_aggr_time"]]
    threads: config["cellranger"]["cellranger_aggr_threads"]
    params:
        out = config["outputs"]["output_dir"] + "CellRanger_aggr",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        id = "Aggregate"
    log: config["outputs"]["output_dir"] + "log/CellRanger_aggr.log"
    shell:
        """
        mkdir {params.out} && cd {params.out} || exit
        singularity exec --bind {params.bind} {params.sif} cellranger aggr \
            --id {params.id} \
            --csv {input.csv} \
            --localcores {resources.threads} \
            --localmem {resources.mem_per_thread_gb}
        """
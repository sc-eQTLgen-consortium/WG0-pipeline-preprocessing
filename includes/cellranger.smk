#!/usr/bin/env python


# Need to remove the sample directory created by Snakemake to prevent 'RuntimeError: {params.out} is not a pipestance directory'
# Also, not defining _<file> files os output so those won't be deleted if the rule crashes but I am not sure if they are all always created
# so I will jut ignore them.
rule CellRanger_count:
    input:
        fastqs = config["inputs"]["fastq_dir"],
        transcriptome = config["refs"]["ref_dir"] + config["refs_extra"]["transcriptome_path"]
    output:
        mri = config["outputs"]["output_dir"] + "CellRanger/{sample}/{sample}.mri.tgz",
        clustering_graph = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/clustering/gene_expression_graphclust/clusters.csv",
        lustering_kmeans = expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/clustering/gene_expression_kmeans_{cluster}_clusters/clusters.csv", cluster=range(2,11), allow_missing=True),
        diffexp_graph = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/diffexp/gene_expression_graphclust/differential_expression.csv",
        diffexp_kmeans = expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/diffexp/gene_expression_kmeans_{cluster}_clusters/differential_expression.csv", cluster=range(2,11), allow_missing=True),
        pca = expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/pca/gene_expression_10_components/{pca_file}.csv", pca_file=["components", "dispersion", "features_selected", "projection", "variance"], allow_missing=True),
        tsne = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/tsne/gene_expression_2_components/projection.csv",
        umap = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/analysis/umap/gene_expression_2_components/projection.csv",
        cloupe = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/cloupe.cloupe",
        filtered_barcodes = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        filtered_features = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/filtered_feature_bc_matrix/features.tsv.gz",
        filtered_matrix = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
        filtered_h5 = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/filtered_feature_bc_matrix.h5",
        metrics_summary = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/metrics_summary.csv",
        molecule_info_h5 = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/molecule_info.h5",
        bam = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/possorted_genome_bam.bam",
        bam_index = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/possorted_genome_bam.bam.bai",
        raw_barcodes = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz",
        raw_features = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/raw_feature_bc_matrix/features.tsv.gz",
        raw_matrix = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz",
        raw_h5 = config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/raw_feature_bc_matrix.h5",
        web_summary = report(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/web_summary.html", category="CellRanger", subcategory="{sample}", caption=config["inputs"]["repo_dir"] + "report_captions/CellRanger.rst"),
        sc_rna_counter_cs = directory(config["outputs"]["output_dir"] + "CellRanger/{sample}/SC_RNA_COUNTER_CS/"),
    resources:
        localmem = lambda wildcards, attempt: (attempt * config["cellranger"]["cellranger_count_memory"] * config["cellranger"]["cellranger_count_threads"] - config["settings_extra"]["memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_count_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_count_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["cellranger"]["cellranger_count_time"]],
        gpu_sbatch_settings = ""
    threads: config["cellranger"]["cellranger_count_threads"]
    params:
        out = config["outputs"]["output_dir"] + "CellRanger",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        sample = lambda wildcards: FASTQ_DICT[wildcards.sample]
    log: config["outputs"]["output_dir"] + "log/CellRanger_count.{sample}.log"
    shell:
        """
        cd {params.out} || exit 
        rm -r {wildcards.sample}
        singularity exec --bind {params.bind} {params.sif} cellranger count \
            --id {wildcards.sample} \
            --fastqs {input.fastqs} \
            --sample {params.sample} \
            --transcriptome {input.transcriptome} \
            --localcores {threads} \
            --localmem {resources.localmem}
        """

# This rules only function is to create a merged html file. If it turns out this html is not useful we can also delete this rule.
rule CellRanger_aggr:
    input:
        csv = config["outputs"]["output_dir"] + "CellRanger/aggregation.csv",
        molecule_info_h5 = expand(config["outputs"]["output_dir"] + "CellRanger/{sample}/outs/molecule_info.h5", sample=SAMPLES),
    output:
        mri = config["outputs"]["output_dir"] + "CellRanger/Aggregate/Aggregate.mri.tgz",
        aggregation = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/aggregation.csv",
        clustering_graph = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/clustering/gene_expression_graphclust/clusters.csv",
        lustering_kmeans = expand(config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/clustering/gene_expression_kmeans_{cluster}_clusters/clusters.csv", cluster=range(2,11)),
        diffexp_graph = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/diffexp/gene_expression_graphclust/differential_expression.csv",
        diffexp_kmeans = expand(config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/diffexp/gene_expression_kmeans_{cluster}_clusters/differential_expression.csv", cluster=range(2,11)),
        pca = expand(config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/pca/gene_expression_10_components/{pca_file}.csv", pca_file=["components", "dispersion", "features_selected", "projection", "variance"]),
        tsne = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/tsne/gene_expression_2_components/projection.csv",
        umap = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/analysis/umap/gene_expression_2_components/projection.csv",
        cloupe = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/cloupe.cloupe",
        filtered_barcodes = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz",
        filtered_features = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/filtered_feature_bc_matrix/features.tsv.gz",
        filtered_matrix = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/filtered_feature_bc_matrix/matrix.mtx.gz",
        filtered_h5 = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/filtered_feature_bc_matrix.h5",
        summary = config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/count/summary.json",
        web_summary = report(config["outputs"]["output_dir"] + "CellRanger/Aggregate/outs/web_summary.html", category="CellRanger", subcategory="Aggregate",  caption=config["inputs"]["repo_dir"] + "report_captions/CellRanger.rst"),
        sc_rna_aggregator_cs = directory(config["outputs"]["output_dir"] + "CellRanger/Aggregate/SC_RNA_AGGREGATOR_CS/")
    resources:
        localmem = lambda wildcards, attempt: (attempt * config["cellranger"]["cellranger_aggr_memory"] * config["cellranger"]["cellranger_aggr_threads"] - config["settings_extra"]["memory_buffer"]),
        mem_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_aggr_memory"],
        disk_per_thread_gb = lambda wildcards, attempt: attempt * config["cellranger"]["cellranger_aggr_memory"],
        time = lambda wildcards, attempt: config["cluster_time"][(attempt - 1) + config["cellranger"]["cellranger_aggr_time"]],
        gpu_sbatch_settings = ""
    threads: config["cellranger"]["cellranger_aggr_threads"]
    params:
        out = config["outputs"]["output_dir"] + "CellRanger/",
        bind = config["inputs"]["bind_path"],
        sif = config["inputs"]["singularity_image"],
        id = "Aggregate"
    log: config["outputs"]["output_dir"] + "log/CellRanger_aggr.log"
    shell:
        """
        cd {params.out} || exit 
        rm -r Aggregate
        singularity exec --bind {params.bind} {params.sif} cellranger aggr \
            --id {params.id} \
            --csv {input.csv} \
            --localcores {threads} \
            --localmem {resources.localmem}
        """

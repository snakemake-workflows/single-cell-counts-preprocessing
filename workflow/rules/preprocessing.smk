rule convert_raw_counts_to_h5ad:
    input:
        matrix=lookup(
            within=samples,
            query="sample_id == '{wildcards.sample}'",
            cols="raw_counts_path",
        ),
    output:
        h5ad="<results>/raw_counts/{sample}/{sample}.raw_counts.h5ad",
    log:
        "<logs>/raw_counts/{sample}/{sample}.convert_raw_counts_to_h5ad.log",
    conda:
        "../envs/scanpy.yaml"
    params:
        input_format=lookup(
            within=samples, query="sample_id == '{wildcards.sample}'", cols="format"
        ),
    script:
        "../scripts/convert_raw_counts_to_h5ad.py"


rule filter_low_quality_barcodes:
    input:
        h5ad="<results>/raw_counts/{sample}/{sample}.raw_counts.h5ad",
    output:
        h5ad="<results>/filtered_counts/{sample}/{sample}.filtered_counts.h5ad",
        raw_total_counts_plot="<results>/raw_counts/{sample}/{sample}.raw_total_counts.html",
        raw_total_counts_plot_png="<results>/raw_counts/{sample}/{sample}.raw_total_counts.png",
        pct_counts_mt_plot="<results>/raw_counts/{sample}/{sample}.pct_counts_mt.html",
        pct_counts_mt_plot_png="<results>/raw_counts/{sample}/{sample}.pct_counts_mt.png",
        log1p_total_counts="<results>/raw_counts/{sample}/{sample}.log1p_total_counts.html",
        log1p_total_counts_png="<results>/raw_counts/{sample}/{sample}.log1p_total_counts.png",
        log1p_n_genes_by_counts="<results>/raw_counts/{sample}/{sample}.log1p_n_genes_by_counts.html",
        log1p_n_genes_by_counts_png="<results>/raw_counts/{sample}/{sample}.log1p_n_genes_by_counts.png",
        pct_counts_in_top_20_genes="<results>/raw_counts/{sample}/{sample}.pct_counts_in_top_20_genes.html",
        pct_counts_in_top_20_genes_png="<results>/raw_counts/{sample}/{sample}.pct_counts_in_top_20_genes.png",
    log:
        "<logs>/filtered_counts/{sample}.filtered_counts.log",
    conda:
        "../envs/scanpy.yaml"
    params:
        min_umis_per_barcode=lookup(
            within=config, dpath="filtering/min_umis_per_barcode"
        ),
        counts_mads=lookup(within=config, dpath="filtering/counts_mads"),
        mt_percent=lookup(within=config, dpath="filtering/mt_percent"),
        mt_percent_mads=lookup(within=config, dpath="filtering/mt_percent_mads"),
    script:
        "../scripts/filter_low_quality_barcodes.py"

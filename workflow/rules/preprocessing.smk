rule convert_raw_counts_to_zarr:
    input:
        matrix=lookup(
            within=samples,
            query="sample_id == '{wildcards.sample}'",
            cols="raw_counts_path",
        ),
    output:
        zarr=directory("<results>/raw_counts/{sample}/{sample}.raw_counts.zarr"),
        raw_total_counts_plot="<results>/raw_counts/{sample}/{sample}.raw_total_counts.html",
        pct_counts_mt_plot="<results>/raw_counts/{sample}/{sample}.pct_counts_mt.html",
    log:
        "<logs>/raw_counts/{sample}/{sample}.raw_counts.to_zarr.log",
    conda:
        "../envs/scanpy.yaml"
    params:
        input_format=lookup(
            within=samples, query="sample_id == '{wildcards.sample}'", cols="format"
        ),
    script:
        "../scripts/convert_raw_counts_to_zarr.py"


rule filter_low_quality_cells:
    input:
        zarr="<results>/raw_counts/{sample}/{sample}.raw_counts.zarr",
    output:
        filtered="<results>/filtered_counts/{sample}/{sample}.filtered_counts.zarr",
        filtered_plot="results/quality_control/{sample}.filtered.html",
    log:
        "logs/filtered/{sample}.quality_control.log",
    conda:
        "../envs/scanpy.yaml"
    params:
        counts_mads=lookup(within=config, dpath="filtering/counts_mads"),
        mt_percent=lookup(within=config, dpath="filtering/mt_percent"),
        mt_percent_mads=lookup(within=config, dpath="filtering/mt_percent_mads"),
    script:
        "../scripts/quality_control.py"

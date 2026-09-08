import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

import numpy as np
import anndata as ad
import scanpy as sc
import altair as alt

from scipy.stats import median_abs_deviation
from math import ceil

adata = ad.read_zarr(snakemake.input["zarr"])

counts_mads = snakemake.params["counts_mads"]
mt_percent = snakemake.params["mt_percent"]
mt_percent_mads = snakemake.params["mt_percent_mads"]

# do very basic filtering of raw data, to get to a reasonable data size
sc.pp.filter_cells(adata, min_counts=30)

# mitochondrial genes
# they can start with different prefixes, depending on the
# species, with human mitochondrial genes usually starting with `MT-` and
# mouse mitochondrial genes usually staritng with `mt-`
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))

# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    inplace=True,
    percent_top=[20],
    log1p=True,
)

def median_abs_deviation_upper_bound(adata, metric, nmads):
    met = adata.obs[metric]
    return np.median(met) + nmads * median_abs_deviation(met)

def median_abs_deviation_lower_bound(adata, metric, nmads):
    met = adata.obs[metric]
    return np.median(met) - nmads * median_abs_deviation(met)

def plot_value_per_barcode(
    adata: AnnData,
    obs_column: str,
    per_barcode_value_description: str | list[str],
    plot_output_name_in_rule: str,
    nmads: int | None = None,
    threshold: float | None = None,
):
    df = adata.obs[[obs_column]]
    df["all_barcodes"] = "all barcodes"
    df.index.name = "barcode"
    df.reset_index(inplace=True)
    ymax = ceil(max(df[obs_column]))
    per_barcode_base_chart = alt.Chart(df)
    p_dot = (
        per_barcode_base_chart.mark_circle(filled=False)
        .encode(
            alt.Y(f"{obs_column}:Q")
            .scale(domain=(0, ymax))
            .title(per_barcode_value_description)
            .axis(labelLimit=200),
            alt.X("all_barcodes:N").axis(labels=False),
            xOffset="jitter:Q",
            tooltip=["barcode", obs_column],
        )
        .transform_calculate(jitter="random()")
        .properties(width=alt.Step(200))
    )
    p_bar = per_barcode_base_chart.mark_bar().encode(
        alt.Y(f"{obs_column}:Q")
        .bin(maxbins=100)
        .scale(domain=(0, ymax))
        .title(per_barcode_value_description),
        alt.X("count()").title("number of barcodes per bin"),
    )
    if nmads:
        p_upper_threshold = per_barcode_base_chart.mark_rule(color='red').encode(
            y=alt.datum(median_abs_deviation_upper_bound(adata, obs_column, nmads)),
        )
        p_lower_threshold = per_barcode_base_chart.mark_rule(color='red').encode(
            y=alt.datum(median_abs_deviation_lower_bound(adata, obs_column, nmads)),
        )
        if threshold:
            p_threshold = per_barcode_base_chart.mark_rule(color='blue').encode(
                y=alt.datum(threshold)
            )
            p = p_dot + p_upper_threshold + p_lower_threshold + p_threshold | p_bar + p_upper_threshold + p_lower_threshold + p_threshold
        else:
            p = p_dot + p_upper_threshold + p_lower_threshold | p_bar + p_upper_threshold + p_lower_threshold
    else:
        p = p_dot | p_bar
    p.save(snakemake.output[plot_output_name_in_rule])
    p.save(snakemake.output[f"{plot_output_name_in_rule}_png"])

plot_value_per_barcode(
    adata=adata,
    obs_column="total_counts",
    per_barcode_value_description="total counts per barcode",
    plot_output_name_in_rule="raw_total_counts_plot",
)

plot_value_per_barcode(
    adata=adata,
    obs_column="pct_counts_mt",
    per_barcode_value_description=[
        "percentage of mitochondrial counts per barcode",
        f"(red: {mt_percent_mads} median absolute deviations thresholds,",
        f"blue: {mt_percent} % absolute maximum threshold)",
    ],
    plot_output_name_in_rule="pct_counts_mt_plot",
    nmads=mt_percent_mads,
    threshold=mt_percent,
)

plot_value_per_barcode(
    adata=adata,
    obs_column="log1p_total_counts",
    per_barcode_value_description=[
        "log1p of total counts per barcode",
        f"(red: {counts_mads} median absolute deviations thresholds)",
    ],
    plot_output_name_in_rule="log1p_total_counts",
    nmads=counts_mads,
)

plot_value_per_barcode(
    adata=adata,
    obs_column="log1p_n_genes_by_counts",
    per_barcode_value_description=[
        "log1p of the number of genes with counts per barcode",
        f"({counts_mads} median absolute deviations thresholds)",
    ],
    plot_output_name_in_rule="log1p_n_genes_by_counts",
    nmads=counts_mads,
)

plot_value_per_barcode(
    adata=adata,
    obs_column="pct_counts_in_top_20_genes",
    per_barcode_value_description=[
        "percentage of counts from the top 20 genes per counts",
        f"({counts_mads} median absolute deviations thresholds)",
    ],
    plot_output_name_in_rule="pct_counts_in_top_20_genes",
    nmads=counts_mads,
)


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < median_abs_deviation_lower_bound(adata, metric, nmads)) | (
        median_abs_deviation_upper_bound(adata, metric, nmads) < M
    )
    return outlier



adata.obs["counts_outlier"] = (
    is_outlier(adata, "log1p_total_counts", counts_mads)
    | is_outlier(adata, "log1p_n_genes_by_counts", counts_mads)
    | is_outlier(adata, "pct_counts_in_top_20_genes", counts_mads)
)

print(
    f"Filtering based on {counts_mads} median absolute standard deviations\n"
    "of log1p_total_counts, log1p_n_genes_by_counts and\n"
    "pct_counts_in_top_20_genes will:\n"
    f"remove: {adata.obs.counts_outlier.value_counts()[True]} barcodes\n"
    f"keep:   {adata.obs.counts_outlier.value_counts()[False]} barcodes.\n",
    file=sys.stderr
)


adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_percent_mads) | (
    adata.obs["pct_counts_mt"] > mt_percent
)

print(
    f"Filtering based on {mt_percent_mads} median absolute standard deviations\n"
    f"of pct_counts_mt and a maximum of {mt_percent} % mitochondrial counts\n"
    "will:\n"
    f"remove: {adata.obs.mt_outlier.value_counts()[True]} barcodes\n"
    f"keep:   {adata.obs.mt_outlier.value_counts()[False]} barcodes.\n",
    file=sys.stderr
)

print(f"Total number of barcodes before filtering of low quality barcodes: {adata.n_obs}\n", file=sys.stderr)
adata = adata[(~adata.obs.counts_outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Total number of barcodes after filtering of low quality barcodes: {adata.n_obs}\n", file=sys.stderr)

adata.write_zarr(snakemake.output["zarr"])

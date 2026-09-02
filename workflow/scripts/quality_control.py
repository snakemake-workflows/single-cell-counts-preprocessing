import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

import warnings
import numpy as np
import scanpy as sc
import altair as alt

from scipy.stats import median_abs_deviation
from math import ceil

# Set figure parameters for clean, minimal plots
sc.settings.set_figure_params(dpi=150, facecolor="white", frameon=False)

input_format = snakemake.params.get("input_format", "")
matrix_file = snakemake.input["matrix"]

match input_format:
    case "10x_h5":
        adata = sc.read_10x_h5(matrix_file)
    case "10x_mtx":
        adata = sc.read_10x_mtx(path=f"{matrix_file}/filtered_feature_bc_matrix/")
    case "visium":
        adata = sc.read_visium(matrix_file)
    case "h5ad":
        adata = sc.read_h5ad(matrix_file)
    case "csv":
        adata = sc.read_csv(matrix_file)
    case "excel":
        adata = sc.read_excel(matrix_file)
    case "hdf":
        adata = sc.read_hdf(matrix_file)
    case "loom":
        adata = sc.read_loom(matrix_file)
    case "mtx":
        adata = sc.read_mtx(matrix_file)
    case "text":
        adata = sc.read_text(matrix_file)
    case "umi_tools":
        adata = sc.read_umi_tools(matrix_file)
    case _:
        warnings.warn(
            f"Unknown input_format specified in `config.yaml`: {input_format}.\n"
            "For known input formats, see the suffixes of scanpy read functions:\n"
            "https://scanpy.scverse.org/en/stable/api/io.html\n"
            "For example, possible options are `10x_mtx`, `visium`, `csv`, or `loom`."
            f"Will try reading file with scanpy.read({matrix_file}), see log file:\n"
            f"{snakemake.log[0]}"
        )
        adata = sc.read(matrix_file)

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



def plot_value_per_barcode(
    adata: AnnData,
    obs_column: str,
    per_barcode_value_description: str,
    plot_output_name_in_rule: str,
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
            .title(f"{per_barcode_value_description} per barcode"),
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
        .title(f"{per_barcode_value_description} per barcode"),
        alt.X("count()").title("number of barcodes per bin"),
    )
    p = p_dot | p_bar
    p.save(snakemake.output[plot_output_name_in_rule])

plot_value_per_barcode(
    adata=adata,
    obs_column="total_counts",
    per_barcode_value_description="total counts",
    plot_output_name_in_rule="raw_total_counts_plot",
)

plot_value_per_barcode(
    adata=adata,
    obs_column="pct_counts_mt",
    per_barcode_value_description="percentage of mitochondrial counts",
    plot_output_name_in_rule="pct_counts_mt_plot",
)

breakpoint()

p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# TODO: write plots to some output file(s)


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


counts_mads = snakemake.params["counts_mads"]

adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", counts_mads)
    | is_outlier(adata, "log1p_n_genes_by_counts", counts_mads)
    | is_outlier(adata, "pct_counts_in_top_20_genes", counts_mads)
)

adata.obs.outlier.value_counts()

mt_percent = snakemake.params["mt_percent"]
mt_percent_mads = snakemake.params["mt_percent_mads"]

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", mt_percent_mads) | (
    adata.obs["pct_counts_mt"] > mt_percent
)

adata.obs.mt_outlier.value_counts()

print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

p4 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

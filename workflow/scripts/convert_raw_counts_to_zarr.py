import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

import warnings
from pathlib import Path
from math import ceil

import scanpy as sc
import altair as alt

input_format = snakemake.params.get("input_format", "")[0]
matrix_file = Path(snakemake.input["matrix"][0])

match input_format:
    case "10x_h5":
        adata = sc.read_10x_h5(matrix_file)
    case "10x_mtx":
        adata = sc.read_10x_mtx(matrix_file)
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
    case "zarr":
        adata = sc.read_zarr(matrix_file)
    case _:
        warnings.warn(
            f"Unknown input_format specified in `sample_sheet`: {input_format}.\n"
            "For known input formats, see the suffixes of scanpy read functions:\n"
            "https://scanpy.scverse.org/en/stable/api/io.html\n"
            "For example, possible options are `10x_mtx`, `visium`, `csv`, or `loom`."
            f"Will try reading file with scanpy.read({matrix_file}), check the log file:\n"
            f"{snakemake.log[0]}"
        )
        adata = sc.read(matrix_file)

adata.var_names_make_unique()

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
#            tooltip=["barcode", obs_column],
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


adata.write_zarr(snakemake.output["zarr"])
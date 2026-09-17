import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

import warnings

import scanpy as sc

input_format = snakemake.params.get("input_format", "")[0]
matrix_file = snakemake.input["matrix"][0]

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

adata.write_h5ad(snakemake.output["h5ad"], compression="gzip")
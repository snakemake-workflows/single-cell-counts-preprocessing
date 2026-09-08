## Workflow overview

This workflow is a best-practice workflow for preprocessing counts from single cell RNA sequencing data.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and follows the [`Preprocessing and visualization` section of the Single Cell Best Practices](https://www.sc-best-practices.org/preprocessing-visualization).

It consists of the following steps:

1. Convert input data to `zarr` format.
2. Filter low-quality barcodes with `scanpy`.
3. Correct for ambient RNA contamination with `SoupX`.
4. Detect doublets with `scDblFinder`.
5. Normalize counts with `scanpy`.

## Running the workflow

To configure the workflow run, go through the provided `config/config.yaml` entry by entry and adjust them where necessary.
After an initial run, we recommend going through the quality control plots mentioned there and double-checking that the provided threshold values for filtering make sense.

### Input data

To specify the input, provide a `config/sample_sheet.tsv` file with the following layout:

| sample_id           | raw_counts_path                   | format   |
| ------------------- | --------------------------------- | -------- |
| cellranger_1        | ../path/to/raw_feature_bc_matrix/ | 10x_mtx  |
| kallisto_bustools_1 | ../path/to/adata.h5ad             | h5ad     |
| alevin_fry_1        | ../path/to/quants_mat.mtx         | mtx      |

Here, the columns are:

* `sample_id`: An arbitrary string identifier of a a sample (or dataset).
* `raw_counts_path`: The path to a file or folder with the raw counts as determined by another tool, for example [`CellRanger`](https://www.10xgenomics.com/support/software/cell-ranger/latest), [`kallisto bustools`](https://kallisto.readthedocs.io/en/latest/sc/pseudoalignment.html) or [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/index.html). We really recommend using the raw counts here (and not any pre-filtered counts), as they are instrumental in the correction for ambient RNA contamination. The workflow filters low-quality barcodes itself, and lets you transparently configure and check the filter thresholds (see the comments in `config/config.yaml`).
* `format`: Format of the input data given in column `raw_counts_path`. Choose a format that scanpy can read, so any suffix in one of the `scanpy.read_` functions. See https://scanpy.scverse.org/en/stable/api/io.html or double-check the code at: https://github.com/scverse/scanpy/blob/a656a33b080a5c1f64b01e841daad76f35f5ec5f/src/scanpy/io/_read.py#L44-L61"
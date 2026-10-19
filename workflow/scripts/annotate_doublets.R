log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

rlang::global_entrace()

library("anndataR")
library("scDblFinder")
library("SingleCellExperiment")

sce <- read_h5ad(
    snakemake@input[["soup_corrected"]],
    as = "SingleCellExperiment"
)

sce_doublets_annotated <- scDblFinder(sce)

write_h5ad(
    object = sce_doublets_annotated,
    path = file.path(snakemake@output["doublets_annotated"]),
    compression = "gzip"
)

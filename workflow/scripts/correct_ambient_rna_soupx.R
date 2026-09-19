log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

rlang::global_entrace()

library("anndataR")
library("SoupX")
library("Seurat")

raw_adata <- read_h5ad(snakemake@input[["raw_h5ad"]])
filtered_adata <- read_h5ad(snakemake@input[["filtered_h5ad"]])

sc <- SoupChannel(
    tod = t(raw_adata$X), # this t() is from Seurat
    toc = t(filtered_adata$X) # this t() is from Seurat
)

sc <- setClusters(sc, unlist(filtered_adata$obs["soupx_groups"]))

pdf_filename <- as.character(snakemake@output["contamination_estimation"])
print(typeof(pdf_filename))
pdf(pdf_filename, width = 8, height = 8)
sc <- autoEstCont(sc, doPlot = FALSE)
dev.off()

out <- adjustCounts(sc, roundToInt = TRUE)

filtered_adata$layers <- list(
    "counts" = filtered_adata$X,
    "soupx_counts" = t(out)
)

write_h5ad(
    object = filtered_adata,
    path = file.path(snakemake@output["soup_corrected"]),
    compression = "gzip"
)

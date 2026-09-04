# import basic packages
import pandas as pd
from snakemake.utils import validate

# read sample sheet
samples = (
    pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
    .set_index("sample_id", drop=False)
    .sort_index()
)


# validate sample sheet and config file
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")


def get_final_results(wildcards):
    final_results = []

    final_results.extend(
        expand(
            "<results>/raw_counts/{sample}/{sample}.raw_counts.zarr",
            sample=samples.sample_id,
        )
    )

    return final_results

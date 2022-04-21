import pandas as pd


def merge(paths_to_datasets: str, path_to_output: str):
    (
        pd
        .concat([pd.read_feather(path_to_dataset) for path_to_dataset in paths_to_datasets], axis=0)
        .reset_index(drop=True)
        .to_feather(path_to_output)
    )


if __name__ == "__main__":
    merge(
        paths_to_datasets=snakemake.input.datasets,
        path_to_output=snakemake.output[0]
    )

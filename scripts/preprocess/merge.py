import pandas as pd


def merge(paths_to_datasets, path_to_output):
    (
        pd
        .concat([pd.read_csv(path_to_dataset, index_col=[0, 1, 2]) for path_to_dataset in paths_to_datasets], axis=0)
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    merge(
        paths_to_datasets=snakemake.input.datasets,
        path_to_output=snakemake.output[0]
    )

import pandas as pd


def to_csv(path_to_feather, path_to_csv):
    (
        pd
        .read_feather(path_to_feather)
        .to_csv(path_to_csv, header=True, index=True)
    )


if __name__ == "__main__":
    to_csv(
        path_to_feather=snakemake.input.feather,
        path_to_csv=snakemake.output[0]
    )

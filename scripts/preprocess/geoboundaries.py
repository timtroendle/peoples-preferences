import geopandas as gpd
import pandas as pd


def merge(paths_to_national_geoboundaries: list[str], path_to_merged_file: str):
    (
        pd
        .concat(
            gpd.read_file(path)
            for path in paths_to_national_geoboundaries
        )
        .to_feather(path_to_merged_file)
    )


if __name__ == "__main__":
    merge(
        paths_to_national_geoboundaries=snakemake.input,
        path_to_merged_file=snakemake.output[0]
    )

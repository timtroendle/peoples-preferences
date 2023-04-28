from pathlib import Path

import pandas as pd
import xarray as xr


def merge_census_features(paths_to_features: list[Path], population_counts: dict[str, str], path_to_output: Path):
    census = (
        xr
        .Dataset({f"frequency_{path.stem}": read_feature(path) for path in paths_to_features})
    )
    census["frequency_countries"] = xr.DataArray( # TODO a bit of a hack, there likely is a better way
        list(population_counts.values()),
        coords={"country": list(population_counts.keys())}
    )
    census["country_of_admin1"] = xr.DataArray(["DEU"] * len(census.admin1), coords={"admin1": census.admin1.values})
    census.to_netcdf(path_to_output)


def read_feature(path: Path) -> pd.DataFrame:
    feature_name = path.stem
    return (
        pd
        .read_feather(path)
        .set_index(feature_name)
        .rename_axis(columns="admin1")
    )


if __name__ == "__main__":
    merge_census_features(
        paths_to_features=[Path(path) for path in snakemake.input.features],
        population_counts=snakemake.params.pop_counts,
        path_to_output=Path(snakemake.output[0])
    )

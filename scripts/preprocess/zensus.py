from io import StringIO

import pandas as pd
import requests


def download_zensus_table(url: str, path_to_output: str):
    r = requests.get(url)
    r.raise_for_status()
    try:
        table_as_str = StringIO(r.json()["Object"]["Content"])
        data = (
            pd
            .read_csv(table_as_str, sep=";", skiprows=4, skipfooter=3, header=[0, 1, 2], engine="python", index_col=0)
            .xs("Anzahl", level=2, axis=1)
            .xs("Personen", level=1, axis=1)
        )
    except KeyError: # some files are a different format
        table_as_str = StringIO(r.json()["Object"]["Content"])
        data = pd.read_csv(
            table_as_str,
            sep=";",
            skiprows=5,
            skipfooter=3,
            header=[0],
            engine="python",
            index_col=0
        )
    data.to_csv(path_to_output, header=True, index=True)


def preprocess_zensus_table(path_to_input: str, dim_name: str, name_mapping: dict[str, str], path_to_output: str):
    (
        pd
        .read_csv(path_to_input, index_col=0)
        .drop(columns="Deutschland", errors="ignore")
        .rename(index=name_mapping, columns=lambda name: name[3:])
        .loc[name_mapping.values(), :]
        .rename_axis(index=dim_name)
        .groupby(level=0, axis=0)
        .sum() # sum over multiple external categories if they match our internal category
        .reset_index()
        .to_feather(path_to_output)
    )


if __name__ == "__main__":
    try:
        download_zensus_table(
            url=snakemake.params.url,
            path_to_output=snakemake.output[0]
        )
    except AttributeError:
        preprocess_zensus_table(
            path_to_input=snakemake.input.data,
            dim_name=snakemake.wildcards.dim,
            name_mapping=snakemake.params.name_mapping,
            path_to_output=snakemake.output[0]
        )

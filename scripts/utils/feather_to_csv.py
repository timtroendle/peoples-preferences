from typing import Union

import pandas as pd


def to_csv(path_to_feather: str, columns: Union[slice, list[str]],
           float_format_per_column: dict[str: str], path_to_csv: str):
    (
        pd
        .read_feather(path_to_feather)
        .loc[:, columns]
        .assign(**{
            col_name: format_float_col(col_name, fmt_str)
            for col_name, fmt_str in float_format_per_column.items()
        })
        .to_csv(path_to_csv, header=True, index=False)
    )


def format_float_col(col: str, fmt_str: str):
    def format_col(df: pd.DataFrame):
        return df[col].map(lambda cell: f"{cell:{fmt_str}}")
    return format_col


def optional_param(name: str, default):
    return snakemake.params[name] if name in snakemake.params.keys() else default


if __name__ == "__main__":
    to_csv(
        path_to_feather=snakemake.input.feather,
        columns=optional_param(name="columns", default=slice(None)),
        float_format_per_column=optional_param(name="float_format_per_column", default={}),
        path_to_csv=snakemake.output[0]
    )

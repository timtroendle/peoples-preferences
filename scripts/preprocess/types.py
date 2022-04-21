import numpy as np
import pandas as pd


LIKERT_TYPE_DEF = {
    "type": "likert",
    "ordered": True,
    "rename": {
        1: "Strongly disagree",
        2: "Disagree",
        3: "Neither agree nor disagree",
        4: "Agree",
        5: "Strongly agree",
    }
}


def data_types(path_to_data: str, types: dict, aggregated_levels: dict, categorised_levels: dict, path_to_output: str):
    (
        pd
        .read_feather(path_to_data)
        .pipe(adjust_types, types)
        .pipe(aggregate_levels, aggregated_levels)
        .pipe(categorise, categorised_levels)
        .to_feather(path_to_output)
    )


def adjust_types(df: pd.DataFrame, types: dict):
    for col, type_def in types.items():
        if type_def["type"] == "factor":
            df[col] = adjust_factor(df[col], type_def)
        elif type_def["type"] == "likert":
            df[col] = adjust_factor(df[col], LIKERT_TYPE_DEF)
        elif type_def["type"] == "numeric":
            df[col] = adjust_numeric(df[col].copy(), min=type_def["min"], max=type_def["max"])
        else:
            raise TypeError(f"Unknown data type: {type_def['type']}")
    return df


def adjust_factor(col: pd.Series, type_def: dict):
    try:
        col = pd.to_numeric(col, errors="raise")
    except ValueError:
        pass # all good, column is not numeric
    if "missing-values" in type_def:
        col.replace(type_def["missing-values"], np.nan, inplace=True)
    col = col.astype("category")
    if "relevel" in type_def:
        col = col.cat.reorder_categories(type_def["relevel"], ordered=type_def["ordered"])
    else:
        col = col.cat.reorder_categories(col.cat.categories, ordered=type_def["ordered"])
    if "rename" in type_def:
        col = col.cat.rename_categories(type_def["rename"])
    return col


def adjust_numeric(col: pd.Series, min: int, max: int):
    invalid_mask = (col < min) | (col > max)
    col.loc[invalid_mask] = None
    return pd.to_numeric(col)


def aggregate_levels(df: pd.DataFrame, aggregated_levels: dict):
    for col, feature_def in aggregated_levels.items():
        new_col_name = col + "_aggregated"
        df[new_col_name] = df[col].astype(str).map(feature_def["mapping"]).astype("category")
        if "ordered" in feature_def:
            df[new_col_name] = df[new_col_name].cat.reorder_categories(feature_def["ordered"], ordered=True)
    return df


def categorise(df: pd.DataFrame, categorised_levels: dict):
    for col, category_def in categorised_levels.items():
        new_col_name = col + "_aggregated"
        df[new_col_name] = pd.cut(
            x=df[col],
            bins=category_def["bins"],
            labels=category_def["labels"],
            ordered=False,
            right=True,
            include_lowest=True
        )
        if "ordered" in category_def:
            df[new_col_name] = df[new_col_name].cat.reorder_categories(category_def["ordered"], ordered=True)
    return df


if __name__ == "__main__":
    data_types(
        path_to_data=snakemake.input.data,
        types=snakemake.params.types,
        aggregated_levels=snakemake.params.aggregated_levels,
        categorised_levels=snakemake.params.categorised_levels,
        path_to_output=snakemake.output[0]
    )

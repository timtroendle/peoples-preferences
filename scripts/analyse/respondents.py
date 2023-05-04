import pandas as pd


def stats(path_to_data: str, country_id: str, path_to_output: str):
    df = pd.read_feather(path_to_data)
    respondents = df.groupby("RESPONDENT_ID").first()
    respondents = respondents[respondents.RESPONDENT_COUNTRY == country_id]

    (
        pd
        .concat([
            shares(respondents, "Q3_GENDER", "Gender"),
            shares(respondents, "Q6_AREA", "Area"),
            shares(respondents, "Q4_BIRTH_YEAR_aggregated", "Age"),
            shares(respondents, "Q9_EDUCATION", "Education")
        ])
        .to_csv(path_to_output, header=True, index=True, float_format="%.1f")
    )


def shares(respondents: pd.DataFrame, feature_col_name: str, feature_name: str) -> pd.Series:
    return (
        respondents[feature_col_name]
        .value_counts(sort=False)
        .div(respondents[feature_col_name].size)
        .mul(100)
        .rename("share (%)")
        .rename_axis(index="level")
        .reset_index()
        .assign(attribute=feature_name)
        .set_index(["attribute", "level"])
    )


if __name__ == "__main__":
    stats(
        path_to_data=snakemake.input.data,
        country_id=snakemake.wildcards.country_id,
        path_to_output=snakemake.output[0]
    )

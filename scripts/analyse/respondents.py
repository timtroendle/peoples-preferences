import pandas as pd


def stats(path_to_data: str, country_id: str, path_to_output: str):
    df = pd.read_feather(path_to_data)
    respondents = df[df.RESPONDENT_COUNTRY == country_id].groupby("RESPONDENT_ID").first()
    gender = respondents["Q3_GENDER"].value_counts(sort=False) / respondents["Q3_GENDER"].size * 100
    area = respondents["Q6_AREA"].value_counts(sort=False) / respondents["Q6_AREA"].size * 100
    income = respondents["Q10_INCOME"].value_counts(sort=False) / respondents["Q10_INCOME"].size * 100

    (
        pd
        .concat([set_stats_index(gender, "Gender"), set_stats_index(area, "Area"), set_stats_index(income, "Income")])
        .to_csv(path_to_output, header=True, index=True, float_format="%.1f")
    )


def set_stats_index(series, attribute_name):
    return (
        series
        .rename("share (%)")
        .rename_axis(index="level")
        .reset_index()
        .assign(attribute=attribute_name)
        .set_index(["attribute", "level"])
    )


if __name__ == "__main__":
    stats(
        path_to_data=snakemake.input.data,
        country_id=snakemake.wildcards.country_id,
        path_to_output=snakemake.output[0]
    )

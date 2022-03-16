import pandas as pd


def stats(path_to_data: str, codes: dict[str, dict[int: str]], path_to_output: str):
    df = pd.read_csv(path_to_data)
    respondents = preprocess_respondents(df, codes)
    gender = respondents["Q3_GENDER"].value_counts(sort=False) / respondents["Q3_GENDER"].size * 100
    area = respondents["Q6_AREA"].value_counts(sort=False) / respondents["Q6_AREA"].size * 100
    income = respondents["Q10_INCOME"].value_counts(sort=False) / respondents["Q10_INCOME"].size * 100

    (
        pd
        .concat([set_stats_index(gender, "Gender"), set_stats_index(area, "Area"), set_stats_index(income, "Income")])
        .to_csv(path_to_output, header=True, index=True, float_format="%.1f")
    )


def preprocess_respondents(df, codes):
    respondents = df.groupby("RESPONDENT_ID").first()
    respondents = respondents.loc[:, codes.keys()]
    for question in codes.keys():
        respondents[question] = (
            respondents[question]
            .astype(int)
            .astype("category")
            .cat
            .rename_categories(codes[question])
        )
    return respondents


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
        codes=snakemake.params.codes,
        path_to_output=snakemake.output[0]
    )

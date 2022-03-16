import pandas as pd

CODES = {
    "Q3_GENDER": {
        1: "female",
        2: "male",
        3: "other"
    },
    "Q6_AREA": {
        1: "rural",
        2: "urban",
        3: "no answer"
    },
    "Q10_INCOME": {
        1: "below 600 EUR",
        2: "600--900 EUR",
        3: "900--1300 EUR",
        4: "1300--1500 EUR",
        5: "1500--2000 EUR",
        6: "2000--2600 EUR",
        7: "2600--3200 EUR",
        8: "3200--4500 EUR",
        9: "4500--6000 EUR",
        10: "6000--10000 EUR",
        11: "above 10000 EUR",
        12: "no answer"
    }
}

def stats(path_to_data, path_to_output):
    df = pd.read_csv(path_to_data)
    respondents = preprocess_respondents(df)
    gender = respondents["Q3_GENDER"].value_counts(sort=False) / respondents["Q3_GENDER"].size * 100
    area = respondents["Q6_AREA"].value_counts(sort=False) / respondents["Q6_AREA"].size * 100
    income = respondents["Q10_INCOME"].value_counts(sort=False) / respondents["Q10_INCOME"].size * 100

    (
        pd
        .concat([set_stats_index(gender, "Gender"), set_stats_index(area, "Area"), set_stats_index(income, "Income")])
        .to_csv(path_to_output, header=True, index=True, float_format="%.1f")
    )


def preprocess_respondents(df):
    respondents = df.groupby("RESPONDENT_ID").first()
    respondents = respondents.rename(columns=language_agnostic_column_name)
    respondents["Q9_EDUCATION_O7"] = replace_text_input_with_1(respondents["Q9_EDUCATION_O7"])
    respondents["Q9_EDUCATION_O8"] = replace_text_input_with_1(respondents["Q9_EDUCATION_O8"]) # FIXME this is not "other"
    respondents["Q12_PARTY_O7"] = replace_text_input_with_1(respondents["Q12_PARTY_O7"])
    if "Q12_PARTY_O8" in respondents.columns:
        respondents["Q12_PARTY_O8"] = replace_text_input_with_1(respondents["Q12_PARTY_O8"])
    if "Q12_PARTY_O9" in respondents.columns:
        respondents["Q12_PARTY_O9"] = replace_text_input_with_1(respondents["Q12_PARTY_O9"])
    if "Q12_PARTY_O10" in respondents.columns:
        respondents["Q12_PARTY_O10"] = replace_text_input_with_1(respondents["Q12_PARTY_O10"])
    respondents = pd.concat([respondents, undummify(respondents[filter(column_is_dummy, respondents.columns)])], axis=1)
    for question in CODES.keys():
        respondents[question] = (
            respondents[question]
            .astype(int)
            .astype("category")
            .cat
            .rename_categories(CODES[question])
        )
    return respondents


def replace_text_input_with_1(series):
    return (
        pd
        .to_numeric(series, errors="coerce")
        .fillna(1)
    )


def language_agnostic_column_name(raw_column_name):
    if not column_is_question_with_options(raw_column_name):
        return raw_column_name
    else:
        elements = raw_column_name.split("_")
        while not (elements[-1].startswith("O") and elements[-1][-1].isnumeric()):
            elements.pop()
        return "_".join(elements)

def column_is_question_with_options(column_name):
    is_question = column_name.startswith("Q")
    has_options = any((element.startswith("O") and element[-1].isnumeric())
                      for element in column_name.split("_"))
    return is_question & has_options


def column_is_dummy(column_name):
    if column_name.startswith("Q15"):
        return False
    else:
        return column_is_question_with_options(column_name)


def undummify(df, prefix_sep="_O"):
    # from https://stackoverflow.com/a/62085741/1856079
    cols2collapse = {
        item.split(prefix_sep)[0]: (prefix_sep in item) for item in df.columns
    }
    series_list = []
    for col, needs_to_collapse in cols2collapse.items():
        if needs_to_collapse:
            undummified = (
                df.filter(like=col)
                .idxmax(axis=1)
                .apply(lambda x: x.split(prefix_sep, maxsplit=1)[1])
                .rename(col)
            )
            series_list.append(undummified)
        else:
            series_list.append(df[col])
    undummified_df = pd.concat(series_list, axis=1)
    return undummified_df


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
        path_to_output=snakemake.output[0]
    )

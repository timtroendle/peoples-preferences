from datetime import date
from itertools import filterfalse

import pandas as pd

INDEX_COLUMNS = ["RESPONDENT_ID", "CHOICE_SET", "LABEL"]
COLUMNS_TO_DROP = [
    "SURVEY_ID", # TODO check what this is
    "RESPONDENT_IP_ADDRESS",
    "RESPONDENT_CITY",
    "RESPONDENT_REGION",
    "RESPONDENT_POSTCODE",
    "RESPONDENT_UNIQUE_CODE",
    "RESPONDENT_TIME_OF_OPENING_SURVEY",
    "RESPONDENT_TIME_OF_COMPLETING_SURVEY",
    "RESPONDENT_LENGTH_OF_INTERVIEW_SECONDS",
    "VARIABLE_LOCALE",
    "PARTICIPANT_STATUSINCLUDEINANALYSIS_TIC",
    "PARTICIPANT_STATUS_INCLUDED_IN_CURRENT_REPORT_TIC",
    "PARTICIPANT_STATUSINCLUDEINANALYSIS_LOCALE",
    "PARTICIPANT_STATUS_INCLUDED_IN_CURRENT_REPORT_LOCALE"
]


def preprocess_conjoint(path_to_conjointly_data: str, path_to_respondi_data: str, country_id: str,
                        population_count: int, pre_test_threshold: date, q12_party_base: int,
                        path_to_output: str):
    respondi = pd.read_excel(path_to_respondi_data, parse_dates=[2])
    (
        pd
        .read_csv(
            path_to_conjointly_data,
            index_col=INDEX_COLUMNS,
            na_values="NULL",
            low_memory=False,
            dtype={"Q5_POSTCODE": str} # FIXME does not work properly for DNK
        )
        .sort_index()
        .drop(columns=COLUMNS_TO_DROP, errors="ignore")
        .fillna({"RESPONDENT_COUNTRY": country_id})
        .pipe(merge_respondi_data, respondi)
        .pipe(undummify_dataset)
        .pipe(filter_pre_test, pre_test_threshold)
        .pipe(shift_q12_party, q12_party_base)
        .assign(WEIGHT=population_count / 1e6) # FIXME must handle the fact that numbers of responses per country vary.
        .to_csv(path_to_output, index=True, header=True)
    )


def merge_respondi_data(df, respondi: pd.DataFrame):
    df = (
        df
        .reset_index() # without this, index gets lost during merge (I do not understand why)
        .merge(respondi, left_on="VARIABLE_TIC", right_on="tic", how="left", validate="many_to_one")
        .drop(columns=["tic"])
        .rename(columns={
            "date_of_last_access": "RESPONDENT_TIME_OF_COMPLETING_SURVEY",
            "duration (in Sek)": "RESPONDENT_LENGTH_OF_INTERVIEW_SECONDS"
        })
        .set_index(INDEX_COLUMNS)
    )
    columns = list(df.columns)
    index_for_duration = columns.index("RESPONDENT_COUNTRY")
    return df.reindex(columns=columns[:index_for_duration + 1] + columns[-2:] + columns[index_for_duration + 1:-2])


def undummify_dataset(df):
    df = df.rename(columns=language_agnostic_column_name)
    df = remove_free_text(df)
    df = pd.concat([df[filterfalse(column_is_dummy, df.columns)], undummify(df[filter(column_is_dummy, df.columns)])], axis=1)
    question_columns = filter(column_is_question, df.columns)
    sorted_question_columns = sorted(question_columns, key=lambda col: int(col.split("_")[0][1:]))
    all_columns_sorted = list(df.columns.drop(sorted_question_columns)) + sorted_question_columns
    return df.reindex(columns=all_columns_sorted)


def filter_pre_test(df, pre_test_threshold):
    return df.loc[df["RESPONDENT_TIME_OF_COMPLETING_SURVEY"] >= pre_test_threshold, :]


def shift_q12_party(df: pd.DataFrame, base: int):
    # Shift values in Q12 column. This is because the meaning of the values in the column is country-specific.
    df["Q12_PARTY"] = pd.to_numeric(df["Q12_PARTY"]) + (base - 1)
    return df


def remove_free_text(df):
    df["Q9_EDUCATION_O7"] = replace_text_input_with_1(df["Q9_EDUCATION_O7"])
    if "Q12_PARTY_O13" in df.columns:
        df["Q12_PARTY_O12"] = replace_text_input_with_1(df["Q12_PARTY_O12"])
    elif "Q12_PARTY_O12" in df.columns:
        df["Q12_PARTY_O11"] = replace_text_input_with_1(df["Q12_PARTY_O11"])
    elif "Q12_PARTY_O11" in df.columns:
        df["Q12_PARTY_O10"] = replace_text_input_with_1(df["Q12_PARTY_O10"])
    elif "Q12_PARTY_O10" in df.columns:
        df["Q12_PARTY_O9"] = replace_text_input_with_1(df["Q12_PARTY_O9"])
    elif "Q12_PARTY_O9" in df.columns:
        df["Q12_PARTY_O8"] = replace_text_input_with_1(df["Q12_PARTY_O8"])
    elif "Q12_PARTY_O8" in df.columns:
        df["Q12_PARTY_O7"] = replace_text_input_with_1(df["Q12_PARTY_O7"])
    return df


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


def column_is_question(column_name):
    return column_name.startswith("Q")


def column_is_question_with_options(column_name):
    has_options = any((element.startswith("O") and element[-1].isnumeric())
                      for element in column_name.split("_"))
    return column_is_question(column_name) & has_options


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


if __name__ == "__main__":
    preprocess_conjoint(
        path_to_conjointly_data=snakemake.input.conjointly,
        path_to_respondi_data=snakemake.input.respondi,
        country_id=snakemake.wildcards.country_id,
        population_count=int(snakemake.params.population),
        pre_test_threshold=snakemake.params.pre_test_threshold,
        q12_party_base=snakemake.params.q12_party_base[snakemake.wildcards.country_id],
        path_to_output=snakemake.output[0]
    )

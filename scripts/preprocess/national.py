from datetime import date
from itertools import filterfalse
from collections.abc import Callable
import re

import pandas as pd

YEAR_OF_STUDY = 2022
MIN_AGE = 18
MAX_AGE = 120
UNREASONABLE_HIGH_AGE = 1000
GERMAN_POSTAL_CODE_AND_CITY = r"^(\d{5})\s([a-zA-ZäöüÄÖÜß\s]+)$"
DANISH_POSTAL_CODE_AND_CITY = r"^(\d{4})\s*([a-zA-ZäöüÄÖÜß\s]+)$"
POLISH_POSTAL_CODE_WITHOUT_DASH = r"^\d\d\d\d\d$"
PORTUGUESE_POSTAL_CODE_WITHOUT_DASH = r"^\d\d\d\d\d\d\d$"
PORTUGUESE_POSTAL_CODE_AND_CITY = r"^(\d{4}-\d{3})\s([a-zA-Záàâãéèêíïóôõöúçñ\s]+)$"
PORTUGUESE_OUTDATED_POSTAL_CODE = r"^\d\d\d\d$"
PORTUGUESE_POSTAL_CODE_NON_DEFAULT_FORMAT = r"^(\d{4})([^-].*|.*[^-])(\d{3})$"
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


def preprocess_conjoint(path_to_conjointly_data: str, path_to_respondi_data: str, path_to_geonames: str,
                        country_id: str,
                        pre_test_threshold: date, q12_party_base: int,
                        path_to_output: str):
    respondi = pd.read_excel(path_to_respondi_data, parse_dates=[2])
    (
        pd
        .read_csv(
            path_to_conjointly_data,
            index_col=INDEX_COLUMNS,
            na_values="NULL",
            low_memory=False,
            dtype={"Q5_POSTCODE": str}
        )
        .sort_index()
        .drop(columns=COLUMNS_TO_DROP, errors="ignore")
        .fillna({"RESPONDENT_COUNTRY": country_id})
        .pipe(merge_respondi_data, respondi)
        .pipe(undummify_dataset)
        .pipe(filter_pre_test, pre_test_threshold)
        .pipe(shift_q12_party, q12_party_base)
        .pipe(fix_broken_entries)
        .pipe(merge_geonames, path_to_geonames, country_id)
        .reset_index()
        .to_feather(path_to_output)
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


def merge_geonames(df: pd.DataFrame, path_to_geonames: str, country_id: str) -> pd.DataFrame:
    geonames = (
        pd
        .read_table(
            path_to_geonames,
            header=None,
            names=[
                "country code", "postal code", "place code", "admin name1", "admin code1",
                "admin name2", "admin code2", "admin name3", "admin code3",
                "latitude", "longitude", "accuracy"
            ],
            index_col=False,
            dtype={
                "postal code": "str",
                "admin name1": "category",
                "admin name2": "category",
                "admin name3": "category"
            }
        )
        .rename(columns={
            "postal code": "RESPONDENT_POSTAL_CODE",
            "admin name1": "RESPONDENT_ADMIN_NAME1",
            "admin name2": "RESPONDENT_ADMIN_NAME2",
            "admin name3": "RESPONDENT_ADMIN_NAME3",
            "latitude": "LATITUDE",
            "longitude": "LONGITUDE"
        })
        .drop(columns={"country code", "admin code1", "admin code2", "admin code3", "place code", "accuracy"})
        .groupby("RESPONDENT_POSTAL_CODE")
        .first() # remove duplicte postal codes (for different place names)
    )
    if country_id == "PRT":
        # ASSUME location precision below 4 digits not necessary in Portugal
        # https://en.wikipedia.org/wiki/Postal_codes_in_Portugal
        geonames_short = geonames.copy()
        geonames_short.index = geonames_short.index.str[:4]
        geonames_short = geonames_short.groupby("RESPONDENT_POSTAL_CODE").first() # use first matching postal code
        geonames = pd.concat([geonames, geonames_short])

        respondents = df.groupby("RESPONDENT_ID").first()
        n_broken = sum(respondents.Q5_POSTCODE.str.match(PORTUGUESE_OUTDATED_POSTAL_CODE))
        print(
            f"Q5_POSTCODE: {n_broken} Portuguese respondentsindicated their outdated, imprecise postal code "
            + "from before 1994. Trying to match with the first matching current postal code (using the first "
            + "four digits of the current postal code only). "
            + "Because current postal codes are much more precise, this will pretend higher precision than "
            + "actually available."
        )
    return (
        df
        .reset_index()
        .merge(geonames, left_on="Q5_POSTCODE", right_on="RESPONDENT_POSTAL_CODE", how="left", validate="many_to_one")
        .set_index(INDEX_COLUMNS)
    )


def undummify_dataset(df):
    df = df.rename(columns=language_agnostic_column_name)
    df = remove_free_text(df)
    df = pd.concat(
        [df[filterfalse(column_is_dummy, df.columns)], undummify(df[filter(column_is_dummy, df.columns)])],
        axis=1
    )
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


class ConjointDataFix:

    def __init__(self, col: str, reason: str, mask: Callable[[pd.Series], pd.Series],
                 fix: Callable[[pd.Series], pd.Series], condition: Callable[[pd.DataFrame], pd.Series] = None):
        self.__col = col
        self.__mask = mask
        self.__fix = fix
        self.__reason = reason
        self.__condition = condition if condition else self.__identity_condition()

    def __identity_condition(self):
        def identity_condition(df: pd.DataFrame):
            return pd.Series(True, index=df.index)
        return identity_condition

    def __call__(self, df):
        if self.sum_broken(df) > 0:
            print(f"{self.__col}: " + self.__reason.format(self.sum_broken(df)))
            data_mask = self.masked(df)
            df.loc[data_mask, self.__col] = self.__fix(df.loc[data_mask, self.__col])
            assert self.sum_broken(df) == 0
        return df

    def masked(self, df):
        condition_mask = self.__condition(df)
        return self.__mask(df.loc[condition_mask, self.__col])

    def sum_broken(self, df):
        respondents = df.groupby("RESPONDENT_ID").first()
        return sum(self.masked(respondents))


def fix_broken_entries(df):
    fixes = [
        ConjointDataFix( # ASSUME negative birth years should be positive.
            col="Q4_BIRTH_YEAR",
            mask=lambda col: col < 0,
            fix=lambda negative: negative.abs(),
            reason="{} respondents indicated negative birth years. These will be considered positive."
        ),
        ConjointDataFix( # ASSUME age entered instead of birth year.
            col="Q4_BIRTH_YEAR",
            mask=lambda col: (col < MAX_AGE) & (col >= MIN_AGE),
            fix=lambda age: YEAR_OF_STUDY - age,
            reason="{} respondents indicated their age instead of birth year. They will be transformed to birth years."
        ),
        ConjointDataFix( # ASSUME negative years in region should be positive.
            col="Q8_YEARS_REGION",
            mask=lambda col: col < 0,
            fix=lambda negative: negative.abs(),
            reason="{} respondents indicated negative years in region. These will be considered positive."
        ),
        ConjointDataFix( # ASSUME calendar years should be number years in region.
            col="Q8_YEARS_REGION",
            mask=lambda col: col > UNREASONABLE_HIGH_AGE,
            fix=lambda calendar_years: YEAR_OF_STUDY - calendar_years,
            reason=(
                "{} respondents indicated calendar years in region instead of number years. "
                + "Will convert to number of years."
            )
        ),
        ConjointDataFix(
            col="Q5_POSTCODE",
            mask=lambda col: col.str.match(GERMAN_POSTAL_CODE_AND_CITY),
            condition=lambda df: df.RESPONDENT_COUNTRY == "DEU",
            fix=lambda col: [re.match(GERMAN_POSTAL_CODE_AND_CITY, postal_code).groups()[0] for postal_code in col],
            reason=("{} German respondents indicated their postal code and city. Retaining code only.")
        ),
        ConjointDataFix(
            col="Q5_POSTCODE",
            mask=lambda col: col.str.match(POLISH_POSTAL_CODE_WITHOUT_DASH),
            condition=lambda df: df.RESPONDENT_COUNTRY == "POL",
            fix=lambda col: [str(postal_code)[:2] + "-" + str(postal_code)[2:] for postal_code in col],
            reason=("{} Polish respondents indicated their postal code without dash. Adding dash.")
        ),
        ConjointDataFix(
            col="Q5_POSTCODE",
            mask=lambda col: col.str.match(PORTUGUESE_POSTAL_CODE_WITHOUT_DASH),
            condition=lambda df: df.RESPONDENT_COUNTRY == "PRT",
            fix=lambda col: [str(postal_code)[:4] + "-" + str(postal_code)[4:] for postal_code in col],
            reason=("{} Portuguese respondents indicated their postal code without dash. Adding dash.")
        ),
        ConjointDataFix(
            col="Q5_POSTCODE",
            mask=lambda col: col.str.match(PORTUGUESE_POSTAL_CODE_AND_CITY),
            condition=lambda df: df.RESPONDENT_COUNTRY == "PRT",
            fix=lambda col: [re.match(PORTUGUESE_POSTAL_CODE_AND_CITY, postal_code).groups()[0]
                             for postal_code in col],
            reason=("{} Portuguese respondents indicated their postal code and city. Retaining code only.")
        ),
        ConjointDataFix(
            col="Q5_POSTCODE",
            mask=lambda col: col.str.match(PORTUGUESE_POSTAL_CODE_NON_DEFAULT_FORMAT),
            condition=lambda df: df.RESPONDENT_COUNTRY == "PRT",
            fix=lambda col: [streamline_non_default_portuguese_postal_code_format(postal_code) for postal_code in col],
            reason=("{} Portuguese respondents indicated their postal code in a non-default format. Streamlining.")
        ),
        ConjointDataFix(
            col="Q5_POSTCODE",
            mask=lambda col: col.str.match(DANISH_POSTAL_CODE_AND_CITY),
            condition=lambda df: df.RESPONDENT_COUNTRY == "DNK",
            fix=lambda col: [re.match(DANISH_POSTAL_CODE_AND_CITY, postal_code).groups()[0] for postal_code in col],
            reason=("{} Danish respondents indicated their postal code and city. Retaining code only.")
        ),
    ]

    for fix in fixes:
        df = fix(df)

    return df


def streamline_non_default_portuguese_postal_code_format(postal_code_non_default_format: str) -> str:
    re_groups = re.match(PORTUGUESE_POSTAL_CODE_NON_DEFAULT_FORMAT, postal_code_non_default_format).groups()
    return re_groups[0] + "-" + re_groups[2]


if __name__ == "__main__":
    preprocess_conjoint(
        path_to_conjointly_data=snakemake.input.conjointly,
        path_to_respondi_data=snakemake.input.respondi,
        path_to_geonames=snakemake.input.geonames,
        country_id=snakemake.wildcards.country_id,
        pre_test_threshold=snakemake.params.pre_test_threshold,
        q12_party_base=snakemake.params.q12_party_base[snakemake.wildcards.country_id],
        path_to_output=snakemake.output[0]
    )

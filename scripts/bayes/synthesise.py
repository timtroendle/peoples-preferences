import itertools

import numpy as np
import pandas as pd


ATTRIBUTES = { # TODO take from config
    "TECHNOLOGY": pd.Series([
        "Rooftop PV",
        "Open-field PV",
        "Wind"
    ]),
    "LAND": pd.Series([
        "0.5%",
        "1%",
        "2%",
        "4%",
        "8%"
    ]),
    "TRANSMISSION": pd.Series([
        "-25.0% .",
        "+0% .",
        "+25% .",
        "+50% .",
        "+75% ."
    ]),
    "SHARE_IMPORTS": pd.Series([
        "0%",
        "10%",
        "50%",
        "90%"
    ]),
    "PRICES": pd.Series([
        "+0%",
        "+15%",
        "+30%",
        "+45%",
        "+60%"
    ]),
    "OWNERSHIP": pd.Series([
        "Public",
        "Community",
        "Private"
    ])
}
COUNTRIES = pd.Series(["DEU", "DNK", "POL", "PRT"])


def create_synthetic_data(path_to_original_data: str, n_respondents: int, n_sets_per_respondent: int,
                          seed: int, path_to_output: str):
    random_generator = np.random.default_rng(seed)
    choice_sets = (
        pd
        .concat([create_choice_sets(n_sets_per_respondent, random_generator) for r in range(n_respondents)])
        .reset_index(drop=True)
    )
    respondents = create_respondents(n_respondents, random_generator)
    respondents = pd.DataFrame({
        column: respondents.loc[:, column].repeat(choice_sets.shape[0] // n_respondents)
        for column in respondents.columns
    }).reset_index(drop=True)
    (
        pd
        .concat([respondents, choice_sets], axis=1)
        .pipe(match_dtypes, original=pd.read_feather(path_to_original_data))
        .to_feather(path_to_output)
    )


def create_random_profiles(n_profiles: int, random_generator: np.random.Generator) -> pd.DataFrame:
    return pd.DataFrame({
        attribute: levels.sample(n_profiles, replace=True, ignore_index=True, random_state=random_generator)
        for attribute, levels in ATTRIBUTES.items()
    })


def create_choice_sets(n_sets: int, random_generator: np.random.Generator) -> pd.DataFrame:
    profiles = create_random_profiles(n_sets * 2, random_generator)
    return (
        profiles
        .assign(
            CHOICE_SET=list(itertools.chain.from_iterable(itertools.repeat(choice_set, 2)
                                                          for choice_set in range(1, n_sets + 1))),
            LABEL=list(itertools.chain.from_iterable(itertools.repeat(["Left", "Right"], n_sets)))
        )
    )


def create_respondents(n_respondents, random_generator: np.random.Generator) -> pd.DataFrame:
    return pd.DataFrame({
        "RESPONDENT_ID": range(n_respondents),
        "RESPONDENT_COUNTRY": COUNTRIES.sample(
            n_respondents,
            replace=True,
            ignore_index=True,
            random_state=random_generator
        )
    })


def match_dtypes(data: pd.DataFrame, original: pd.DataFrame) -> pd.DataFrame:
    for column in data.columns:
        if column != "RESPONDENT_ID":
            data[column] = data[column].astype(original[column].dtype)
        else:
            data[column] = data[column].astype("category")
    return data


if __name__ == "__main__":
    create_synthetic_data(
        path_to_original_data=snakemake.input.original,
        n_respondents=snakemake.params.n_respondents,
        n_sets_per_respondent=snakemake.params.n_sets_per_respondent,
        seed=snakemake.params.seed,
        path_to_output=snakemake.output[0]
    )

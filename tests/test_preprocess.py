import pandas as pd

MINIMAL_PARTISANSHIP_CODE = {
    "DEU": 1,
    "DNK": 9,
    "POL": 22,
    "PRT": 34
}

MAXIMAL_PARTISANSHIP_CODE = {
    "DEU": 8,
    "DNK": 21,
    "POL": 33,
    "PRT": 46
}


def test_minimal_partisanship_code(country_id, national_conjoint):
    assert national_conjoint["Q12_PARTY"].min() >= MINIMAL_PARTISANSHIP_CODE[country_id]


def test_maximal_partisanship_code(country_id, national_conjoint):
    assert national_conjoint["Q12_PARTY"].max() <= MAXIMAL_PARTISANSHIP_CODE[country_id]


def test_imputation_not_destructive(conjoint_imputed: pd.DataFrame, imputed_feature_name: str):
    feature = conjoint_imputed[imputed_feature_name]
    feature_imputed = conjoint_imputed[f"{imputed_feature_name}_imputed"]
    not_missing = ~feature.isnull()
    assert all(feature[not_missing] == feature_imputed[not_missing])

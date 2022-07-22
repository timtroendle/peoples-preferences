import numpy as np
from datatest import register_accessors

register_accessors()


def test_gender_covariate(covariate_model, respondents):
    model_gender = (
        covariate_model
        .constant_data
        .g
        .to_dataframe()
        .reset_index()
        .pivot(index="respondent", columns="gender", values="g")
        .assign(Female=0)
    )
    model_gender.loc[(model_gender.Male == 0) & (model_gender.Other == 0), "Female"] = 1
    model_gender = model_gender.idxmax(axis=1)

    expected_gender = respondents.reindex(model_gender.index).Q3_GENDER

    model_gender.validate(expected_gender)


def test_education_covariate(covariate_model, respondents):
    levels = ["None"] + list(covariate_model.posterior.education.values)
    level_mapping = {i: level for i, level in enumerate(levels)}
    model_education = covariate_model.constant_data.edu.to_series().map(level_mapping)

    expected_education = respondents.reindex(model_education.index).Q9_EDUCATION

    model_education.validate(expected_education)


def test_age_covariate(covariate_model, respondents):
    model_age = covariate_model.constant_data.age.to_series()

    expected_age = 2022 - respondents.reindex(model_age.index).Q4_BIRTH_YEAR

    model_age.validate(expected_age)

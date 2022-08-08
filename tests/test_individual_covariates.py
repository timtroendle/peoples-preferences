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


def test_country_covariate(covariate_model, respondents):
    model_country = (
        covariate_model
        .constant_data
        .c
        .to_dataframe()
        .reset_index()
        .pivot(index="respondent", columns="country", values="c")
        .assign(DEU=0)
    )
    model_country.loc[(model_country.DNK == 0) & (model_country.PRT == 0) & (model_country.POL == 0), "DEU"] = 1
    model_country = model_country.idxmax(axis=1)

    expected_country = respondents.reindex(model_country.index).RESPONDENT_COUNTRY

    model_country.validate(expected_country)


def test_area_covariate(covariate_model, respondents):
    model_area = (
        covariate_model
        .constant_data
        .a
        .to_dataframe()
        .reset_index()
        .pivot(index="respondent", columns="area", values="a")
        .assign(Rural=0)
    )
    model_area.loc[(model_area.Urban == 0), "Rural"] = 1
    model_area = model_area.idxmax(axis=1)

    expected_area = respondents.reindex(model_area.index).Q6_AREA

    model_area.validate(expected_area)


def test_renewables_covariate(covariate_model, respondents):
    model_renewables = (
        covariate_model
        .constant_data
        .re
        .to_dataframe()
        .reset_index()
        .pivot(index="respondent", columns="renewables", values="re")
    )
    model_renewables["Open-field PV"] = 0
    pv = (model_renewables.Both == 0) & (model_renewables.Wind == 0) % model_renewables.Neither == 0
    model_renewables.loc[pv, "Open-field PV"] = 1
    model_renewables = model_renewables.idxmax(axis=1)

    expected_renewables = respondents.reindex(model_renewables.index).Q7_RENEWABLES

    model_renewables.validate(expected_renewables)


def test_party_covariate(covariate_model, respondents):
    model_party = (
        covariate_model
        .constant_data
        .p
        .to_dataframe()
        .reset_index()
        .pivot(index="respondent", columns="party", values="p")
        .assign(Other=0)
    )
    other = (model_party.Green == 0) & (model_party.Right == 0)
    model_party.loc[other, "Other"] = 1
    model_party = model_party.idxmax(axis=1)

    expected_party = respondents.reindex(model_party.index).Q12_PARTY_aggregated

    model_party.validate(expected_party)


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

# TODO add income
# TODO add years
# TODO add climate concern

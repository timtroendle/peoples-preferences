import numpy as np
import pandas as pd
import pymc as pm

ATTRIBUTES = [
    "TECHNOLOGY",
    "LAND",
    "TRANSMISSION",
    "SHARE_IMPORTS",
    "PRICES",
    "OWNERSHIP"
]
OPTIONS_PER_RESPONDENT = 16
YEAR_OF_SURVEY = 2022


def hierarchical_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, limit_respondents: bool,
                       n_respondents_per_country: int, random_seed: int, individual_covariates: bool,
                       path_to_output: str):
    conjoint = (
        pd
        .read_feather(path_to_data)
        .set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
        .pipe(filter_respondents, limit_respondents, n_respondents_per_country)
        .pipe(prepare_respondent_age)
        .pipe(prepare_years_region)
    )
    dummies = pd.get_dummies(conjoint.loc[:, ATTRIBUTES], drop_first=True, prefix_sep=":")

    model = pm.Model(coords={
        "level": dummies.columns.values,
        "level_repeat": dummies.columns.values,
        "respondent": conjoint.index.get_level_values("RESPONDENT_ID").remove_unused_categories().categories,
        "education": pd.get_dummies(conjoint.Q9_EDUCATION.cat.remove_unused_categories(), drop_first=True).columns.values,
        "gender": pd.get_dummies(conjoint.Q3_GENDER, drop_first=True).columns.values,
        "country": conjoint.RESPONDENT_COUNTRY.cat.categories,
        "area": pd.get_dummies(conjoint.Q6_AREA, drop_first=True).columns.values,
        "renewables": pd.get_dummies(conjoint.Q7_RENEWABLES, drop_first=True).columns.values,
        "party": pd.get_dummies(conjoint.Q12_PARTY_aggregated, drop_first=True).columns.values,
        "income": pd.get_dummies(conjoint.Q10_INCOME.cat.remove_unused_categories(), drop_first=True).columns.values,
        "concern": pd.get_dummies(conjoint.Q11_CLIMATE_CONCERN.cat.remove_unused_categories(), drop_first=True).columns.values
    })

    n_levels = len(model.coords["level"])
    n_respondents = len(model.coords["respondent"])
    n_educations = len(model.coords["education"])
    n_incomes = len(model.coords["income"])
    n_concerns = len(model.coords["concern"])

    with model:
        # data
        r = pm.ConstantData(
            "r",
            conjoint.xs("Left", level="LABEL").index.get_level_values("RESPONDENT_ID").remove_unused_categories().codes,
            dims="choice_situation"
        )
        dummies_left = pm.ConstantData(
            "dummies_left",
            dummies.xs("Left", level="LABEL").values,
            dims=["choice_situation", "level"]
        )
        dummies_right = pm.ConstantData(
            "dummies_right",
            dummies.xs("Right", level="LABEL").values,
            dims=["choice_situation", "level"]
        )
        choice_left = pm.ConstantData(
            "choice_left",
            conjoint.xs("Left", level="LABEL").CHOICE_INDICATOR.values,
            dims="choice_situation"
        )
        age = pm.ConstantData(
            "age",
            conjoint.groupby("RESPONDENT_ID").RESPONDENT_AGE.first().values,
            dims="respondent"
        )
        age_normed = pm.ConstantData(
            "age_normed",
            conjoint.groupby("RESPONDENT_ID").RESPONDENT_AGE_NORM.first().values.repeat(n_levels).reshape(n_respondents, n_levels).T,
            dims=["level", "respondent"] # FIXME this should be respondent only
        )
        years = pm.ConstantData(
            "years",
            conjoint.groupby("RESPONDENT_ID").Q8_YEARS_REGION.first().values,
            dims="respondent"
        )
        years_normed = pm.ConstantData(
            "years_normed",
            conjoint.groupby("RESPONDENT_ID").Q8_YEARS_REGION_NORM.first().values.repeat(n_levels).reshape(n_respondents, n_levels).T,
            dims=["level", "respondent"] # FIXME this should be respondent only
        )
        edu = pm.ConstantData(
            "edu",
            conjoint.groupby("RESPONDENT_ID").Q9_EDUCATION.first().cat.remove_unused_categories().cat.codes.values,
            dims="respondent"
        )
        i = pm.ConstantData(
            "i",
            conjoint.groupby("RESPONDENT_ID").Q10_INCOME.first().cat.remove_unused_categories().cat.codes.values,
            dims="respondent"
        )
        cc = pm.ConstantData(
            "cc",
            conjoint.groupby("RESPONDENT_ID").Q11_CLIMATE_CONCERN.first().cat.remove_unused_categories().cat.codes.values,
            dims="respondent"
        )
        g = pm.ConstantData(
            "g",
            pd.get_dummies(conjoint.groupby("RESPONDENT_ID").Q3_GENDER.first(), drop_first=True).values,
            dims=["respondent", "gender"]
        )
        c = pm.ConstantData(
            "c",
            conjoint.groupby("RESPONDENT_ID").RESPONDENT_COUNTRY.first().cat.codes,
            dims=["respondent"]
        )
        a = pm.ConstantData(
            "a",
            pd.get_dummies(conjoint.groupby("RESPONDENT_ID").Q6_AREA.first(), drop_first=True).values,
            dims=["respondent", "area"]
        )
        re = pm.ConstantData(
            "re",
            pd.get_dummies(conjoint.groupby("RESPONDENT_ID").Q7_RENEWABLES.first(), drop_first=True).values,
            dims=["respondent", "renewables"]
        )
        p = pm.ConstantData(
            "p",
            pd.get_dummies(conjoint.groupby("RESPONDENT_ID").Q12_PARTY_aggregated.first(), drop_first=True).values,
            dims=["respondent", "party"]
        )

        # parameters
        alpha = pm.Normal('alpha', 0, sigma=4, dims="level")

        chol_country, rho_country, sigma_country = pm.LKJCholeskyCov(
            "chol_country",
            n=len(model.coords["level"]),
            eta=4,
            sd_dist=pm.Exponential.dist(1.0),
            compute_corr=True,
            store_in_trace=False
        )
        pm.Deterministic("sigma_country", sigma_country, dims="level")
        pm.Deterministic("rho_country", rho_country, dims=["level", "level_repeat"])

        z_country = pm.Normal("z_country", 0.0, 1.0, dims=["level", "country"])
        countries = pm.Deterministic("countries", pm.math.dot(chol_country, z_country), dims=["level", "country"])

        country_effect = pm.Deterministic(
            'country_effect',
            countries[:, c].T,
            dims=["respondent", "level"]
        )

        chol_individuals, rho_individuals, sigma_individuals = pm.LKJCholeskyCov( # TODO should these be per-country?
            "chol_individuals",
            n=len(model.coords["level"]),
            eta=4,
            sd_dist=pm.Exponential.dist(1.0),
            compute_corr=True,
            store_in_trace=False
        )
        pm.Deterministic("sigma_individuals", sigma_individuals, dims="level")
        pm.Deterministic("rho_individuals", rho_individuals, dims=["level", "level_repeat"])

        z_individuals = pm.Normal("z_individuals", 0.0, 1.0, dims=["level", "respondent"])
        individuals = pm.Deterministic("individuals", pm.math.dot(chol_individuals, z_individuals), dims=["level", "respondent"])

        if individual_covariates:
            alpha_gender = pm.Normal('alpha_gender', mu=0, sigma=1, dims=["gender", "level"]) # TODO add covariation?
            alpha_area = pm.Normal('alpha_area', mu=0, sigma=1, dims=["area", "level"]) # TODO add covariation?
            alpha_renewables = pm.Normal('alpha_renewables', mu=0, sigma=1, dims=["renewables", "level"]) # TODO add covariation?
            alpha_party = pm.Normal('alpha_party', mu=0, sigma=1, dims=["party", "level"]) # TODO add covariation?
            beta_age_normed = pm.Normal('beta_age_normed', mu=0, sigma=1, dims="level") # TODO add covariation?
            beta_years_normed = pm.Normal('beta_years_normed', mu=0, sigma=1, dims="level") # TODO add covariation?
            beta_edu = pm.Normal("beta_edu", mu=0, sigma=1, dims="level") # TODO add covariation?
            beta_income = pm.Normal("beta_income", mu=0, sigma=1, dims="level") # TODO add covariation?
            beta_concern = pm.Normal("beta_concern", mu=0, sigma=1, dims="level") # TODO add covariation?
            d_edu = pm.Dirichlet("d_edu", a=np.ones([n_levels, n_educations]) * 2, transform=pm.distributions.transforms.simplex, dims=["level", "education"])
            d_edu_cumsum = pm.math.stack([pm.math.sum(d_edu[:, :i], axis=1) for i in range(n_educations + 1)])[edu, :]
            d_income = pm.Dirichlet("d_income", a=np.ones([n_levels, n_incomes]) * 2, transform=pm.distributions.transforms.simplex, dims=["level", "income"])
            d_income_cumsum = pm.math.stack([pm.math.sum(d_income[:, :i], axis=1) for i in range(n_incomes + 1)])[i, :]
            d_concern = pm.Dirichlet("d_concern", a=np.ones([n_levels, n_concerns]) * 2, transform=pm.distributions.transforms.simplex, dims=["level", "concern"])
            d_concern_cumsum = pm.math.stack([pm.math.sum(d_concern[:, :i], axis=1) for i in range(n_concerns + 1)])[cc, :]

            gender_effect = pm.Deterministic(
                'gender_effect',
                pm.math.dot(g, alpha_gender),
                dims=["respondent", "level"]
            )
            area_effect = pm.Deterministic(
                'area_effect',
                pm.math.dot(a, alpha_area),
                dims=["respondent", "level"]
            )
            renewables_effect = pm.Deterministic(
                'renewables_effect',
                pm.math.dot(re, alpha_renewables),
                dims=["respondent", "level"]
            )
            party_effect = pm.Deterministic(
                'party_effect',
                pm.math.dot(p, alpha_party),
                dims=["respondent", "level"]
            )
            age_effect = pm.Deterministic(
                'age_effect',
                beta_age_normed * age_normed.T,
                dims=["respondent", "level"]
            )
            years_effect = pm.Deterministic(
                'years_effect',
                beta_years_normed * years_normed.T,
                dims=["respondent", "level"]
            )
            edu_effect = pm.Deterministic(
                'edu_effect',
                beta_edu * d_edu_cumsum,
                dims=["respondent", "level"]
            )
            income_effect = pm.Deterministic(
                'income_effect',
                beta_income * d_income_cumsum,
                dims=["respondent", "level"]
            )
            concern_effect = pm.Deterministic(
                'concern_effect',
                beta_concern * d_concern_cumsum,
                dims=["respondent", "level"]
            )

            partworths = pm.Deterministic(
                "partworths",
                alpha + gender_effect + country_effect + area_effect + renewables_effect
                + party_effect + age_effect + years_effect + edu_effect + income_effect + concern_effect
                + individuals.T,
                dims=["respondent", "level"]
            )
        else:
            partworths = pm.Deterministic(
                "partworths",
                alpha + country_effect + individuals.T,
                dims=["respondent", "level"]
            )

        mu_left_intercept = pm.Normal('mu_left_intercept', 0, sigma=4)
        sigma_left_intercept = pm.Exponential('sigma_left_intercept', 1)
        z_left_intercept = pm.Normal("z_left_intercept", 0, sigma=1, dims="respondent")
        left_intercept = pm.Deterministic( # TODO let covar with partsworths
            'left_intercept',
            mu_left_intercept + z_left_intercept * sigma_left_intercept,
            dims="respondent"
        )

        u_left = left_intercept[r] + pm.math.sum(partworths[r, :] * dummies_left, axis=1)
        u_right = pm.math.sum(partworths[r, :] * dummies_right, axis=1)
        p_left = pm.Deterministic(
            "p_left",
            pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right)),
            dims="choice_situation"
        )

        pm.Bernoulli(
            f"choice",
            p=p_left,
            observed=choice_left,
            dims="choice_situation"
        )

        inference_data = pm.sample(
            draws=n_draws,
            tune=n_tune,
            cores=n_cores,
            random_seed=random_seed,
            return_inferencedata=True,
            target_accept=0.9
        )
    inference_data.to_netcdf(path_to_output)


def filter_respondents(df, limit_respondents, n_respondents_per_country):
    df.dropna( # FIXME don't do this
        axis="index",
        subset=["Q3_GENDER", "Q4_BIRTH_YEAR", "Q6_AREA", "Q7_RENEWABLES", "Q8_YEARS_REGION",
                "Q9_EDUCATION", "Q10_INCOME", "Q11_CLIMATE_CONCERN", "Q12_PARTY_aggregated"],
        inplace=True
    )
    if limit_respondents:
        df = df.groupby("RESPONDENT_COUNTRY").head(n_respondents_per_country * OPTIONS_PER_RESPONDENT)
    return df


def prepare_respondent_age(df):
    df["RESPONDENT_AGE"] = YEAR_OF_SURVEY - df["Q4_BIRTH_YEAR"].astype("int")
    df["RESPONDENT_AGE_NORM"] = (df["RESPONDENT_AGE"] - df["RESPONDENT_AGE"].mean()) / df["RESPONDENT_AGE"].std()
    return df


def prepare_years_region(df):
    df["Q8_YEARS_REGION_NORM"] = (df["Q8_YEARS_REGION"] - df["Q8_YEARS_REGION"].mean()) / df["Q8_YEARS_REGION"].std()
    return df


if __name__ == "__main__":
    hierarchical_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        limit_respondents=bool(snakemake.params.limit_respondents),
        n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
        n_cores=snakemake.threads,
        random_seed=snakemake.params.random_seed,
        individual_covariates=snakemake.params.individual_covariates,
        path_to_output=snakemake.output[0]
    )

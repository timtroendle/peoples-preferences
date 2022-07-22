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
    )
    dummies = pd.get_dummies(conjoint.loc[:, ATTRIBUTES], drop_first=True, prefix_sep=":")

    model = pm.Model(coords={
        "level": dummies.columns.values,
        "level_repeat": dummies.columns.values,
        "respondent": conjoint.index.get_level_values("RESPONDENT_ID").remove_unused_categories().categories,
        "education": pd.get_dummies(conjoint.Q9_EDUCATION.cat.remove_unused_categories(), drop_first=True).columns.values,
        "gender": pd.get_dummies(conjoint.Q3_GENDER, drop_first=True).columns.values
    })

    n_levels = len(model.coords["level"])
    n_respondents = len(model.coords["respondent"])
    n_educations = len(model.coords["education"])


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
        edu = pm.ConstantData(
            "edu",
            conjoint.groupby("RESPONDENT_ID").Q9_EDUCATION.first().cat.remove_unused_categories().cat.codes.values,
            dims="respondent"
        )
        g = pm.ConstantData(
            "g",
            pd.get_dummies(conjoint.groupby("RESPONDENT_ID").Q3_GENDER.first(), drop_first=True).values,
            dims=["respondent", "gender"]
        )

        # parameters
        alpha = pm.Normal('alpha', 0, sigma=4, dims="level")
        chol_individuals, rho_individuals, sigma_individuals = pm.LKJCholeskyCov(
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
            alpha_gender = pm.Normal('alpha_gender', mu=0, sigma=1, dims=["gender", "level"])
            beta_age_normed = pm.Normal('beta_age_normed', mu=0, sigma=1, dims="level") # TODO add covariation?
            beta_edu = pm.Normal("beta_edu", mu=0, sigma=1, dims="level") # TODO add covariation?
            d_edu = pm.Dirichlet("d_edu", a=np.ones([n_levels, n_educations]) * 2, transform=pm.distributions.transforms.simplex, dims=["level", "education"])
            d_edu_cumsum = pm.math.stack([pm.math.sum(d_edu[:, :i], axis=1) for i in range(n_educations + 1)])[edu, :]

            partworths = pm.Deterministic(
                "partworths",
                alpha + pm.math.dot(g, alpha_gender) + beta_age_normed * age_normed.T + beta_edu * d_edu_cumsum + individuals.T,
                dims=["respondent", "level"]
            )
        else:
            partworths = pm.Deterministic(
                "partworths",
                alpha + individuals.T,
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
        p_left = pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right))

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
    df.dropna(axis="index", subset=["Q3_GENDER", "Q4_BIRTH_YEAR", "Q9_EDUCATION"], inplace=True) # FIXME don't do this
    if limit_respondents:
        df = df.groupby("RESPONDENT_COUNTRY").head(n_respondents_per_country * OPTIONS_PER_RESPONDENT)
    return df


def prepare_respondent_age(df):
    df["RESPONDENT_AGE"] = YEAR_OF_SURVEY - df["Q4_BIRTH_YEAR"].astype("int")
    df["RESPONDENT_AGE_NORM"] = (df["RESPONDENT_AGE"] - df["RESPONDENT_AGE"].mean()) / df["RESPONDENT_AGE"].std()
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

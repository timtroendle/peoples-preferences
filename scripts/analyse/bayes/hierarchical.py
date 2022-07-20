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


def hierarchical_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, limit_respondents: bool,
                       n_respondents_per_country: int, random_seed: int, path_to_output: str):
    conjoint = (
        pd
        .read_feather(path_to_data)
        .set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
    )
    if limit_respondents:
        conjoint = conjoint.groupby("RESPONDENT_COUNTRY").head(n_respondents_per_country * OPTIONS_PER_RESPONDENT)
    dummies = pd.get_dummies(conjoint.loc[:, ATTRIBUTES], drop_first=True, prefix_sep=":")

    model = pm.Model(coords={
        "level": dummies.columns.values,
        "level_repeat": dummies.columns.values,
        "respondent": conjoint.index.get_level_values("RESPONDENT_ID").remove_unused_categories().categories
    })

    with model:
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
        partworths = pm.Deterministic("partworths", (alpha + individuals.T).T, dims=["level", "respondent"])

        mu_left_intercept = pm.Normal('mu_left_intercept', 0, sigma=4)
        sigma_left_intercept = pm.Exponential('sigma_left_intercept', 1)
        z_left_intercept = pm.Normal("z_left_intercept", 0, sigma=1, dims="respondent")
        left_intercept = pm.Deterministic( # TODO let covar with partsworths
            'left_intercept',
            mu_left_intercept + z_left_intercept * sigma_left_intercept,
            dims="respondent"
        )

        r = pm.ConstantData(
            "r",
            conjoint.xs("Left", level="LABEL").index.get_level_values("RESPONDENT_ID").remove_unused_categories().codes,
            dims="choice_situations"
        )
        dummies_left = pm.ConstantData(
            "dummies_left",
            dummies.xs("Left", level="LABEL").values,
            dims=["choice_situations", "level"]
        )
        dummies_right = pm.ConstantData(
            "dummies_right",
            dummies.xs("Right", level="LABEL").values,
            dims=["choice_situations", "level"]
        )
        choice_left = pm.ConstantData(
            "choice_left",
            conjoint.xs("Left", level="LABEL").CHOICE_INDICATOR.values,
            dims="choice_situations"
        )

        u_left = left_intercept[r] + pm.math.sum(partworths[:, r] * dummies_left.T, axis=0)
        u_right = pm.math.sum(partworths[:, r] * dummies_right.T, axis=0)
        p_left = pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right))

        pm.Bernoulli(
            f"choice",
            p=p_left,
            observed=choice_left,
            dims="choice_situations"
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


if __name__ == "__main__":
    hierarchical_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        limit_respondents=bool(snakemake.params.limit_respondents),
        n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
        n_cores=snakemake.threads,
        random_seed=snakemake.params.random_seed,
        path_to_output=snakemake.output[0]
    )

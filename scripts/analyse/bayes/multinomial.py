import pandas as pd
import pymc as pm
import arviz as az

ATTRIBUTES = [
    "CHOICE_INDICATOR",
    "TECHNOLOGY",
    "LAND",
    "PRICES",
    "TRANSMISSION",
    "OWNERSHIP",
    "SHARE_IMPORTS"
]

OPTIONS_PER_RESPONDENT = 16


def multinomial_logit_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, limit_respondents: bool,
                            n_respondents_per_country: int, random_seed: int, path_to_output: str):
    conjoint = (
        pd
        .read_feather(path_to_data)
        .set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
    )
    if limit_respondents:
        conjoint = conjoint.groupby("RESPONDENT_COUNTRY").head(n_respondents_per_country)
    pure_conjoint = pd.get_dummies(conjoint.loc[:, ATTRIBUTES], drop_first=True, prefix_sep=":")
    left_dummies = pure_conjoint.xs("Left", level="LABEL").drop(columns="CHOICE_INDICATOR")
    right_dummies = pure_conjoint.xs("Right", level="LABEL").drop(columns="CHOICE_INDICATOR")
    choice_left = pure_conjoint.xs("Left", level="LABEL").CHOICE_INDICATOR

    def indices_of_country(country, label):
        return conjoint.RESPONDENT_COUNTRY.xs(label, level="LABEL") == country

    model = pm.Model(coords={
        "level": pure_conjoint.drop(columns="CHOICE_INDICATOR").columns.values,
        "country": conjoint.RESPONDENT_COUNTRY.unique()
    })

    with model:
        partworths = pm.Normal("partworths", 0, 4, dims=["level", "country"]) # TODO N or MVN?
        intercept_left = pm.Normal("intercept", 0, 4, dims=["country"])

        for country_id, country in enumerate(model.coords["country"]):
            u_left = intercept_left[country_id] + pm.math.sum(
                left_dummies[indices_of_country(country, "Left")].values * partworths[:, country_id], axis=1
            )
            u_right = pm.math.sum(
                right_dummies[indices_of_country(country, "Right")].values * partworths[:, country_id], axis=1
            )
            p_left = pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right))

            choices = pm.Bernoulli(
                f"choice_{country}",
                p=p_left,
                observed=choice_left[indices_of_country(country, "Left")]
            )

        inference_data = pm.sample(
            draws=n_draws,
            tune=n_tune,
            cores=n_cores,
            random_seed=random_seed,
            return_inferencedata=True
        )

    inference_data.to_netcdf(path_to_output)


if __name__ == "__main__":
    multinomial_logit_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        n_cores=snakemake.threads,
        limit_respondents=bool(snakemake.params.limit_respondents),
        n_respondents_per_country=int(snakemake.params.limit_respondents),
        random_seed=snakemake.params.random_seed,
        path_to_output=snakemake.output[0]
    )

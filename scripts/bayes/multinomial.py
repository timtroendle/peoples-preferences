import pandas as pd
import pymc as pm

ATTRIBUTES = [
    "TECHNOLOGY",
    "LAND",
    "PRICES",
    "TRANSMISSION",
    "OWNERSHIP",
    "SHARE_IMPORTS"
]

OPTIONS_PER_RESPONDENT = 16


def multinomial_logit_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, limit_respondents: bool,
                            n_respondents_per_country: int, random_seed: int, sample_type: str,
                            path_to_output: str):
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
        "country": conjoint.RESPONDENT_COUNTRY.cat.categories
    })

    with model:
        partworths = pm.Normal("partworths", 0, 4, dims=["country", "level"]) # TODO N or MVN?
        intercept_left = pm.Normal("intercept", 0, 4, dims=["country"])

        c = pm.ConstantData(
            "c",
            conjoint.xs("Left", level="LABEL").RESPONDENT_COUNTRY.cat.codes,
            dims="choice_situations"
        )
        dummies_left = pm.ConstantData(
            "dummies_left",
            dummies.xs("Left", level="LABEL").values,
            dims=["choice_situations", "levels"]
        )
        dummies_right = pm.ConstantData(
            "dummies_right",
            dummies.xs("Right", level="LABEL").values,
            dims=["choice_situations", "levels"]
        )
        choice_left = pm.ConstantData(
            "choice_left",
            conjoint.xs("Left", level="LABEL").CHOICE_INDICATOR.values,
            dims="choice_situations"
        )

        u_left = intercept_left[c] + pm.math.sum(partworths[c] * dummies_left, axis=1)
        u_right = pm.math.sum(partworths[c] * dummies_right, axis=1)
        p_left = pm.Deterministic(
            "p_left",
            pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right)),
            dims="choice_situations"
        )

        pm.Bernoulli(
            "choice",
            p=p_left,
            observed=choice_left,
            dims="choice_situations"
        )

        match sample_type:
            case "prior":
                inference_data = pm.sample_prior_predictive(samples=n_draws, random_seed=random_seed)
            case "posterior":
                inference_data = pm.sample(
                    draws=n_draws,
                    tune=n_tune,
                    cores=n_cores,
                    random_seed=random_seed,
                    return_inferencedata=True
                )
            case _:
                raise ValueError(f"Unknown sample type {sample_type}.")

    inference_data.to_netcdf(path_to_output)


if __name__ == "__main__":
    multinomial_logit_model(
        path_to_data=snakemake.input.data,
        sample_type=snakemake.wildcards.sample,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        n_cores=snakemake.threads,
        limit_respondents=bool(snakemake.params.limit_respondents),
        n_respondents_per_country=int(snakemake.params.limit_respondents) // 4,
        random_seed=snakemake.params.random_seed,
        path_to_output=snakemake.output[0]
    )

import numpy as np
import pandas as pd
import pymc as pm
import aesara.tensor as at


def hierarchical_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, n_respondents: int,
                       random_seed: int, path_to_output: str):
    data = pd.read_feather(path_to_data)
    pure_conjoint = (
        data
        .iloc[:16 * n_respondents, :10]
        .set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
    )
    pure_conjoint = pd.get_dummies(pure_conjoint, drop_first=True)

    left_dummies = pure_conjoint.xs("Left", level="LABEL").drop(columns="CHOICE_INDICATOR")
    right_dummies = pure_conjoint.xs("Right", level="LABEL").drop(columns="CHOICE_INDICATOR")
    choice_left = pure_conjoint.xs("Left", level="LABEL").CHOICE_INDICATOR

    model = pm.Model(coords={
        "level": left_dummies.columns.values,
        "respondent": left_dummies.index.get_level_values("RESPONDENT_ID").unique()
    })

    with model:
        alpha = pm.Normal('alpha', 0, sigma=10, dims="level")
        partworths = pm.MvNormal( # TODO N or MVN?
            "partworths",
            alpha,
            tau=np.eye(len(model.coords["level"])),
            dims=["respondent", "level"]
        )

        u_left = at.concatenate([
            pm.math.dot(
                left_dummies.xs(respondent, level="RESPONDENT_ID").values,
                partworths[i, :]
            )
            for i, respondent in enumerate(model.coords["respondent"])
        ], axis=0)
        u_right = at.concatenate([
            pm.math.dot(
                right_dummies.xs(respondent, level="RESPONDENT_ID").values,
                partworths[i, :]
            )
            for i, respondent in enumerate(model.coords["respondent"])
        ], axis=0)

        p_left = pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right))

        choices = pm.Bernoulli(
            "choice",
            p=p_left,
            observed=choice_left.values
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
    hierarchical_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        n_respondents=snakemake.params.n_respondents,
        n_cores=snakemake.threads,
        random_seed=snakemake.params.random_seed,
        path_to_output=snakemake.output[0]
    )

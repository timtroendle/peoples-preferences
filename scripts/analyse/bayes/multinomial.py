import pandas as pd
import pymc3 as pm
import arviz as az


def multinomial_logit_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, path_to_output: str):
    data = pd.read_feather(path_to_data)
    pure_conjoint = data.iloc[:, :10]
    dummies = pd.get_dummies(pure_conjoint.iloc[:, 4:], drop_first=True)
    left_dummies = dummies[pure_conjoint["LABEL"] == "Left"]
    right_dummies = dummies[pure_conjoint["LABEL"] == "Right"]
    choice_left = (
        pure_conjoint
        .set_index(["RESPONDENT_ID", "CHOICE_SET", "LABEL"])
        .to_xarray()
        .sel(LABEL="Left")
        .CHOICE_INDICATOR
        .to_series()
    )

    model = pm.Model(coords={
        "level": dummies.columns.values
    })

    with model:
        partworths = pm.Normal("partworths", 0, 10, dims="level") # TODO N or MVN?
        u_left = pm.math.dot(left_dummies.values, partworths)
        u_right = pm.math.dot(right_dummies.values, partworths)

        p_left = pm.math.exp(u_left) / (pm.math.exp(u_left) + pm.math.exp(u_right))

        choices = pm.Bernoulli("choice", p=p_left, observed=choice_left.values)

        inference_data = pm.sample(
            draws=n_draws,
            tune=n_tune,
            cores=n_cores,
            return_inferencedata=True
        )

    inference_data.to_netcdf(path_to_output)


if __name__ == "__main__":
    multinomial_logit_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        n_cores=snakemake.threads,
        path_to_output=snakemake.output[0]
    )

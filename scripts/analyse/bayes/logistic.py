import pandas as pd
import pymc3 as pm
import arviz as az


def logistic_regression_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int, path_to_output: str):
    data = pd.read_feather(path_to_data)
    pure_conjoint = data.iloc[:, :10]
    dummies = pd.get_dummies(pure_conjoint.iloc[:, 4:], drop_first=True)

    model = pm.Model(coords={
        "partsworths": dummies.columns.values
    })

    with model:
        α = pm.Normal('α', mu=0, sd=10)
        partworths = pm.Normal("partworths", 0, 10, dims="partsworths") # TODO or MVN?
        u = α + pm.math.dot(dummies.values, partworths)

        θ = pm.math.sigmoid(u)
        y_1 = pm.Bernoulli('y_1', p=θ, observed=pure_conjoint.CHOICE_INDICATOR.values)

        inference_data = pm.sample(
            draws=n_draws,
            tune=n_tune,
            cores=n_cores,
            return_inferencedata=True
        )

    inference_data.to_netcdf(path_to_output)


if __name__ == "__main__":
    logistic_regression_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        n_cores=snakemake.threads,
        path_to_output=snakemake.output[0]
    )

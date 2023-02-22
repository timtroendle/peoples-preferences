import pandas as pd
import pymc as pm


def logistic_regression_model(path_to_data: str, n_tune: int, n_draws: int, n_cores: int,
                              random_seed: int, path_to_output: str):
    data = pd.read_feather(path_to_data)
    pure_conjoint = data.iloc[:, :10]
    dummies = pd.get_dummies(pure_conjoint.iloc[:, 4:], drop_first=True)

    model = pm.Model(coords={
        "partsworths": dummies.columns.values
    })

    with model:
        α = pm.Normal('α', mu=0, sigma=10)
        partworths = pm.Normal("partworths", 0, 10, dims="partsworths") # TODO or MVN?
        u = α + pm.math.dot(dummies.values, partworths)

        θ = pm.math.sigmoid(u)
        pm.Bernoulli('y_1', p=θ, observed=pure_conjoint.CHOICE_INDICATOR.values)

        inference_data = pm.sample(
            draws=n_draws,
            tune=n_tune,
            cores=n_cores,
            random_seed=random_seed,
            return_inferencedata=True
        )

    inference_data.to_netcdf(path_to_output)


if __name__ == "__main__":
    logistic_regression_model(
        path_to_data=snakemake.input.data,
        n_tune=snakemake.params.n_tune,
        n_draws=snakemake.params.n_draws,
        n_cores=snakemake.threads,
        random_seed=snakemake.params.random_seed,
        path_to_output=snakemake.output[0]
    )

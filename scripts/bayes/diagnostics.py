from pathlib import Path

import arviz as az
import pandas as pd
import seaborn as sns


def diagnostics(path_to_inference_data: str, path_to_trace_plot: str, path_to_pop_means_plot: str,
                path_to_forest_plot: str, path_to_summary: str, path_to_rhos_individual_plot: str,
                path_to_rhos_country_plot: str,
                path_to_confusion_matrix: str, path_to_accuracy: str,
                hdi_prob: float, path_to_individuals_plot: str):
    inference_data = az.from_netcdf(path_to_inference_data)
    inference_data = retransform_normalised(inference_data)

    trace_plot(inference_data, path_to_trace_plot)
    pop_means_plot(inference_data, hdi_prob, path_to_pop_means_plot)
    forest_plot(inference_data, hdi_prob, path_to_forest_plot)
    rhos_plot(inference_data, "rho_individuals", path_to_rhos_individual_plot)
    rhos_plot(inference_data, "rho_country", path_to_rhos_country_plot)
    individuals_plot(inference_data, path_to_individuals_plot)
    summary(inference_data, hdi_prob, path_to_summary)
    prediction_accuracy(inference_data, path_to_confusion_matrix, path_to_accuracy)


def retransform_normalised(inference_data: az.InferenceData):
    age_std = inference_data.constant_data.age.std()
    years_std = inference_data.constant_data.years.std()
    if "beta_age_normed" in inference_data.posterior:
        inference_data.posterior["beta_age"] = inference_data.posterior["beta_age_normed"] / age_std
    if "beta_years_normed" in inference_data.posterior:
        inference_data.posterior["beta_years"] = inference_data.posterior["beta_years_normed"] / years_std
    return inference_data


def trace_plot(inference_data: az.InferenceData, path_to_plot: str):
    var_names = ["alpha", "beta", "sigma_individuals", "sigma_country", "mu_left_intercept", "sigma_left_intercept"]
    axes = az.plot_trace(
        inference_data,
        var_names=var_names,
        filter_vars="like"
    )
    fig = axes[0, 0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def pop_means_plot(inference_data: az.InferenceData, hdi_prob: float, path_to_plot: str):
    axes = az.plot_forest(inference_data, var_names="alpha", combined=True, hdi_prob=hdi_prob)
    ax = axes[0]
    ax.vlines(x=0, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color="black", linestyles="dotted")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def forest_plot(inference_data: az.InferenceData, hdi_prob: float, path_to_plot: str):
    var_names = ["alpha", "beta", "sigma_individuals", "mu_left_intercept", "sigma_left_intercept"]
    axes = az.plot_forest(
        inference_data,
        var_names=var_names,
        filter_vars="like",
        combined=False,
        hdi_prob=hdi_prob
    )
    ax = axes[0]
    ax.vlines(x=0, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color="black", linestyles="dotted")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def rhos_plot(inference_data: az.InferenceData, parameter: str, path_to_plot: str):
    rhos = (
        inference_data
        .posterior[parameter]
        .mean(["draw"])
        .to_dataframe()[parameter]
        .reset_index()
    )
    level_mapper = inference_data.constant_data.level.to_series().reset_index(drop=True).to_dict()
    def heatmap(df):
        return (
            df
            .pivot(index="level", columns="level_repeat", values=parameter)
            .rename(index=level_mapper, columns=level_mapper)
        )

    g = sns.FacetGrid(rhos, col='chain', col_wrap=2, height=6)
    g.map_dataframe(lambda data, color: sns.heatmap(heatmap(data), square=True, vmin=-1, vmax=1))
    g.savefig(path_to_plot)


def draw_and_chain_mean_covariates(data):
    return (
        draw_and_chain_mean_nocovariates(data)
        .sel(gender="Male") # remove "other" gender because sample is very small
    )


def draw_and_chain_mean_nocovariates(data):
    return (
        data
        .mean(["draw", "chain"])
        .expand_dims("chain")
    )


def individuals_plot(inference_data: az.InferenceData, path_to_plot: str):
    if "gender_effect" in inference_data.posterior:
        var_names=[
            "partworths", "gender_effect", "country_effect", "area_effect", "renewables_effect",
            "party_effect", "age_effect", "years_effect", "edu_effect", "income_effect", "concern_effect",
            "individuals"
        ]
        draw_and_chain_mean = draw_and_chain_mean_covariates
    else:
        var_names=["partworths", "individuals"]
        draw_and_chain_mean = draw_and_chain_mean_nocovariates
    axes = az.plot_forest(
        inference_data,
        var_names=var_names,
        combine_dims=set(["respondent"]),
        combined=True,
        transform=draw_and_chain_mean
    )
    ax = axes[0]
    ax.vlines(x=0, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color="black", linestyles="dotted")
    ax.set_title("Range of average individual-level partworths")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def summary(inference_data: az.InferenceData, hdi_prob: float, path_to_summary: str):
    (
        az
        .summary(
            inference_data,
            var_names=["alpha", "beta", "sigma_country", "sigma_individuals",
                       "mu_left_intercept", "sigma_left_intercept",
                       "d_edu", "rho_individuals", "rho_country"],
            filter_vars="like",
            hdi_prob=hdi_prob,
        )
        .reset_index()
        .to_feather(path_to_summary)
    )


def prediction_accuracy(inference_data: az.InferenceData, path_to_confusion_matrix: str, path_to_accuracy: str):
    best_estimate = (
        inference_data
        .posterior
        .p_left
        .mean(["chain", "draw"])
        .to_series()
        .map(lambda x: 1 if x >= 0.5 else 0)
    )
    choices = inference_data.observed_data.choice.to_series()
    confusion_matrix = pd.crosstab(choices, best_estimate)
    accuracy = (confusion_matrix.loc[0, 0] + confusion_matrix.loc[1, 1]) / confusion_matrix.sum().sum()

    confusion_matrix.to_csv(path_to_confusion_matrix)
    with Path(path_to_accuracy).open("w") as f_acc:
        f_acc.write(f"{accuracy:.2f}")


if __name__ == "__main__":
    diagnostics(
        path_to_inference_data=snakemake.input.inference_data,
        path_to_trace_plot=snakemake.output.trace,
        path_to_pop_means_plot=snakemake.output.pop_means,
        path_to_forest_plot=snakemake.output.forest,
        path_to_rhos_individual_plot=snakemake.output.rhos_individual,
        path_to_rhos_country_plot=snakemake.output.rhos_country,
        path_to_individuals_plot=snakemake.output.individuals,
        path_to_summary=snakemake.output.summary,
        path_to_confusion_matrix=snakemake.output.confusion,
        path_to_accuracy=snakemake.output.accuracy,
        hdi_prob=snakemake.params.hdi_prob
    )

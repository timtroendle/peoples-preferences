from functools import cache
from pathlib import Path

import arviz as az
import pandas as pd
import xarray as xr
import seaborn as sns


def diagnostics(path_to_inference_data: str, sample_type: str,
                path_to_trace_plot: str, path_to_pop_means_plot: str,
                path_to_forest_plot: str, path_to_summary: str, path_to_rhos_individual_plot: str,
                path_to_rhos_country_plot: str, path_to_choice_probability_plot: str,
                path_to_confusion_matrix: str, path_to_accuracy: str, path_to_utility_plot: str,
                hdi_prob: float, path_to_individuals_plot: str):
    inference_data = az.from_netcdf(path_to_inference_data)
    observed_data = inference_data.observed_data.load() # This avoids a weird segmentation fault on my machine.
                                                        # It seems this segmentation fault is caused by too many
                                                        # accesses to the data. Therefore, I am loading it into
                                                        # memory here.
    constant_data = inference_data.constant_data
    inference_data = inference_data[sample_type]
    inference_data = retransform_normalised(inference_data, constant_data)

    trace_plot(inference_data, path_to_trace_plot)
    pop_means_plot(inference_data, hdi_prob, path_to_pop_means_plot)
    forest_plot(inference_data, hdi_prob, path_to_forest_plot)
    rhos_plot(inference_data, rho="rho_individuals", var="individuals", path_to_plot=path_to_rhos_individual_plot)
    rhos_plot(inference_data, rho="rho_country", var="countries", path_to_plot=path_to_rhos_country_plot)
    individuals_plot(inference_data, path_to_individuals_plot)
    summary(inference_data, hdi_prob, path_to_summary)
    prediction_accuracy(inference_data, observed_data, path_to_confusion_matrix, path_to_accuracy)
    choice_probability_plot(inference_data, hdi_prob, path_to_choice_probability_plot)
    utility_plot(inference_data, hdi_prob, path_to_utility_plot)


def retransform_normalised(inference_data: xr.Dataset, constant_data: xr.Dataset):
    if "beta_age_normed" in inference_data:
        age_std = constant_data.age.std()
        inference_data["beta_age"] = inference_data["beta_age_normed"] / age_std
    if "beta_years_normed" in inference_data:
        years_std = inference_data.constant_data.years.std()
        inference_data["beta_years"] = inference_data["beta_years_normed"] / years_std
    return inference_data


def trace_plot(inference_data: xr.Dataset, path_to_plot: str):
    var_names = ["alpha", "beta", "sigma_individuals", "sigma_country", "mu_left_intercept", "sigma_left_intercept"]
    axes = az.plot_trace(
        inference_data,
        var_names=var_names,
        filter_vars="like"
    )
    fig = axes[0, 0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def pop_means_plot(inference_data: xr.Dataset, hdi_prob: float, path_to_plot: str):
    axes = az.plot_forest(inference_data, var_names="alpha", combined=True, hdi_prob=hdi_prob)
    ax = axes[0]
    ax.vlines(x=0, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], color="black", linestyles="dotted")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def forest_plot(inference_data: xr.Dataset, hdi_prob: float, path_to_plot: str):
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


def rhos_plot(inference_data: xr.Dataset, rho: str, var: str, path_to_plot: str):
    level_mapper = inference_data.level.to_series().reset_index(drop=True).to_dict()
    if rho in inference_data:
        data = (
            inference_data[rho]
            .mean(["draw"])
            .to_dataframe()[rho]
            .reset_index()
        )

        def heatmap(df):
            return (
                df
                .pivot(index="level", columns="level_repeat", values=rho)
                .rename(index=level_mapper, columns=level_mapper)
            )
    else:
        # calculate rhos yourself
        data = (
            inference_data[var]
        )
        data = pd.concat([correlation_across_level(data, chain) for chain in data.chain])

        def heatmap(df: pd.DataFrame):
            return (
                df["correlation"]
                .unstack()
            )

    g = sns.FacetGrid(data, col='chain', col_wrap=2, height=6)
    g.map_dataframe(lambda data, color: sns.heatmap(heatmap(data), square=True, vmin=-1, vmax=1, cmap="vlag"))
    g.savefig(path_to_plot)


def correlation_across_level(data: xr.DataArray, chain: int) -> pd.DataFrame:
    levels = [level.item() for level in data.level]
    data = data.sel(chain=chain)

    @cache
    def correlation(level1: str, level2: str) -> float:
        if level1 == level2:
            return 1.0
        else:
            return xr.corr(data.sel(level=level1), data.sel(level=level2)).item()

    return (
        pd
        .DataFrame(
            index=levels,
            data={
                level1: [correlation(*sorted([level1, level2])) for level2 in levels]
                for level1 in levels
            }
        )
        .stack() # from heatmap-form to long-form to be usable within Seaborn's FacetGrid
        .to_frame("correlation")
        .assign(chain=chain.item())
    )


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


def individuals_plot(inference_data: xr.Dataset, path_to_plot: str):
    if "gender_effect" in inference_data:
        var_names = [
            "partworths", "gender_effect", "country_effect", "area_effect", "renewables_effect",
            "party_effect", "age_effect", "years_effect", "edu_effect", "income_effect", "concern_effect",
            "individuals"
        ]
        draw_and_chain_mean = draw_and_chain_mean_covariates
    else:
        var_names = ["partworths", "individuals"]
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


def summary(inference_data: xr.Dataset, hdi_prob: float, path_to_summary: str):
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


def prediction_accuracy(inference_data: xr.Dataset, observed_data: xr.Dataset, path_to_confusion_matrix: str,
                        path_to_accuracy: str):
    best_estimate = (
        inference_data
        .p_left
        .mean(["chain", "draw"])
        .to_series()
        .map(lambda x: 1 if x >= 0.5 else 0)
    )
    choices = observed_data.choice.to_series()
    confusion_matrix = pd.crosstab(choices, best_estimate)
    accuracy = (confusion_matrix.loc[0, 0] + confusion_matrix.loc[1, 1]) / confusion_matrix.sum().sum()

    confusion_matrix.to_csv(path_to_confusion_matrix)
    with Path(path_to_accuracy).open("w") as f_acc:
        f_acc.write(f"{accuracy:.2f}")


def choice_probability_plot(inference_data: xr.Dataset, hdi_prob: float, path_to_plot: str):
    axes = az.plot_density(
        inference_data,
        var_names=["p_left"],
        combine_dims=set(["chain", "draw", "choice_situation"]),
        hdi_prob=hdi_prob
    )
    fig = axes[0, 0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def utility_plot(inference_data: xr.Dataset, hdi_prob: float, path_to_plot: str):
    axes = az.plot_density(
        inference_data,
        var_names=["u_left", "u_right"],
        combine_dims=set(["chain", "draw", "choice_situation"]),
        hdi_prob=hdi_prob
    )
    fig = axes[0, 0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    diagnostics(
        path_to_inference_data=snakemake.input.inference_data,
        sample_type=snakemake.wildcards.sample,
        path_to_trace_plot=snakemake.output.trace,
        path_to_pop_means_plot=snakemake.output.pop_means,
        path_to_forest_plot=snakemake.output.forest,
        path_to_rhos_individual_plot=snakemake.output.rhos_individual,
        path_to_rhos_country_plot=snakemake.output.rhos_country,
        path_to_individuals_plot=snakemake.output.individuals,
        path_to_summary=snakemake.output.summary,
        path_to_confusion_matrix=snakemake.output.confusion,
        path_to_accuracy=snakemake.output.accuracy,
        path_to_choice_probability_plot=snakemake.output.probability,
        path_to_utility_plot=snakemake.output.utility,
        hdi_prob=snakemake.params.hdi_prob
    )

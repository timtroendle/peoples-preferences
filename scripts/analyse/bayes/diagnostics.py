import arviz as az
import seaborn as sns


def diagnostics(path_to_inference_data: str, path_to_trace_plot: str, path_to_pop_means_plot: str,
                path_to_forest_plot: str, path_to_summary: str, path_to_rhos_plot: str,
                hdi_prob: float, path_to_individuals_plot: str):
    inference_data = az.from_netcdf(path_to_inference_data)
    inference_data = retransform_normalised(inference_data)

    trace_plot(inference_data, path_to_trace_plot)
    pop_means_plot(inference_data, hdi_prob, path_to_pop_means_plot)
    forest_plot(inference_data, hdi_prob, path_to_forest_plot)
    rhos_plot(inference_data, path_to_rhos_plot)
    individuals_plot(inference_data, path_to_individuals_plot)
    summary(inference_data, hdi_prob, path_to_summary)


def retransform_normalised(inference_data: az.InferenceData):
    age_std = inference_data.constant_data.age.std()
    years_std = inference_data.constant_data.years.std()
    if "beta_age_normed" in inference_data.posterior:
        inference_data.posterior["beta_age"] = inference_data.posterior["beta_age_normed"] / age_std
    if "beta_years_normed" in inference_data.posterior:
        inference_data.posterior["beta_years"] = inference_data.posterior["beta_years_normed"] / years_std
    return inference_data


def trace_plot(inference_data: az.InferenceData, path_to_plot: str):
    var_names = ["alpha", "beta", "sigma_individuals", "mu_left_intercept", "sigma_left_intercept"]
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
    fig = axes[0].get_figure()
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
    fig = axes[0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def rhos_plot(inference_data: az.InferenceData, path_to_plot: str):
    rhos = (
        inference_data
        .posterior
        .rho_individuals
        .mean(["draw"])
        .to_dataframe()["rho_individuals"]
        .reset_index()
    )
    level_mapper = inference_data.constant_data.level.to_series().reset_index(drop=True).to_dict()
    def heatmap(df):
        return (
            df
            .pivot(index="level", columns="level_repeat", values="rho_individuals")
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
    axes[0].set_title("Range of average individual-level partworths")
    fig = axes[0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def summary(inference_data: az.InferenceData, hdi_prob: float, path_to_summary: str):
    (
        az
        .summary(
            inference_data,
            var_names=["alpha", "beta", "sigma_individuals", "mu_left_intercept", "sigma_left_intercept",
                       "d_edu", "rho_individuals"],
            filter_vars="like",
            hdi_prob=hdi_prob,
        )
        .reset_index()
        .to_feather(path_to_summary)
    )


if __name__ == "__main__":
    diagnostics(
        path_to_inference_data=snakemake.input.inference_data,
        path_to_trace_plot=snakemake.output.trace,
        path_to_pop_means_plot=snakemake.output.pop_means,
        path_to_forest_plot=snakemake.output.forest,
        path_to_rhos_plot=snakemake.output.rhos,
        path_to_individuals_plot=snakemake.output.individuals,
        path_to_summary=snakemake.output.summary,
        hdi_prob=snakemake.params.hdi_prob
    )

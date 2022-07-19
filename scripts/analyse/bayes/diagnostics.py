import arviz as az
import seaborn as sns


def diagnostics(path_to_inference_data: str, path_to_trace_plot: str, path_to_pop_means_plot: str, path_to_forest_plot: str, path_to_summary: str, path_to_rhos_plot: str):
    inference_data = az.from_netcdf(path_to_inference_data)

    trace_plot(inference_data, path_to_trace_plot)
    pop_means_plot(inference_data, path_to_pop_means_plot)
    forest_plot(inference_data, path_to_forest_plot)
    rhos_plot(inference_data, path_to_rhos_plot)
    summary(inference_data, path_to_summary)


def trace_plot(inference_data, path_to_plot):
    axes = az.plot_trace(
        inference_data,
        var_names=["alpha", "sigma_partworths", "mu_left_intercept", "sigma_left_intercept"]
    )
    fig = axes[0, 0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def pop_means_plot(inference_data, path_to_plot):
    axes = az.plot_forest(inference_data, var_names="alpha", combined=True, hdi_prob=0.9)
    fig = axes[0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def forest_plot(inference_data, path_to_plot):
    axes = az.plot_forest(
        inference_data,
        var_names=["alpha", "sigma_partworths", "mu_left_intercept", "sigma_left_intercept", "Rho_partworths"],
        combined=False,
        hdi_prob=0.9
    )
    fig = axes[0].get_figure()
    fig.tight_layout()
    fig.savefig(path_to_plot)


def rhos_plot(inference_data, path_to_plot):
    rhos = (
        inference_data
        .posterior
        .Rho_partworths
        .mean(["draw"])
        .to_dataframe()["Rho_partworths"]
        .reset_index()
        .rename(columns={"chol_partworths_corr_dim_0": "dim0", "chol_partworths_corr_dim_1": "dim1"})
    )
    level_mapper = inference_data.constant_data.level.to_series().reset_index(drop=True).to_dict()
    def heatmap(df):
        return (
            df
            .pivot(index="dim0", columns="dim1", values="Rho_partworths")
            .rename(index=level_mapper, columns=level_mapper)
        )

    g = sns.FacetGrid(rhos, col='chain', col_wrap=2, height=6)
    g.map_dataframe(lambda data, color: sns.heatmap(heatmap(data), square=True, vmin=-1, vmax=1))
    g.savefig(path_to_plot)


def summary(inference_data, path_to_summary):
    (
        az
        .summary(
            inference_data,
            var_names=["alpha", "sigma_partworths", "mu_left_intercept", "sigma_left_intercept", "Rho_partworths"],
            hdi_prob=0.9,
            round_to=2
        )
        .to_csv(path_to_summary)
    )


if __name__ == "__main__":
    diagnostics(
        path_to_inference_data=snakemake.input.inference_data,
        path_to_trace_plot=snakemake.output.trace,
        path_to_pop_means_plot=snakemake.output.pop_means,
        path_to_forest_plot=snakemake.output.forest,
        path_to_rhos_plot=snakemake.output.rhos,
        path_to_summary=snakemake.output.summary
    )

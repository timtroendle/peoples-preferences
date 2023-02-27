import arviz as az
import altair as alt
import pandas as pd
import xarray as xr


DARK_GREY = "#424242"


def plot_varying_effects(path_to_inference_data: str, varying_variable_name: str, pop_mean_variable_name: str,
                         sample_type: str,
                         level_choice: dict[str, str], narrow_hdi: float, wide_hdi: str, path_to_plot: str):
    data = az.from_netcdf(path_to_inference_data)[sample_type]
    varying = uncertainty_areas(
        varying=data[varying_variable_name].sel(level_choice, drop=True),
        wide_hdi=wide_hdi,
        narrow_hdi=narrow_hdi
    )
    pop_means = (
        data[pop_mean_variable_name]
        .sel(level_choice, drop=True)
        .mean("chain")
        .to_dataframe()
    )
    varying_effects_chart = create_chart_of_varying_effects(varying)
    pop_means_chart = create_chart_of_pop_means(pop_means, pop_mean_variable_name)

    chart = alt.hconcat(
        varying_effects_chart,
        pop_means_chart
    ).resolve_scale(
        y='shared'
    )

    (
        chart
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )


def create_chart_of_varying_effects(varying: pd.DataFrame):
    base = (
        alt
        .Chart(
            varying.reset_index(),
            title="A"
        ).encode(
            x=alt.X("respondent:N", title="Respondent", axis=alt.Axis(labels=False, ticks=False), sort=alt.EncodingSortField(field="mean")),
            y=alt.Y("mean:Q", title="Individual partworths")
        ).properties(
            width=250,
            height=200
        )
    )
    interval_narrow = (
        base
        .mark_area(opacity=0.5)
        .encode(y="lower", y2="higher")
    )
    interval_wide = (
        base
        .mark_area(opacity=0.4)
        .encode(y="lowest", y2="highest")
    )
    line = base.mark_line(size=1)
    return interval_wide + interval_narrow + line


def create_chart_of_pop_means(pop_means: pd.DataFrame, variable_name: str):
    return (
        alt
        .Chart(
            pop_means.reset_index(),
            title="B"
        ).transform_density(
            variable_name,
            as_=[variable_name, 'density']
        ).encode(
            y=alt.Y(f"{variable_name}:Q", title="Population-mean partworth"),
            x=alt.X("density:Q", title="Density", sort="ascending")
        ).properties(
            width=50,
            height=200
        ).mark_line(
            orient=alt.Orientation("horizontal")
        )
    )


def uncertainty_areas(varying: xr.DataArray, wide_hdi: float, narrow_hdi: float):
    varying_wide_hdi = az.hdi(varying, wide_hdi).to_array().isel(variable=0, drop=True)
    varying_narrow_hdi = az.hdi(varying, narrow_hdi).to_array().isel(variable=0, drop=True)
    return pd.DataFrame({
        "mean": varying.mean(["chain", "draw"]).to_series(),
        "higher": varying_narrow_hdi.sel(hdi="higher", drop=True),
        "lower": varying_narrow_hdi.sel(hdi="lower", drop=True),
        "highest": varying_wide_hdi.sel(hdi="higher", drop=True),
        "lowest": varying_wide_hdi.sel(hdi="lower", drop=True),
    })


if __name__ == "__main__":
    plot_varying_effects(
        path_to_inference_data=snakemake.input.data,
        sample_type=snakemake.wildcards.sample,
        varying_variable_name=snakemake.params.varying_variable_name,
        pop_mean_variable_name=snakemake.params.pop_means_variable_name,
        level_choice=snakemake.params.level_choice,
        narrow_hdi=snakemake.params.narrow_hdi,
        wide_hdi=snakemake.params.wide_hdi,
        path_to_plot=snakemake.output[0]
    )

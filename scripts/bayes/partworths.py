from typing import Union

import pandas as pd
import xarray as xr
import arviz as az

import altair as alt


DARK_GREY = "#424242"


BASELINE_LEVELS = [ # TODO ADD FROM CONFIG
    "TECHNOLOGY:Rooftop PV",
    "LAND:0.5%",
    "PRICES:+0%",
    "TRANSMISSION:-25.0% .",
    "OWNERSHIP:Public",
    "SHARE_IMPORTS:0%"
]

EXPECTED_VALUE = "expected"
ESTIMAND_NAME = {
    "prior": "Partworth utility",
    "posterior": "Partworth utility",
    "amce": "Average marginal component effect",
    "poststratify": "Partworth utility",
    "subgroups": "Partworth utility"
}
WIDTH_SINGLE_PLOT = 400 - 7
WIDTH_SINGLE_PLOT_IN_FACETS = 96.5


def visualise_partworths(path_to_inference_data: str, facet_by_country: bool, aggregate_individuals: bool,
                         variable_names: Union[str, list[str]], hdi_prob: float, sample_type: str,
                         nice_names: dict[str, dict[str, str]], estimand_name: str,
                         path_to_plot: str):
    data = read_data(
        path_to_inference_data=path_to_inference_data,
        sample_type=sample_type,
        variable_names=variable_names,
        facet_by_country=facet_by_country,
        aggregate_individuals=aggregate_individuals,
        hdi_prob=hdi_prob,
        nice_names=nice_names
    )

    base = alt.Chart(
        data
    ).encode(
        y=alt.Y("level", sort=list(nice_names["levels"].values()), title="Level"),
        x=alt.X(EXPECTED_VALUE, title=estimand_name),
        color=alt.Color("attribute", sort=list(nice_names["attributes"].values()), legend=alt.Legend(title="Attribute")),
    ).properties(
        width=WIDTH_SINGLE_PLOT if not facet_by_country else WIDTH_SINGLE_PLOT_IN_FACETS,
    )

    base_line = alt.Chart(data).mark_rule(color=DARK_GREY, strokeDash=[4], opacity=0.8).encode(
        x='zero:Q'
    )

    area = base.mark_area(opacity=1, filled=True, fillOpacity=0.7, line=False).encode(
        x=alt.value(0),
        x2=alt.value(4),
    )

    if aggregate_individuals:
        main_marks = base.mark_boxplot(opacity=1, size=10)
    else:
        point = base.mark_circle(opacity=1)
        interval = base.mark_rule(strokeWidth=1.5, opacity=0.6).encode(
            x="lower",
            x2="higher"
        )
        main_marks = (interval + point)

    entire_chart = (area + base_line + main_marks)

    if facet_by_country:
        entire_chart = entire_chart.facet(alt.Column('Country', title=None), columns=4)
    (
        entire_chart
        .configure(font="Lato", facet=alt.CompositionConfig(spacing=0))
        .configure_view(step=16)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom", columns=3)
        .save(path_to_plot)
    )


def read_data(path_to_inference_data: str, variable_names: Union[str, list[str]], facet_by_country: bool,
              sample_type: str,
              aggregate_individuals: True, nice_names: dict[str, dict[str, str]], hdi_prob: float):
    full = az.from_netcdf(path_to_inference_data)[sample_type]
    attr_levels = full.level.to_series().to_list()
    all_attr_levels = attr_levels + BASELINE_LEVELS
    fill_value = pd.NA if aggregate_individuals else 0

    data = (
        full[variable_names]
        .pipe(sum_variables_if_more_than_one)
        .reindex(level=all_attr_levels)
        .fillna(fill_value)
    )

    expected_value = (
        data
        .mean(["chain", "draw"])
        .rename(EXPECTED_VALUE)
        .to_series()
        .reset_index()
    )
    expected_value["attribute"] = expected_value.level.str.split(":").str[0].map(nice_names["attributes"])
    expected_value["level"] = expected_value.level.str.split(":").str[1].map(nice_names["levels"])

    if not aggregate_individuals:
        # interval is highest density interval
        interval = (
            az
            .hdi(data, hdi_prob=hdi_prob)
            .to_array()
            .isel(variable=0)
            .to_series()
            .unstack("hdi")
            .reset_index()
        )
        interval["attribute"] = interval.level.str.split(":").str[0].map(nice_names["attributes"])
        interval["level"] = interval.level.str.split(":").str[1].map(nice_names["levels"])

        index_cols = ["attribute", "level", "country"] if facet_by_country else ["attribute", "level"]
        data = pd.merge(left=interval, right=expected_value, on=index_cols, validate="1:1")
    else:
        data = expected_value

    data = (
        data
        .pipe(preprocess_country_if_necessary, facet_by_country, nice_names["countries"])
        .assign(zero=0)
    )
    return data


def sum_variables_if_more_than_one(data: Union[xr.DataArray, xr.Dataset]) -> xr.DataArray:
    try:
        return data.to_array().sum("variable")
    except AttributeError:
        # not a Dataset
        return data


def preprocess_country_if_necessary(df: pd.DataFrame, facet_by_country: bool, nice_country_names):
    if facet_by_country:
        return df.assign(Country=df.country.map(nice_country_names))
    else:
        return df


def optional_param(name: str, default):
    return snakemake.params[name] if name in snakemake.params.keys() else default


if __name__ == "__main__":
    visualise_partworths(
        path_to_inference_data=snakemake.input.data,
        sample_type=snakemake.wildcards.sample,
        nice_names=snakemake.params.nice_names,
        estimand_name=ESTIMAND_NAME[snakemake.wildcards.sample],
        facet_by_country=optional_param("facet_by_country", False),
        aggregate_individuals=optional_param("aggregate_individuals", False),
        variable_names=snakemake.params.variable_names,
        hdi_prob=snakemake.params.hdi_prob,
        path_to_plot=snakemake.output[0]
    )

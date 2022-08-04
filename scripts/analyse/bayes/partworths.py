import pandas as pd
import arviz as az

import altair as alt


DARK_GREY = "#424242"


BASELINE_LEVELS = [ # TODO ADD FROM CONFIG
    "TECHNOLOGY:Rooftop PV",
    "LAND:0.5%",
    "PRICES:+0%",
    "TRANSMISSION:-25% .",
    "OWNERSHIP:Public",
    "SHARE_IMPORTS:0%"
]


def visualise_partworths(path_to_posterior: str, facet_by_country: bool, aggregate_individuals: bool,
                         variable_name: str, hdi_prob: float, title: str, nice_names: dict[str, dict[str, str]],
                         path_to_plot: str):
    data = read_data(
        path_to_posterior=path_to_posterior,
        variable_name=variable_name,
        facet_by_country=facet_by_country,
        aggregate_individuals=aggregate_individuals,
        hdi_prob=hdi_prob,
        nice_names=nice_names
    )

    base = alt.Chart(
        data,
        title=title,
    ).encode(
        y=alt.Y("level", sort=list(nice_names["levels"].values()), title="Level"),
        x=alt.X(variable_name, title="Partworth utility"),
        color=alt.Color("attribute", sort=list(nice_names["attributes"].values()), legend=alt.Legend(title="Attribute")),
    ).properties(
        width=300,
        height=300
    )

    interval = base.mark_rule(strokeWidth=1.5, opacity=0.6).encode(
        x="lower",
        x2="higher"
    )

    point = base.mark_circle(opacity=1)

    base_line = alt.Chart(data).mark_rule(color=DARK_GREY, strokeDash=[4], opacity=0.8).encode(
        x='zero:Q'
    )

    area = base.mark_area(opacity=1, filled=True, fillOpacity=0.7, line=False).encode(
        x=alt.value(0),
        x2=alt.value(4),
    )

    entire_chart = (area + base_line + interval + point)
    if facet_by_country:
        entire_chart = entire_chart.facet("Country", columns=2, title=title)
    (
        entire_chart
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )


def read_data(path_to_posterior: str, variable_name: str, facet_by_country: bool, aggregate_individuals: True,
              nice_names: dict[str, dict[str, str]], hdi_prob: float):
    full = az.from_netcdf(path_to_posterior)
    attr_levels = full.posterior.level.to_series().to_list()
    all_attr_levels = attr_levels + BASELINE_LEVELS

    if aggregate_individuals:
        # interval is mean across individuals
        individual_partworths = (
            full
            .posterior
            .mean(["draw", "chain"])
            .data_vars[variable_name]
            .reindex(level=all_attr_levels)
            .fillna(0)
        )
        interval = pd.DataFrame({
            "lower": individual_partworths.min("respondent").to_series(),
            "higher": individual_partworths.max("respondent").to_series()
        }).reset_index()
    else:
        # interval is highest density interval
        interval = (
            az
            .hdi(full.posterior.reindex(level=all_attr_levels), hdi_prob=hdi_prob)
            .data_vars[variable_name]
            .to_series()
            .unstack("hdi")
            .reset_index()
            .fillna(0)
        )
    interval["attribute"] = interval.level.str.split(":").str[0].map(nice_names["attributes"])
    interval["level"] = interval.level.str.split(":").str[1].map(nice_names["levels"])

    mean_vars = ["chain", "draw", "respondent"] if aggregate_individuals else ["chain", "draw"]

    mean = (
        full
        .posterior
        .data_vars[variable_name]
        .reindex(level=all_attr_levels)
        .fillna(0)
        .mean(mean_vars)
        .to_series()
        .reset_index()
    )
    mean["attribute"] = mean.level.str.split(":").str[0].map(nice_names["attributes"])
    mean["level"] = mean.level.str.split(":").str[1].map(nice_names["levels"])

    index_cols = ["attribute", "level", "country"] if facet_by_country else ["attribute", "level"]

    data = (
        pd
        .merge(left=interval, right=mean, on=index_cols, validate="1:1")
        .pipe(preprocess_country_if_necessary, facet_by_country, nice_names["countries"])
        .assign(zero=0)
    )
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
        path_to_posterior=snakemake.input.posterior,
        title=snakemake.params.title,
        nice_names=snakemake.params.nice_names,
        facet_by_country=optional_param("facet_by_country", False),
        aggregate_individuals=optional_param("aggregate_individuals", False),
        variable_name=snakemake.params.variable_name,
        hdi_prob=snakemake.params.hdi_prob,
        path_to_plot=snakemake.output[0]
    )

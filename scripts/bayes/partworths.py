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
                         variable_name: str, hdi_prob: float, nice_names: dict[str, dict[str, str]],
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
        data
    ).encode(
        y=alt.Y("level", sort=list(nice_names["levels"].values()), title="Level"),
        x=alt.X(variable_name, title="Partworth utility"),
        color=alt.Color("attribute", sort=list(nice_names["attributes"].values()), legend=alt.Legend(title="Attribute")),
    ).properties(
        width=300,
        height=400
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
        main_marks = (point + interval)

    entire_chart = (area + base_line + main_marks)

    if facet_by_country:
        entire_chart = entire_chart.facet("Country", columns=2)
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
    fill_value = pd.NA if aggregate_individuals else 0

    expected_value = (
        full
        .posterior
        .data_vars[variable_name]
        .reindex(level=all_attr_levels)
        .fillna(fill_value)
        .mean(["chain", "draw"])
        .to_series()
        .reset_index()
    )
    expected_value["attribute"] = expected_value.level.str.split(":").str[0].map(nice_names["attributes"])
    expected_value["level"] = expected_value.level.str.split(":").str[1].map(nice_names["levels"])

    if not aggregate_individuals:
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
        nice_names=snakemake.params.nice_names,
        facet_by_country=optional_param("facet_by_country", False),
        aggregate_individuals=optional_param("aggregate_individuals", False),
        variable_name=snakemake.params.variable_name,
        hdi_prob=snakemake.params.hdi_prob,
        path_to_plot=snakemake.output[0]
    )

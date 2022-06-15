import pandas as pd
import arviz as az

import altair as alt


DARK_GREY = "#424242"

ATTRIBUTE_ORDER = [
    "Technology",
    "Land",
    "Transmission",
    "Share of imports",
    "Prices",
    "Ownership"
]

LEVEL_ORDER = [
    "Rooftop PV",
    "Open-field PV",
    "Wind",
    "0.5%",
    "1%",
    "2%",
    "4%",
    "8%",
    '-25% .',
    "+0% .",
    '+25% .',
    '+50% .',
    '+75% .',
    "0%",
    '10%',
    '50%',
    '90%',
    "+0%",
    "+15%",
    "+30%",
    "+45%",
    "+60%",
    "Public",
    "Community",
    "Private"
]

BASELINE_LEVELS = [ # TODO ADD FROM CONFIG
    "TECHNOLOGY:Rooftop PV",
    "LAND:0.5%",
    "PRICES:+0%",
    "TRANSMISSION:-25% .",
    "OWNERSHIP:Public",
    "SHARE_IMPORTS:0%"
]

NICE_NAME_ATTRIBUTES = {
    "TECHNOLOGY": "Technology",
    "LAND": "Land",
    "PRICES": "Prices",
    "TRANSMISSION": "Transmission",
    "OWNERSHIP": "Ownership",
    "SHARE_IMPORTS": "Share of imports"
}

NICE_NAME_COUNTRIES = {
    "DEU": "Germany",
    "DNK": "Denmark",
    "POL": "Poland",
    "PRT": "Portugal"
}


def visualise_partworths(path_to_posterior: str, path_to_plot: str):
    data = read_data(path_to_posterior)

    base = alt.Chart(data).encode(
        y=alt.Y("level", sort=LEVEL_ORDER, title="Level"),
        x=alt.X("partworths", title="Partworth utility"),
        color=alt.Color("attribute", sort=ATTRIBUTE_ORDER, legend=alt.Legend(title="Attribute")),
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

    (
        (area + base_line + interval + point)
        .facet(facet="Country", columns=2)
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )


def read_data(path_to_posterior: str):
    full = az.from_netcdf(path_to_posterior)
    attr_levels = full.posterior.level.to_series().to_list()
    all_attr_levels = attr_levels + BASELINE_LEVELS
    all_attr_levels

    hdi = (
        az
        .hdi(full.posterior.reindex(level=all_attr_levels), hdi_prob=0.94)
        .partworths
        .to_series()
        .unstack("hdi")
        .reset_index()
        .fillna(0)
    )
    hdi["attribute"] = hdi.level.str.split(":").str[0].map(NICE_NAME_ATTRIBUTES)
    hdi["level"] = hdi.level.str.split(":").str[1]

    mean = (
        full
        .posterior
        .partworths
        .reindex(level=all_attr_levels)
        .fillna(0)
        .mean(["chain", "draw"])
        .to_series()
        .reset_index()
    )
    mean["attribute"] = mean.level.str.split(":").str[0].map(NICE_NAME_ATTRIBUTES)
    mean["level"] = mean.level.str.split(":").str[1]

    data = (
        pd
        .merge(left=hdi, right=mean, on=["attribute", "level", "country"], validate="1:1")
        .rename(columns={"country": "Country"}) # TODO ensure upper-case in a more smart way
        .assign(zero=0)
    )
    data.Country = data.Country.map(NICE_NAME_COUNTRIES)
    return data


if __name__ == "__main__":
    visualise_partworths(
        path_to_posterior=snakemake.input.posterior,
        path_to_plot=snakemake.output[0]
    )

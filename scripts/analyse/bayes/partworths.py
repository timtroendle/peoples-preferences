import pandas as pd
import arviz as az

import altair as alt


DARK_GREY = "#424242"

NICE_NAME_LEVELS = {
    "Rooftop PV": "Rooftop PV",
    "Open-field PV": "Open-field PV",
    "Wind": "Wind turbines",
    "0.5%": "Very low (0.5%)",
    "1%": "Low (1%)",
    "2%": "Medium (2%)",
    "4%": "High (4%)",
    "8%": "Very high (8%)",
    '-25% .': "Slight decrease (-25%)",
    "+0% .": "Today's level (0%)",
    '+25% .': "Slight increase (+25%)",
    '+50% .': "Moderate increase (+50%)",
    '+75% .': "Strong increase (+75%)",
    "0%": "None",
    '10%': "Low (10%)",
    '50%': "Medium (50%)",
    '90%': "High (90%)",
    "+0%": "Today's level",
    "+15%": "Slight increase (+15%)",
    "+30%": "Moderate increase (+30%)",
    "+45%": "Strong increase (+45%)",
    "+60%": "Very strong increase (+60%)",
    "Public": "Public",
    "Community": "Community",
    "Private": "Private"
}

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
    "TRANSMISSION": "Transmission",
    "SHARE_IMPORTS": "Share of imports",
    "PRICES": "Prices",
    "OWNERSHIP": "Ownership",
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
        y=alt.Y("level", sort=list(NICE_NAME_LEVELS.values()), title="Level"),
        x=alt.X("partworths", title="Partworth utility"),
        color=alt.Color("attribute", sort=list(NICE_NAME_ATTRIBUTES.values()), legend=alt.Legend(title="Attribute")),
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
    hdi["level"] = hdi.level.str.split(":").str[1].map(NICE_NAME_LEVELS)

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
    mean["level"] = mean.level.str.split(":").str[1].map(NICE_NAME_LEVELS)

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

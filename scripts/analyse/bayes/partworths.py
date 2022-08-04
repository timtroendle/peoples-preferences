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


def visualise_partworths(path_to_posterior: str, facet_by_country: bool, aggregate_individuals: bool,
                         variable_name: str, hdi_prob: float, title: str, path_to_plot: str):
    data = read_data(
        path_to_posterior=path_to_posterior,
        variable_name=variable_name,
        facet_by_country=facet_by_country,
        aggregate_individuals=aggregate_individuals,
        hdi_prob=hdi_prob
    )

    base = alt.Chart(
        data,
        title=title,
    ).encode(
        y=alt.Y("level", sort=list(NICE_NAME_LEVELS.values()), title="Level"),
        x=alt.X(variable_name, title="Partworth utility"),
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
              hdi_prob: float):
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
    interval["attribute"] = interval.level.str.split(":").str[0].map(NICE_NAME_ATTRIBUTES)
    interval["level"] = interval.level.str.split(":").str[1].map(NICE_NAME_LEVELS)

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
    mean["attribute"] = mean.level.str.split(":").str[0].map(NICE_NAME_ATTRIBUTES)
    mean["level"] = mean.level.str.split(":").str[1].map(NICE_NAME_LEVELS)

    index_cols = ["attribute", "level", "country"] if facet_by_country else ["attribute", "level"]

    data = (
        pd
        .merge(left=interval, right=mean, on=index_cols, validate="1:1")
        .pipe(preprocess_country_if_necessary, facet_by_country)
        .assign(zero=0)
    )
    return data


def preprocess_country_if_necessary(df: pd.DataFrame, facet_by_country: bool):
    if facet_by_country:
        return df.assign(Country=df.country.map(NICE_NAME_COUNTRIES))
    else:
        return df


def optional_param(name: str, default):
    return snakemake.params[name] if name in snakemake.params.keys() else default

if __name__ == "__main__":
    visualise_partworths(
        path_to_posterior=snakemake.input.posterior,
        title=snakemake.params.title,
        facet_by_country=optional_param("facet_by_country", False),
        aggregate_individuals=optional_param("aggregate_individuals", False),
        variable_name=snakemake.params.variable_name,
        hdi_prob=snakemake.params.hdi_prob,
        path_to_plot=snakemake.output[0]
    )

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


def visualise_covariates(path_to_posterior: str, interval: float, nice_names: dict[str, dict[str, str]],
                         path_to_plot: str):
    shares = read_shares(path_to_posterior, interval, nice_names)
    bars = alt.Chart(
        data=shares,
        title="Heterogeneity of individuals explained by covariates"
    ).mark_bar().encode(
        x=alt.X(
            'share',
            title="Covariate-induced variation across individuals relative to total variation",
            scale=alt.Scale(domain=(-0.1, 2)),
            axis=alt.Axis(format='.0%')
        ),
        y=alt.Y("level", sort=list(nice_names["levels"].values()), title="Level"),
        color=alt.Color("covariate", legend=alt.Legend(title="Covariate"), scale=alt.Scale(scheme='pastel1'))
    )
    area = bars.mark_area(opacity=1, filled=False, fillOpacity=0, line=False, strokeWidth=2).encode(
        x=alt.value(0),
        x2=alt.value(4),
        stroke=alt.Stroke("attribute", sort=list(nice_names["attributes"].values()), legend=alt.Legend(title="Attribute")),
    )
    (
        (bars + area)
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )

def read_shares(path_to_posterior: str, interval: float, nice_names: dict[str, dict[str, str]]):
    full = az.from_netcdf(path_to_posterior)
    total_range = range_across_individuals(full, "partworths", interval)

    effect_ranges = {
        name: range_across_individuals(full, f"{name.lower()}_effect", interval)
        for name in ["Age", "Edu", "Gender", "Country", "Area", "Renewables", "Party", "Income", "Concern", "Years"]
    }

    shares = pd.DataFrame({
        name: effect_range / total_range
        for name, effect_range in effect_ranges.items()
    })
    shares = (
        shares
        .unstack()
        .reset_index()
        .rename(columns={"level_0": "covariate", 0: "share"})
    )
    shares = shares.assign(
        attribute=shares.level.str.split(":").str[0].map(nice_names["attributes"]),
        level=shares.level.str.split(":").str[1].map(nice_names["levels"])
    )
    return shares


def range_across_individuals(idata: az.InferenceData, variable_name: str, interval: float):
    individual_partworths = (
        idata
        .posterior
        .data_vars[variable_name]
        .mean(["draw", "chain"])
    )
    return (
        individual_partworths.quantile(q=1 - (1 - interval) / 2, dim="respondent")
        - individual_partworths.quantile(q=(1 - interval) / 2, dim="respondent")
    ).to_series()


if __name__ == "__main__":
    visualise_covariates(
        path_to_posterior=snakemake.input.posterior,
        interval=snakemake.params.interval,
        nice_names=snakemake.params.nice_names,
        path_to_plot=snakemake.output[0]
    )

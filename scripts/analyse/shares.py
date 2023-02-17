import altair as alt
import pandas as pd


DARK_GREY = "#424242"


def plot_likert_items(path_to_data: str, items: dict[str, str], by_country: bool,
                      colors: list[str], path_to_plot: str):
    dtype = infer_dtype(path_to_data, items)
    data = read_data(path_to_data, items, dtype=dtype, by_country=by_country)

    base = alt.Chart(
        data
    ).transform_calculate(
        signed_percentage="if(datum.type === 'Disagree' || datum.type === 'Strongly disagree', datum.percentage, 0) + if(datum.type === 'Neither agree nor disagree', datum.percentage / 2, 0)",
        Question="split(datum.item, '||')"
    ).transform_stack(
        stack="percentage",
        as_=["v1", "v2"],
        groupby=["item", "country"]
    ).transform_joinaggregate(
        offset="sum(signed_percentage)",
        groupby=["item", "country"]
    ).transform_calculate(
        nx1="datum.v1 - datum.offset",
        nx2="datum.v2 - datum.offset"
    ).encode(
        x=alt.X("nx1", type="quantitative", title="Percentage", axis=alt.Axis(format="%"), scale={"domain": [-1, 1]}),
        x2=alt.X2("nx2"),
        color=alt.Color("type", type="nominal", title="Response", scale=alt.Scale(domain=dtype.categories.values, range=colors)),
    ).properties(
        width=300
    )
    finalise_and_write_plot(base, by_country, facet_by="Question", path_to_plot=path_to_plot)


def plot_agreement_items(path_to_data: str, items: dict[str, str], by_country: bool, path_to_plot: str):
    dtype = infer_dtype(path_to_data, items)
    data = read_data(path_to_data, items, dtype, by_country=by_country)

    base = alt.Chart(
        data
    ).transform_calculate(
        Question="split(datum.item, '||')"
    ).encode(
        x=alt.X("sum(percentage)", type="quantitative", title="Percentage", axis=alt.Axis(format="%"), scale={"domain": [0, 1]}, sort=dtype.categories.values),
        color=alt.Color("type", type="ordinal", title="Response", sort=dtype.categories.values),
        order=alt.Order('type_order', sort="ascending")
    ).properties(
        width=300
    )
    finalise_and_write_plot(base, by_country, facet_by="Question", path_to_plot=path_to_plot)


def plot_demographics_item(path_to_data: str, items: dict[str, str], category_colors: list[str],
                           by_country: bool, path_to_plot: str):
    dtype = infer_dtype(path_to_data, items)
    name = list(items.values())[0]
    data = read_data(path_to_data, items, dtype=dtype, by_country=by_country)
    if dtype.ordered:
        chart_colors = alt.Color("type", type="ordinal", title=name, sort=dtype.categories.values)
    else:
        chart_colors = alt.Color("type", type="ordinal", title=name, scale=alt.Scale(domain=dtype.categories.values, range=category_colors))

    base = alt.Chart(
        data
    ).encode(
        x=alt.X("sum(percentage)", type="quantitative", title="Percentage", axis=alt.Axis(format="%"), scale={"domain": [0, 1]}, sort=dtype.categories.values),
        color=chart_colors,
        order=alt.Order('type_order', sort="ascending")
    ).properties(
        width=300
    )
    finalise_and_write_plot(base, by_country, facet_by=None, path_to_plot=path_to_plot)


def finalise_and_write_plot(chart: alt.Chart, by_country: bool, facet_by: str, path_to_plot: str):
    if by_country:
        chart = chart.encode(y=alt.Y("country", type="nominal", title="Country"))
    chart = chart.mark_bar()
    if facet_by:
        chart = chart.facet(f"{facet_by}:N", columns=2)
    (
        chart
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )


def infer_dtype(path_to_data, items):
    return pd.read_feather(path_to_data).loc[:, items.keys()].dtypes[0]


def read_data(path_to_data: str, items: dict[str, str], dtype: pd.CategoricalDtype, by_country: bool) -> pd.DataFrame:
    data = pd.read_feather(path_to_data)
    counts = [
        data.loc[:, [q, "RESPONDENT_COUNTRY"]].groupby("RESPONDENT_COUNTRY").value_counts(normalize=False)
        for q in items.keys()
    ]
    df_counts = (
        pd
        .concat(counts, axis=1)
        .rename(columns={i: q for i, q in enumerate(items.values())})
        .stack()
        .reset_index()
        .rename(columns={0: "percentage", "level_2": "item", "level_1": "type", "RESPONDENT_COUNTRY": "country"})
    )
    if len(items) == 1:
        df_counts = df_counts.rename(columns={list(items.keys())[0]: "type"})
    df_counts = (
        df_counts
        .assign(
            type=df_counts.type.astype(dtype),
            type_order=df_counts.type.astype(dtype).cat.codes
        )
        .sort_values(by=["item", "country", "type"])
    )
    if by_country:
        df_shares = (
            df_counts
            .groupby(["country", "item", "type"])
            .sum("percentage")
            .groupby(["country", "item"])
            .transform(lambda x: x / x.sum())
            .reset_index()
        )
    else:
        df_shares = (
            df_counts
            .groupby(["item", "type"])
            .sum("percentage")
            .groupby(["item"])
            .transform(lambda x: x / x.sum())
            .reset_index()
        )
    return df_shares


if __name__ == "__main__":
    match (snakemake.params.type):
        case ("likert"):
            plot_likert_items(
                path_to_data=snakemake.input.data,
                path_to_plot=snakemake.output[0],
                items=snakemake.params.plot_items,
                by_country=snakemake.params.by_country,
                colors=snakemake.params.colors
            )
        case ("agreement"):
            plot_agreement_items(
                path_to_data=snakemake.input.data,
                path_to_plot=snakemake.output[0],
                by_country=snakemake.params.by_country,
                items=snakemake.params.plot_items,
            )
        case ("demographics"):
            plot_demographics_item(
                path_to_data=snakemake.input.data,
                path_to_plot=snakemake.output[0],
                by_country=snakemake.params.by_country,
                category_colors=snakemake.params.category_colors,
                items=snakemake.params.plot_items,
            )
        case _:
            raise ValueError(f"Unknown plot type {snakemake.params.type}.")

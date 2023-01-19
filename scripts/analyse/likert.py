import altair as alt
import pandas as pd
from pandas.api.types import CategoricalDtype


DARK_GREY = "#424242"

likert = CategoricalDtype(
    categories=["Strongly disagree", "Disagree", "Neither agree nor disagree", "Agree", "Strongly agree"],
    ordered=True
)


def plot_likert_items(path_to_data: str, items: dict[str, str], colors: list[str], path_to_plot: str):
    data = read_data(path_to_data, items)

    base = alt.Chart(
        data
    ).transform_calculate(
        signed_percentage="if(datum.type === 'Disagree' || datum.type === 'Strongly disagree', datum.percentage, 0) + if(datum.type === 'Neither agree nor disagree', datum.percentage / 2, 0)",
        Question="split(datum.question, '||')"
    ).transform_stack(
        stack="percentage",
        as_=["v1", "v2"],
        groupby=["question", "country"]
    ).transform_joinaggregate(
        offset="sum(signed_percentage)",
        groupby=["question", "country"]
    ).transform_calculate(
        nx1="datum.v1 - datum.offset",
        nx2="datum.v2 - datum.offset"
    ).encode(
        x=alt.X("nx1", type="quantitative", title="Percentage", axis=alt.Axis(format="%"), scale={"domain": [-1, 1]}),
        x2=alt.X2("nx2"),
        y=alt.Y("country", type="nominal", title="Country"),
        color=alt.Color("type", type="nominal", title="Response", scale=alt.Scale(domain=likert.categories.values, range=colors)),
    ).properties(
        width=300
    )

    (
        base
        .mark_bar()
        .facet("Question:N", columns=2)
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )


def read_data(path_to_data: str, items: dict[str, str]) -> pd.DataFrame:
    data = pd.read_feather(path_to_data)
    percentages = [
        data.loc[:, [q, "RESPONDENT_COUNTRY"]].groupby("RESPONDENT_COUNTRY").value_counts(normalize=True)
        for q in items.keys()
    ]
    df_percentages = (
        pd
        .concat(percentages, axis=1)
        .rename(columns={i: q for i, q in enumerate(items.values())})
        .stack()
        .reset_index()
        .rename(columns={0: "percentage", "level_2": "question", "level_1": "type", "RESPONDENT_COUNTRY": "country"})
    )
    df_percentages = (
        df_percentages
        .assign(
            type=df_percentages.type.astype(likert),
        )
        .sort_values(by=["question", "country", "type"])
    )
    return df_percentages


if __name__ == "__main__":
    plot_likert_items(
        path_to_data=snakemake.input.data,
        path_to_plot=snakemake.output[0],
        items=snakemake.params.likert_items,
        colors=snakemake.params.colors
    )

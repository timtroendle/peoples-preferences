import altair as alt
import pandas as pd


DARK_GREY = "#424242"
LABEL_LIMIT = 300
WIDTH_SINGLE = 393


def design_plot(path_to_data: str, level_nice_names: dict[str, str], path_to_plot: str):
    df = pd.read_feather(path_to_data)
    df["level_0"] = df["level_0"].map(level_nice_names)
    df["level_1"] = df["level_1"].map(level_nice_names)

    chart = visualise_single_sector(
        df,
        level_order=list(level_nice_names.values()),
        y_orientation="left",
    )
    (
        chart
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom")
        .save(path_to_plot)
    )


def visualise_single_sector(data: pd.DataFrame, level_order: list[str], y_orientation: str):
    return (
        alt
        .Chart(data)
        .properties(
            width=WIDTH_SINGLE,
            height=WIDTH_SINGLE
        )
        .mark_rect()
        .encode(
            x=alt.X('level_0:O', sort=level_order, title="Level"),
            y=alt.Y('level_1:O', sort=level_order, title="Level",
                    axis=alt.Axis(orient=y_orientation, labelLimit=LABEL_LIMIT)),
            color=alt.Color('frequency:Q', title="Occurrence", scale=alt.Scale(scheme="greys"), legend=alt.Legend(format='%'))
        )
    )


if __name__ == "__main__":
    design_plot(
        path_to_data=snakemake.input.data,
        level_nice_names=snakemake.params.level_nice_names,
        path_to_plot=snakemake.output[0]
    )

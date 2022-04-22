import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def design_validation(path_to_data, path_to_plot):
    df = pd.read_feather(path_to_data)
    probs = pd.concat(axis="index", objs=[
        pd.crosstab(columns=df["LABEL"], index=df["TECHNOLOGY"], normalize="index").assign(Attribute="Technology"),
        pd.crosstab(columns=df["LABEL"], index=df["TRANSMISSION"], normalize="index").assign(Attribute="Transmission"),
        pd.crosstab(columns=df["LABEL"], index=df["LAND"], normalize="index").assign(Attribute="Land"),
        pd.crosstab(columns=df["LABEL"], index=df["SHARE_IMPORTS"], normalize="index").assign(Attribute="Share imports"),
        pd.crosstab(columns=df["LABEL"], index=df["PRICES"], normalize="index").assign(Attribute="Prices"),
        pd.crosstab(columns=df["LABEL"], index=df["OWNERSHIP"], normalize="index").assign(Attribute="Ownership"),
    ]).reset_index().set_index(["Attribute", "index"]).mul(100)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots()

    sns.heatmap(
        probs,
        ax=ax,
        vmin=0,
        vmax=100,
        center=50,
        cmap="vlag",
        annot=True,
        fmt=".0f",
        cbar_kws={"label": "Probability (%)"}
    )
    fig.tight_layout()
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    design_validation(
        path_to_data=snakemake.input.data,
        path_to_plot=snakemake.output[0]
    )

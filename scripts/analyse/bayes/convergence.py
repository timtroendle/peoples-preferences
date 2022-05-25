from importlib.resources import path


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def convergence_plot(path_to_betas: str, burn_in_length: int, path_to_plot: str):
    betas = pd.read_feather(path_to_betas)

    avg_respondent = betas.groupby(["iteration", "parameter"]).mean()

    g = sns.relplot(
        data=avg_respondent.reset_index(),
        x="iteration",
        y="value",
        col="parameter",
        hue="parameter",
        col_wrap=4,
        height=3,
        legend=False,
        kind="line"
    )
    for ax in g.axes:
        ax.axvline(burn_in_length, color="k", linestyle=":")

    g.fig.tight_layout()
    g.fig.savefig(path_to_plot)


if __name__ == "__main__":
    convergence_plot(
        path_to_betas=snakemake.input.betas,
        burn_in_length=snakemake.params.burn_in,
        path_to_plot=snakemake.output[0]
    )

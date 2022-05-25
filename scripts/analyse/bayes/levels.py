import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def level_part_worths(path_to_betas: str, burn_in_length: int, path_to_plot: str):
    betas = pd.read_feather(path_to_betas)
    avg_parameter = betas.set_index("iteration").iloc[burn_in_length:].groupby("parameter").mean()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    sns.barplot(
        data=avg_parameter.reset_index(),
        y="parameter",
        x="value",
        ax=ax
    )
    ax.set_xlabel("Part-worth utility")
    ax.set_ylabel("Level")

    fig.tight_layout()
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    level_part_worths(
        path_to_betas=snakemake.input.betas,
        burn_in_length=snakemake.params.burn_in,
        path_to_plot=snakemake.output[0]
    )

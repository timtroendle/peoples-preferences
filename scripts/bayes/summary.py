import arviz as az


def summary(path_to_inference_data: str, path_to_summary: str, hdi_prob: float):
    inference_data = az.from_netcdf(path_to_inference_data)["poststratify"]

    (
        az
        .summary(
            inference_data,
            hdi_prob=hdi_prob,
        )
        .rename_axis(index="parameter")
        .reset_index()
        .to_feather(path_to_summary)
    )


if __name__ == "__main__":
    summary(
        path_to_inference_data=snakemake.input.inference_data,
        path_to_summary=snakemake.output.feather,
        hdi_prob=snakemake.params.hdi_prob
    )

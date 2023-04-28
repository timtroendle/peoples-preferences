import arviz as az
import numpy as np
import xarray as xr

EFFECT_PREFIX = "effect_"


def max_subgroup_effect(inference_data: az.InferenceData, features_to_exclude: list[str]) -> az.InferenceData:
    """Retrieve largest effect across subgroups of covariate models.

    Expects effects to be named `effect_{feature}`.
    """

    effect_names = list(filter(
        lambda var: var.startswith(EFFECT_PREFIX) and feature_of_effect(var) not in features_to_exclude,
        inference_data.posterior.data_vars)
    )
    feature_names = list(map(feature_of_effect, effect_names))

    total_effect = (
        inference_data
        .posterior[effect_names]
        .to_array("effect")
        .sum("effect")
        .stack(subgroup=feature_names)
    )
    mean_total_effect = total_effect.mean(["draw", "chain"])

    subgroups_with_max_effect = np.fabs(mean_total_effect).argmax("subgroup")

    max_effect = xr.Dataset({
        level.item(): total_effect.isel(subgroup=subgroup, drop=True).sel(level=level, drop=True)
        for level, subgroup in zip(subgroups_with_max_effect.level, subgroups_with_max_effect.values)
    }).to_array("level")

    return az.convert_to_inference_data(
        max_effect.rename("max_subgroup_effect"),
        group="subgroups"
    )


def feature_of_effect(effect_name: str):
    return effect_name.split(EFFECT_PREFIX)[-1]


if __name__ == "__main__":
    (
        max_subgroup_effect(
            inference_data=az.from_netcdf(snakemake.input.inference_data),
            features_to_exclude=snakemake.params.exclude
        )
        .to_netcdf(snakemake.output[0])
    )

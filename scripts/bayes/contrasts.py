import itertools

import arviz as az
import pandas as pd
import xarray as xr


EFFECT_PREFIX = "effect_"


def all_contrasts(inference_data: az.InferenceData, features_to_exclude: list[str], hdi_prob: float) -> pd.DataFrame:
    effect_names = list(filter(
        lambda var: var.startswith(EFFECT_PREFIX) and feature_of_effect(var) not in features_to_exclude,
        inference_data.posterior.data_vars)
    )

    return pd.concat([
        contrast_for_effect(inference_data, effect_name, hdi_prob)
        for effect_name in effect_names
    ]).reset_index(drop=True)


def contrast_for_effect(inference_data: az.InferenceData, effect_name: str, hdi_prob: float) -> pd.DataFrame:
    effect = az.extract(inference_data, group="posterior", var_names=effect_name)
    mean = effect.rename({effect_dim(effect): "subgroup"})

    type_ind = pd.MultiIndex.from_tuples(
        filter(lambda x: x[1] > x[0], itertools.product(range(len(mean.subgroup)), range(len(mean.subgroup)))),
        names=["subgroup1", "subgroup2"]
    )
    type_ind1 = xr.DataArray(type_ind.to_frame().subgroup1.values, dims="contrast")
    type_ind2 = xr.DataArray(type_ind.to_frame().subgroup2.values, dims="contrast")

    contrast = (
        (mean[:, type_ind1].rename(subgroup="subgroup1")
         - mean[:, type_ind2].rename(subgroup="subgroup2"))
    ).rename("contrast")

    hdi = (
        az
        .hdi(contrast, input_core_dims=[["sample"]], hdi_prob=hdi_prob)["contrast"]
        .to_dataframe()
        .set_index(["subgroup1", "subgroup2"], append=True)
        .unstack("hdi")
        .loc[:, "contrast"]
    )
    mean = (
        contrast
        .mean("sample")
        .to_dataframe()
        .set_index(["subgroup1", "subgroup2"], append=True)["contrast"]
        .rename("expected")
    )

    return pd.DataFrame({
        "expected": mean,
        "lower": hdi["lower"],
        "higher": hdi["higher"]
    }).assign(effect=effect_name).reset_index()


def effect_dim(data: xr.DataArray) -> str:
    dim = set(data.dims) - {"level", "sample"}
    assert len(dim) == 1
    return dim.pop()


def feature_of_effect(effect_name: str):
    return effect_name.split(EFFECT_PREFIX)[-1]


if __name__ == "__main__":
    contrasts = all_contrasts(
        inference_data=az.from_netcdf(snakemake.input.inference_data),
        features_to_exclude=snakemake.params.exclude,
        hdi_prob=snakemake.params.hdi_prob
    )
    contrasts.to_feather(snakemake.output[0])

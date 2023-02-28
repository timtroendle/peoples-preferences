import arviz as az
import xarray as xr


NONE_BASE_LEVELS = { # TODO take from config
    "TECHNOLOGY": [
        "Open-field PV",
        "Wind"
    ],
    "LAND": [
        "1%",
        "2%",
        "4%",
        "8%"
    ],
    "TRANSMISSION": [
        "+0% .",
        "+25% .",
        "+50% .",
        "+75% ."
    ],
    "SHARE_IMPORTS": [
        "10%",
        "50%",
        "90%"
    ],
    "PRICES": [
        "+15%",
        "+30%",
        "+45%",
        "+60%"
    ],
    "OWNERSHIP": [
        "Community",
        "Private"
    ]
}


def AMCE(path_to_prediction_data: str, path_to_output: str):
    data = az.from_netcdf(path_to_prediction_data)
    amces = calculate_AMCE(data) # TODO handle left-hand side intercept
    amces.to_netcdf(path_to_output, group="amce")


def calculate_AMCE(data: az.InferenceData) -> xr.DataArray:
    all_amces = [
        calculate_AMCE_for_attribute(data, [f"{attribute}:{level}" for level in levels])
        for attribute, levels in NONE_BASE_LEVELS.items()
    ]
    return xr.concat(all_amces, dim="level")


def calculate_AMCE_for_attribute(data: az.InferenceData, attribute_levels: list[str]) -> xr.DataArray:
    all_amces = [
        calculate_AMCE_for_level(data, level, attribute_levels).assign_coords(level=level)
        for level in attribute_levels
    ]
    return xr.concat(all_amces, dim="level")


def calculate_AMCE_for_level(data: az.InferenceData, attribute_level: str,
                             all_none_base_levels: list[str]) -> xr.DataArray:
    # This function implements the difference-in-means estimator from [@Hainmueller:2014]
    return (
        probability_of_choice_of_level(data, attribute_level)
        - probability_of_choice_of_base_level(data, all_none_base_levels)
    )


def probability_of_choice_of_level(data: az.InferenceData, level: str) -> xr.DataArray:
    p_left = data.predictions.p_left

    level_is_left = data.predictions_constant_data.dummies_left.sel(level=level) == 1
    level_is_right = data.predictions_constant_data.dummies_right.sel(level=level) == 1

    return probability_of_choice(p_left, level_is_left, level_is_right)


def probability_of_choice_of_base_level(data: az.InferenceData, other_attribute_levels: list[str]) -> xr.DataArray:
    p_left = data.predictions.p_left

    level_is_left = (data.predictions_constant_data.dummies_left.sel(level=other_attribute_levels) == 0).all("level")
    level_is_right = (data.predictions_constant_data.dummies_right.sel(level=other_attribute_levels) == 0).all("level")

    return probability_of_choice(p_left, level_is_left, level_is_right)


def probability_of_choice(p_left: xr.DataArray, is_left: xr.DataArray, is_right: xr.DataArray):
    p_right = 1 - p_left
    return (
        xr
        .concat([p_left.where(is_left), p_right.where(is_right)], dim="choice_situation")
        .mean("choice_situation")
        .rename("p")
    )


if __name__ == "__main__":
    AMCE(
        path_to_prediction_data=snakemake.input.data,
        path_to_output=snakemake.output[0]
    )

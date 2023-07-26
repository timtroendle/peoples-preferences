import altair as alt
import arviz as az
import geopandas as gpd
import pandas as pd


MAP_MIN_X = 2550000
MAP_MIN_Y = 1550000
MAP_MAX_X = 5400000
MAP_MAX_Y = 3950000
EXTENT_ROI_FEATURE = {
    "type": "Feature",
    "geometry": {"type": "Polygon",
                 "coordinates": [[
                     [MAP_MAX_X, MAP_MAX_Y],
                     [MAP_MAX_X, MAP_MIN_Y],
                     [MAP_MIN_X, MAP_MIN_Y],
                     [MAP_MIN_X, MAP_MAX_Y],
                     [MAP_MAX_X, MAP_MAX_Y]]]},
    "properties": {}
}
WORLD_COLOR = "#E5E5E5"
WIDTH = 470


def map_regional_heterogeinity(inference_data: az.InferenceData, regions: gpd.GeoDataFrame, world: gpd.GeoDataFrame,
                               level: str) -> alt.Chart:
    regions = regions.to_crs(epsg=3035)
    world = world.to_crs(epsg=3035).cx[MAP_MIN_X:MAP_MAX_X, MAP_MIN_Y:MAP_MAX_Y]
    mean_partworths = (
        inference_data["poststratify"]
        .mean(["chain", "draw"])
        .regional_partworths
        .sel(level=level)
        .to_series()
        .reset_index()
        .merge(regions[["shapeID", "shapeGroup"]], left_on="region", right_on="shapeID")
    )
    mean_partworths = mean_partworths[mean_partworths.country == mean_partworths.shapeGroup]
    return world_chart(world) + regions_chart(regions, mean_partworths)


def world_chart(world: gpd.GeoDataFrame) -> alt.Chart:
    return (
        alt
        .Chart(world.cx[MAP_MIN_X:MAP_MAX_X, MAP_MIN_Y:MAP_MAX_Y], width=WIDTH)
        .mark_geoshape(filled=True, fill=WORLD_COLOR, stroke="white", strokeWidth=1, clip=True)
        .project(type="identity", reflectY=True, fit=EXTENT_ROI_FEATURE)
    )


def regions_chart(regions: gpd.GeoDataFrame, mean_partworths: pd.DataFrame) -> alt.Chart:
    return (
        alt
        .Chart(regions.to_crs(epsg=3035), width=WIDTH)
        .mark_geoshape(clip=True)
        .project(type="identity", reflectY=True, fit=EXTENT_ROI_FEATURE)
        .transform_lookup(
            lookup='shapeID',
            from_=alt.LookupData(data=mean_partworths, key='region', fields=['regional_partworths'])
        ).encode(
            alt.Color('regional_partworths:Q', title="Partworth", sort="descending")
        )
    )


if __name__ == "__main__":
    chart = map_regional_heterogeinity(
        inference_data=az.from_netcdf(snakemake.input.data),
        regions=gpd.read_feather(snakemake.input.regions),
        world=gpd.read_file(snakemake.input.world),
        level=snakemake.wildcards.level.replace("___", ":").replace("__", " ")
    )
    chart.save(snakemake.output[0])

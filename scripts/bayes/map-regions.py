import altair as alt
import arviz as az
import geopandas as gpd


def map_regional_heterogeinity(inference_data: az.InferenceData, regions: gpd.GeoDataFrame, level: str) -> alt.Chart:
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

    return (
        alt
        .Chart(regions.to_crs(epsg=3035), width=600, height=500)
        .mark_geoshape()
        .project(type="identity", reflectY=True)
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
        level=snakemake.wildcards.level.replace("___", ":").replace("__", " ")
    )
    chart.save(snakemake.output[0])

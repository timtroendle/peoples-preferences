import arviz as az
import altair as alt


DARK_GREY = "#424242"


def plot_density(path_to_inference_data: str, variable_name: str, nice_variable_name: str, path_to_plot: str):
    data = az.from_netcdf(path_to_inference_data)

    df = data.posterior[variable_name].to_dataframe().reset_index()

    base = (
        alt
        .Chart(
            df
        ).transform_density(
            variable_name,
            groupby=['chain'],
            as_=[variable_name, 'density']
        ).encode(
            x=alt.X(f"{variable_name}:Q", title=nice_variable_name),
            color=alt.Color("chain:N", title="Chain"),
            y=alt.Y("density:Q", title="Density")
        ).properties(
            width=300,
            height=200
        ).mark_area(
            opacity=0,
            line=True
        )
        .configure(font="Lato")
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .save(path_to_plot)
    )




if __name__ == "__main__":
    plot_density(
        path_to_inference_data=snakemake.input.data,
        path_to_plot=snakemake.output[0],
        variable_name=snakemake.params.variable_name,
        nice_variable_name=snakemake.params.nice_variable_name
    )

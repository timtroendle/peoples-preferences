require(cluster)
library(ggplot2)
library(dplyr)
require(arrow)
library(uwot)


cluster <- function(path_data, all_features, feature, path_plot) {
    conjoint <- read_feather(path_data)
    respondents <- conjoint  %>%
        group_by(RESPONDENT_ID) %>%
        summarise_all(first)
    dist <- daisy(respondents[, all_features], metric = "gower")

    plot_umap(dist, respondents, feature, path_plot)
}

plot_umap <- function(dist, respondents, feature, path_plot) {
    set.seed(0815)
    d_umap <- as.data.frame(umap(dist, n_neighbors = 10))

    p <- ggplot(data = d_umap, aes(x = V1, y = V2, col = respondents[[feature]])) +
        geom_point() +
        ggtitle("UMAP") +
        scale_color_discrete(feature)
    ggsave(path_plot, p)
}

cluster(
    snakemake@input[["data"]],
    snakemake@params[["features"]],
    snakemake@wildcards[["feature"]],
    snakemake@output[[1]]
)

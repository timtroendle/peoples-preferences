require(cluster)
library(dplyr)
library(ggplot2)
require(arrow)
library(uwot)


cluster <- function(path_data, features, n_cluster, path_cluster, path_plot_tree, path_plot_silh, path_plot_umap) {
    conjoint <- read_feather(path_data)
    respondents <- conjoint  %>%
        group_by(RESPONDENT_ID) %>%
        summarise_all(first)
    dist <- daisy(respondents[, features], metric = "gower")

    tree <- calc_and_plot_tree(dist, path_plot_tree)
    cluster_tree <- calc_and_plot_cluster(tree, n_cluster, dist, path_plot_silh)
    plot_umap(dist, cluster_tree, path_plot_umap)

    clustered <- data.frame(
        RESPONDENT_ID = respondents$RESPONDENT_ID,
        cluster = factor(cluster_tree)
    )
    all_data <- conjoint %>% left_join(clustered, by = "RESPONDENT_ID")
    write_feather(all_data, path_cluster)
}

calc_and_plot_tree <- function(dist, path_plot) {
    tree <- hclust(dist, method = "ward.D2")
    png(path_plot)
    plot(tree)
    dev.off()
    tree
}

calc_and_plot_cluster <- function(tree, n_cluster, dist, path_plot) {
    cluster_tree <- cutree(tree, k = n_cluster)
    png(path_plot)
    plot(silhouette(cluster_tree, dist), col = "black")
    dev.off()
    cluster_tree
}

plot_umap <- function(dist, cluster_tree, path_plot) {
    set.seed(0815)
    d_umap <- as.data.frame(umap(dist, n_neighbors = 10))
    d_umap$cluster <- factor(cluster_tree)

    p <- ggplot(data = d_umap, aes(x = V1, y = V2, col = cluster)) +
        geom_point() +
        ggtitle("UMAP")
    ggsave(path_plot, p)
}

cluster(
    snakemake@input[["data"]],
    snakemake@params[["features"]],
    snakemake@params[["n_cluster"]],
    snakemake@output[["data"]],
    snakemake@output[["tree"]],
    snakemake@output[["silh"]],
    snakemake@output[["umap"]]
)

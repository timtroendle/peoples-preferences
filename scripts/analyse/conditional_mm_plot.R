library("cregg")
library("ggplot2")

conditional_mm_plot <- function(path.data, path.plot, estimate, by, factors, codes, cuts) {
    d.conjoint <- read.csv(path.data)

    for (factor_name in factors) {
        d.conjoint[[factor_name]] <- factor(d.conjoint[[factor_name]])
    }
    levels(d.conjoint$TRANSMISSION) <- paste(levels(d.conjoint$TRANSMISSION), ".")
    for (attribute in names(codes)) {
        d.conjoint[[attribute]] <- factor(
            d.conjoint[[attribute]],
            labels=codes[[attribute]]
        )
    }
    for (attribute in names(cuts)) {
        d.conjoint[[attribute]] <- cut(
            d.conjoint[[attribute]],
            breaks = cuts[[attribute]]$breaks,
            labels = cuts[[attribute]]$labels
        )
    }

    d.conjoint <- preprocess_conditional_attribut(d.conjoint, by)

    f1 <- CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + PRICES + TRANSMISSION + OWNERSHIP

    by.sym <- as.symbol(by)
    mm_diff.conjoint <- eval(bquote(cj( # eval and bquote necessary because of dynamic var "by"
        data = d.conjoint,
        formula = f1,
        estimate = estimate,
        id = ~RESPONDENT_ID,
        weights = ~WEIGHT,
        by = ~.(by.sym)
    )))
    if (estimate == "mm") {
        p <- plot(mm_diff.conjoint, group = by, vline = 0.5)
    }
    else {
        p <- plot(mm_diff.conjoint) +  xlab(paste(c("Estimated Difference of ", mm_diff.conjoint$BY[1]), collapse=" "))
    }
    ggsave(path.plot, p)
}


preprocess_conditional_attribut <- function(data, by) {
    data <- droplevels(data[data[[by]] != "other", ])
    data <- droplevels(data[data[[by]] != "do not know", ])
    data <- droplevels(data[data[[by]] != "no answer", ])
}


conditional_mm_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["estimate"]],
    snakemake@params[["by"]],
    snakemake@params[["factors"]],
    snakemake@params[["codes"]],
    snakemake@params[["cuts"]]
)

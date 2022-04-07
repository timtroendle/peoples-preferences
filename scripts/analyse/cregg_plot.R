library("cregg")
library("ggplot2")

cregg_plot <- function(path.data, path.plot, estimate, factors) {
    d.conjoint <- read.csv(path.data)

    for (factor_name in factors) {
        d.conjoint[[factor_name]] <- factor(d.conjoint[[factor_name]])
    }
    levels(d.conjoint$TRANSMISSION) <- paste(levels(d.conjoint$TRANSMISSION), ".")

    f1 <- CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + PRICES + TRANSMISSION + OWNERSHIP

    cj.conjoint <- cj(
        data=d.conjoint,
        formula=f1,
        id = ~RESPONDENT_ID,
        estimate = estimate,
        weights = ~WEIGHT,
    )
    if (estimate == "mm") {
        p <- plot(cj.conjoint, vline = 0.5, xlim=c(0.2,.8))
    }
    else {
        p <- plot(cj.conjoint)
    }
    ggsave(path.plot, p)
}

cregg_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["estimate"]],
    snakemake@params[["factors"]]
)

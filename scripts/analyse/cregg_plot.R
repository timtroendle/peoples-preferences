library("cregg")
library("ggplot2")

cregg_plot <- function(path.data, path.plot, estimate) {
    d.conjoint <- read.csv(path.data)
    d.conjoint$TECHNOLOGY <- factor(
        d.conjoint$TECHNOLOGY,
    )
    d.conjoint$SHARE_IMPORTS <- factor(
        d.conjoint$SHARE_IMPORTS,
    )
    d.conjoint$LAND <- factor(
        d.conjoint$LAND,
    )
    d.conjoint$PRICES <- factor(
        d.conjoint$PRICES,
    )
    d.conjoint$TRANSMISSION <- factor(
        d.conjoint$TRANSMISSION,
        levels=c("-25%", "+0%", "+25%", "+50%", "+75%")
    )
    levels(d.conjoint$TRANSMISSION) <- c("-25.0%", "+0.0%", "+25.0%", "+50.0%", "+75.0%")
    d.conjoint$OWNERSHIP <- factor(
        d.conjoint$OWNERSHIP,
    )

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
    snakemake@params[["estimate"]]
)

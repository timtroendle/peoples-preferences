library(arrow)
library(cjoint)

amce_plot <- function(path.data, path.plot) {
    d.conjoint <- read_feather(path.data)

    results <- amce(CHOICE_INDICATOR ~  TECHNOLOGY + SHARE_IMPORTS + LAND
        + PRICES + TRANSMISSION + OWNERSHIP,
        data=d.conjoint,
        cluster=TRUE,
        respondent.id="RESPONDENT_ID",
        weight="WEIGHT") # TODO verify this is correct
    png(path.plot)
    plot(
        results,
        xlab="Change in Pr(Design preferred)",
        xlim=c(-.5,.5),
        breaks=c(-.4, 0, .4),
        labels=c("-.4","0",".4"),
        text.size=13
    )
    dev.off()

}

amce_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]]
)

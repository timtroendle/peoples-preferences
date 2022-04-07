library("cjoint")

amce_plot <- function(path.data, path.plot, factors) {
    d.conjoint <- read.csv(path.data)

    for (factor_name in factors) {
        d.conjoint[[factor_name]] <- factor(d.conjoint[[factor_name]])
    }
    levels(d.conjoint$TRANSMISSION) <- paste(levels(d.conjoint$TRANSMISSION), ".")
    levels(d.conjoint$TRANSMISSION)[1] <- "-25.0% ."
    d.conjoint$TECHNOLOGY <- relevel(d.conjoint$TECHNOLOGY, "Rooftop PV")
    d.conjoint$TRANSMISSION <- relevel(d.conjoint$TRANSMISSION, "+0% .")
    d.conjoint$OWNERSHIP <- relevel(d.conjoint$OWNERSHIP, "Public")

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
    snakemake@output[[1]],
    snakemake@params[["factors"]]
)

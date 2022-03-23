library("cjoint")

acie_plot <- function(path.data, path.plot, formula) {
    d.conjoint <- read.csv(path.data)
    d.conjoint$TECHNOLOGY <- factor(
        d.conjoint$TECHNOLOGY,
        levels=c("Rooftop PV", "Open-field PV", "Wind")
    )
    d.conjoint$SHARE_IMPORTS <- factor(
        d.conjoint$SHARE_IMPORTS,
        levels=c("0%", "10%", "50%", "90%")
    )
    d.conjoint$LAND <- factor(
        d.conjoint$LAND,
        levels=c("0.5%", "1%", "2%", "4%", "8%")
    )
    d.conjoint$PRICES <- factor(
        d.conjoint$PRICES,
        levels=c("+0%", "+15%", "+30%", "+45%", "+60%")
    )
    d.conjoint$TRANSMISSION <- factor(
        d.conjoint$TRANSMISSION,
        levels=c("+0%", "-25%", "+25%", "+50%", "+75%")
    )
    levels(d.conjoint$TRANSMISSION)[2] <- "-25.0%"
    d.conjoint$OWNERSHIP <- factor(
        d.conjoint$OWNERSHIP,
        levels=c("Public", "Community", "Private")
    )
    results <- amce(
        as.formula(formula),
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

acie_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["formula"]]
)

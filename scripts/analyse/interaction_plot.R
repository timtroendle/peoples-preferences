library(arrow)
library(cjoint)

acie_plot <- function(path.data, path.plot, formula) {
    d.conjoint <- read_feather(path.data)

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
        xlim=c(-.35,.35),
        breaks=c(-.2, 0, .2),
        labels=c("-.2","0",".2"),
        text.size=13
    )
    dev.off()

}

acie_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["formula"]]
)

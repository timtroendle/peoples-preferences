library(arrow)
library(cregg)
library(ggplot2)

conditional_mm_plot <- function(path.data, path.plot, estimate, by) {
    d.conjoint <- read_feather(path.data)

    f1 <- CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + PRICES + TRANSMISSION + OWNERSHIP

    by.sym <- as.symbol(by)
    mm_diff.conjoint <- eval(bquote(cj( # eval and bquote necessary because of dynamic var "by"
        data = d.conjoint,
        formula = f1,
        estimate = estimate,
        id = ~RESPONDENT_ID,
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

conditional_mm_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["estimate"]],
    snakemake@params[["by"]]
)

library("cregg")
library("ggplot2")

H6_plot <- function(path.data, path.plot) {
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
    d.conjoint$Q6_AREA <- factor(
        d.conjoint$Q6_AREA,
        labels=c("urban", "rural", "no_answer")
    )
    d.conjoint.urban_rural <- droplevels(d.conjoint[d.conjoint$Q6_AREA != "no_answer", ])

    f1 <- CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + PRICES + TRANSMISSION + OWNERSHIP

    mm_diff.conjoint.urban_rural <- mm_diffs(
        data=d.conjoint.urban_rural,
        formula=f1,
        id = ~RESPONDENT_ID,
        weights = ~WEIGHT,
        by = ~Q6_AREA
    )
    p <- plot(mm_diff.conjoint.urban_rural)
    ggsave(path.plot, p)
}

H6_plot(snakemake@input[["data"]], snakemake@output[[1]])

library("cregg")
library("ggplot2")

conditional_mm_plot <- function(path.data, path.plot, estimate, by, codes) {
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
    d.conjoint$RESPONDENT_COUNTRY <- factor(d.conjoint$RESPONDENT_COUNTRY)
    for (attribute in names(codes)) {
        d.conjoint[[attribute]] <- factor(
            d.conjoint[[attribute]],
            labels=codes[[attribute]]
        )
    }
    if (by == "Q3_GENDER") {
        d.conjoint <- preprocess_Q3_GENDER(d.conjoint)
    }
    else if (by == "Q6_AREA") {
        d.conjoint <- preprocess_Q6_AREA(d.conjoint)
    }

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


preprocess_Q3_GENDER <- function(data) {
    droplevels(data[data$Q3_GENDER != "other", ])
}


preprocess_Q6_AREA <- function(data) {
    droplevels(data[data$Q6_AREA != "no_answer", ])
}


conditional_mm_plot(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["estimate"]],
    snakemake@params[["by"]],
    snakemake@params[["codes"]]
)

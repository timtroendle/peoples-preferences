library(arrow)
library(dplyr)
library(missForest)


impute <- function(path_to_data, seed, features_to_impute, path_to_output_data, path_to_output_error) {
    data <- read_feather(path_to_data)
    respondents <- data %>%
        group_by(RESPONDENT_ID) %>%
        summarise_all(first) %>%
        select(starts_with(c("RESPONDENT_ID", "RESPONDENT_COUNTRY", "Q"))) %>%
        select(-c(
            Q2_AGE_TRESHOLD, # factor with single level not supported
            Q5_POSTCODE, # string not supported
            Q10_INCOME, # leads to failure, unsure why
            Q12_PARTY, # leads to failure, unsure why
            Q26_FEEDBACK # string not supported
        ))

    set.seed(seed)
    respondents_imputed <- missForest(as.data.frame(respondents %>% select(-RESPONDENT_ID)), variablewise = TRUE)

    # write data
    imputed <- respondents_imputed$ximp[features_to_impute] %>%
        rename_with(~ paste0(.x, "_imputed")) %>%
        mutate(RESPONDENT_ID = respondents$RESPONDENT_ID)
    data <- data %>% left_join(imputed, by = "RESPONDENT_ID")
    write_feather(data, path_to_output_data)

    # write error
    error <- data.frame(
        feature = colnames(respondents_imputed$ximp),
        OOBerror = as.numeric(as.list(respondents_imputed$OOBerror))
    )
    write_feather(error, path_to_output_error)
}


impute(
    snakemake@input[["data"]],
    snakemake@params[["seed"]],
    snakemake@params[["features"]],
    snakemake@output[["data"]],
    snakemake@output[["error"]]
)

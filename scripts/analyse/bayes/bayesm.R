library(arrow)
library(bayesm)
library(dplyr)
library(tidyr)

sink(snakemake@log[[1]], type = c("output", "message"))


bayesm <- function(path_to_data, path_to_betas, formula, n_iterations, keep) {
    conjoint <- read_feather(path_to_data)
    conjoint$RESPONDENT_ID <- factor(conjoint$RESPONDENT_ID)

    dat <- preprocess_data(conjoint, formula)

    data  <- list(lgtdata = dat, p = 2)
    prior <- list(ncomp = 1)
    mcmc  <- list(R = n_iterations, nprint = 0, keep = keep)

    out <- rhierMnlRwMixture(Data = data, Prior = prior, Mcmc = mcmc)

    betas <- excavate_betas(
        out = out,
        respondents = levels(factor(conjoint$RESPONDENT_ID)),
        parameter_names = colnames(dat[[1]]$X),
        keep = keep
    )

    write_feather(betas, path_to_betas)

}


preprocess_data <- function(conjoint, formula) {
    dummy <- model.matrix(as.formula(formula), data = conjoint)
    dummy <- dummy[, colnames(dummy)[2:length(colnames(dummy))]] # remove intercept

    N <- nlevels(conjoint$RESPONDENT_ID)
    N <- 1000 # FIXME remove
    dat <- vector(mode = "list", length = N)
    for (i in 1:N) {
        respondent_mask <- conjoint$RESPONDENT_ID == levels(conjoint$RESPONDENT_ID)[i]
        dat[[i]]$y <- (conjoint[respondent_mask, ] %>%
            group_by(CHOICE_SET) %>%
            summarise(choice = last(CHOICE_INDICATOR) + 1))$choice
        dat[[i]]$X <- dummy[respondent_mask, ]
    }
    return(dat)
}


excavate_betas <- function(out, respondents, parameter_names, keep) {
    betas <- NULL
    for (i in 1:dim(out$betadraw)[2]) {
        beta <- data.frame(out$betadraw[1:dim(out$betadraw)[1], i, 1:dim(out$betadraw)[3]])
        colnames(beta) <- lapply(colnames(beta), integer_iteration)
        beta$respondent <- respondents[1:dim(out$betadraw)[1]]
        beta <- beta %>% pivot_longer(!respondent, names_to = "iteration", values_to = "value")
        beta$iteration <- as.integer(beta$iteration) * keep
        beta$parameter <- parameter_names[i]
        if (is.null(betas)) {
            betas <- beta
        }
        else {
            betas <- rbind(betas, beta)
        }
    }
    betas$respondent <- as.factor(betas$respondent)
    betas$parameter <- as.factor(betas$parameter)
    return(betas)
}


integer_iteration <- function(iteration_name) {
    return(as.integer(substr(iteration_name, 2, nchar(iteration_name))))
}


bayesm(
    snakemake@input[["data"]],
    snakemake@output[[1]],
    snakemake@params[["formula"]],
    snakemake@params[["n_iterations"]],
    snakemake@params[["keep"]]
)

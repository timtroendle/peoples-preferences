rule logistic_regression:
    message: "Fit a simplistic logistic regression model."
    input:
        data = rules.global_conjoint.output[0]
    params:
        n_tune = 2000,
        n_draws = 2000,
        random_seed = 4000,
    resources:
        runtime = 60
    threads: 4
    output: "build/results/models/logistic-regression/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/logistic.py"


rule multinomial_logit:
    message: "Sample {wildcards.dist} from a multinomial logit model."
    input: data = rules.global_conjoint.output[0]
    params:
        n_tune = config["models"]["multinomial"]["n-tune"],
        n_draws = config["models"]["multinomial"]["n-draws"],
        limit_respondents = config["models"]["multinomial"]["limit-respondents"],
        random_seed = config["models"]["multinomial"]["random-seed"],
    resources:
        runtime = 60,
        memory = 4000
    threads: 4
    output: "build/results/models/multinomial-logit/{dist}/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/multinomial.py"


def hierarchical_model_config(param_name):
    def hierarchical_model_config(wildcards):
        name = f"hierarchical-{wildcards.name}"
        return config["models"][name][param_name]
    return hierarchical_model_config


rule hierarchical:
    message: "Sample {wildcards.dist} from hierarchical Bayes model '{wildcards.name}' using PyMC."
    input: data = rules.global_conjoint.output[0]
    params:
        n_tune = hierarchical_model_config("n-tune"),
        n_draws = hierarchical_model_config("n-draws"),
        limit_respondents = hierarchical_model_config("limit-respondents"),
        random_seed = hierarchical_model_config("random-seed"),
        individual_covariates = hierarchical_model_config("individual-covariates"),
        covariances = hierarchical_model_config("covariances")
    resources:
        runtime = hierarchical_model_config("runtime"),
        memory = lambda wildcards, threads: hierarchical_model_config("memory")(wildcards) // threads
    threads: 4
    output: "build/results/models/hierarchical-{name}/{dist}/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/hierarchical.py"


rule diagnostics:
    message: "Run diagnostics for hierarchical model {wildcards.name}."
    input:
        inference_data = rules.hierarchical.output[0]
    params: hdi_prob = 0.94
    output:
        trace = "build/results/models/hierarchical-{name}/{dist}/diagnostics/trace.png",
        pop_means = "build/results/models/hierarchical-{name}/{dist}/diagnostics/pop-means.png",
        forest = "build/results/models/hierarchical-{name}/{dist}/diagnostics/forest.png",
        summary = "build/results/models/hierarchical-{name}/{dist}/diagnostics/summary.feather",
        rhos_individual = "build/results/models/hierarchical-{name}/{dist}/diagnostics/rhos-individual.png",
        rhos_country = "build/results/models/hierarchical-{name}/{dist}/diagnostics/rhos-country.png",
        individuals = "build/results/models/hierarchical-{name}/{dist}/diagnostics/individuals.png",
        confusion = "build/results/models/hierarchical-{name}/{dist}/diagnostics/confusion-matrix.csv",
        accuracy = "build/results/models/hierarchical-{name}/{dist}/diagnostics/in-sample-prediction-accuracy.txt",
        probability = "build/results/models/hierarchical-{name}/{dist}/diagnostics/choice-probability.png"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/diagnostics.py"


rule visualise_partworths:
    message: "Visualise population mean partworths of multinomial logit model."
    input: data = rules.multinomial_logit.output[0]
    params:
        facet_by_country = True,
        variable_names = "partworths",
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/multinomial-logit/{dist}/pop-means.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_population_means:
    message: "Visualise population mean partworths of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "alpha",
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/hierarchical-{name}/{dist}/pop-means.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_country_differences:
    message: "Visualise country differences of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "countries",
        facet_by_country = True,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/hierarchical-{name}/{dist}/country-differences.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_country_means:
    message: "Visualise country means of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = ["alpha", "countries"],
        facet_by_country = True,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/hierarchical-{name}/{dist}/country-means.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_partworths_heterogeneity:
    message: "Visualise heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "partworths",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/hierarchical-{name}/{dist}/individual-partworths.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_unexplained_heterogeneity:
    message: "Visualise unexplained heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "individuals",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/hierarchical-{name}/{dist}/unexplained-heterogeneity.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_covariates:
    message: "Visualise the explaining power of covariates."
    input: data = "build/models/hierarchical-covariates/{dist}/inference-data.nc"
    params:
        interval = 0.9, # show only this share of the total interval
        nice_names = config["report"]["nice-names"],
    output: "build/results/models/hierarchical-covariates/{dist}/covariates.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/covariates.py"


rule visualise_mu_left_intercept:
    message: "Density plot of left intercept."
    input:
        data = rules.hierarchical.output[0]
    params:
        variable_name = "mu_left_intercept",
        nice_variable_name = "Mean additional partworth utility of left option"
    output: "build/results/models/hierarchical-{name}/{dist}/left-option.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/density.py"


rule visualise_varying_left_option_effect:
    message: "Visualise heterogeneity of left-hand option partworth of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        varying_variable_name = "left_intercept",
        pop_means_variable_name = "mu_left_intercept",
        level_choice = {},
        narrow_hdi = config["report"]["hdi-prob"]["narrow"],
        wide_hdi = config["report"]["hdi-prob"]["default"]
    output: "build/results/models/hierarchical-{name}/{dist}/varying/left-intercept.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/varying.py"


rule visualise_varying_attribute_levels:
    message: "Visualise heterogeneity of attribute level partworth of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        varying_variable_name = "partworths",
        pop_means_variable_name = "alpha",
        level_choice = lambda wildcards, output: {"level": wildcards["level"].replace("__", " ")},
        narrow_hdi = config["report"]["hdi-prob"]["narrow"],
        wide_hdi = config["report"]["hdi-prob"]["default"]
    output: "build/results/models/hierarchical-{name}/{dist}/varying/{level}.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/varying.py"

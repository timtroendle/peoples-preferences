def memory_requirements_bayes(wildcards):
    n_iterations = config["models"][wildcards.model]["n-iterations"]
    keep = config["models"][wildcards.model]["keep"]
    n_respondents = 4000
    n_parameters = 20
    n_data_points = n_iterations * n_respondents * n_parameters / keep
    mega_byte = n_data_points * 64 / 8 / 1_000_000
    requested_memory = mega_byte * 5
    return requested_memory if requested_memory > 16000 else 16000


def memory_requirements_plots(wildcards):
    memory_bayes = memory_requirements_bayes(wildcards)
    return memory_bayes if memory_bayes > 16000 else 16000


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
    output: "build/models/logistic-regression/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/logistic.py"


rule multinomial_logit:
    message: "Fit a multinomial logit model."
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
    output: "build/models/multinomial-logit/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/multinomial.py"


def hierarchical_model_config(param_name):
    def hierarchical_model_config(wildcards):
        name = f"hierarchical-{wildcards.name}"
        return config["models"][name][param_name]
    return hierarchical_model_config


rule hierarchical:
    message: "Fit hierarchical Bayes model '{wildcards.name}' using PyMC."
    input: data = rules.global_conjoint.output[0]
    params:
        n_tune = hierarchical_model_config("n-tune"),
        n_draws = hierarchical_model_config("n-draws"),
        limit_respondents = hierarchical_model_config("limit-respondents"),
        random_seed = hierarchical_model_config("random-seed"),
        individual_covariates = hierarchical_model_config("individual-covariates")
    resources:
        runtime = hierarchical_model_config("runtime"),
        memory = lambda wildcards, threads: hierarchical_model_config("memory")(wildcards) // threads
    threads: 4
    output: "build/models/hierarchical-{name}/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/hierarchical.py"


rule diagnostics:
    message: "Run diagnostics for hierarchical model {wildcards.name}."
    input:
        inference_data = rules.hierarchical.output[0]
    params: hdi_prob = 0.94
    output:
        trace = "build/models/hierarchical-{name}/diagnostics/trace.png",
        pop_means = "build/models/hierarchical-{name}/diagnostics/pop-means.png",
        forest = "build/models/hierarchical-{name}/diagnostics/forest.png",
        summary = "build/models/hierarchical-{name}/diagnostics/summary.feather",
        rhos_individual = "build/models/hierarchical-{name}/diagnostics/rhos-individual.png",
        rhos_country = "build/models/hierarchical-{name}/diagnostics/rhos-country.png",
        individuals = "build/models/hierarchical-{name}/diagnostics/individuals.png",
        confusion = "build/models/hierarchical-{name}/diagnostics/confusion-matrix.csv",
        accuracy = "build/models/hierarchical-{name}/diagnostics/in-sample-prediction-accuracy.txt"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/diagnostics.py"


rule visualise_partworths:
    message: "Visualise population mean partworths of multinomial logit model."
    input: posterior = rules.multinomial_logit.output[0]
    params:
        facet_by_country = True,
        variable_name = "partworths",
        hdi_prob = config["report"]["hdi_prob"],
        nice_names = config["report"]["nice-names"],
    output: "build/models/multinomial-logit/pop-means.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_population_means:
    message: "Visualise population mean partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        variable_name = "alpha",
        hdi_prob = config["report"]["hdi_prob"],
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/pop-means.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_country_differences:
    message: "Visualise heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        variable_name = "countries",
        facet_by_country = True,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi_prob"],
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/country-differences.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_partworths_heterogeneity:
    message: "Visualise heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        variable_name = "partworths",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/individual-partworths.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_unexplained_heterogeneity:
    message: "Visualise unexplained heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        variable_name = "individuals",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/unexplained-heterogeneity.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_covariates:
    message: "Visualise the explaining power of covariates."
    input: posterior = "build/models/hierarchical-covariates/inference-data.nc"
    params:
        interval = 0.9, # show only this share of the total interval
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-covariates/covariates.vega.json"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/covariates.py"

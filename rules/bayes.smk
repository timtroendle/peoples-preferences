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


rule bayesm_model:
    message: "Fit hierarchical Bayes model {wildcards.model}."
    input:
        data = rules.global_conjoint.output[0]
    params:
        formula = lambda wildcards: config["models"][f"{wildcards.model}"]["formula"],
        n_iterations = lambda wildcards: config["models"][f"{wildcards.model}"]["n-iterations"],
        keep = lambda wildcards: config["models"][f"{wildcards.model}"]["keep"],
    log: "build/logs/{model}/bayesm.log"
    output: "build/models/{model}/betas.feather"
    resources:
        runtime = 240,
        memory = memory_requirements_bayes
    conda: "../envs/bayesm.yaml"
    script: "../scripts/analyse/bayes/bayesm.R"


rule convergence_plot:
    message: "Plot convergence of {wildcards.model} parameters."
    input:
        betas = rules.bayesm_model.output[0]
    params: burn_in = lambda wildcards: config["models"][f"{wildcards.model}"]["burn-in-length"],
    output: "build/models/{model}/convergence.png"
    resources:
        memory = memory_requirements_plots
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/bayes/convergence.py"


rule level_part_worth_utility:
    message: "Plot part-worth utility of all {wildcards.model} levels."
    input:
        betas = rules.bayesm_model.output[0]
    params: burn_in = lambda wildcards: config["models"][f"{wildcards.model}"]["burn-in-length"],
    output: "build/models/{model}/level-part-worths.png"
    resources:
        memory = memory_requirements_plots
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/bayes/levels.py"


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
    script: "../scripts/analyse/bayes/logistic.py"


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
    script: "../scripts/analyse/bayes/multinomial.py"


def hierarchical_model_config(param_name):
    def hierarchical_model_config(wildcards):
        name = f"hierarchical-{wildcards.name}"
        return config["models"][name][param_name]
    return hierarchical_model_config


def memory_requirements_hierarchical(wildcards, threads):
    n_draws = hierarchical_model_config("n-draws")(wildcards)
    if hierarchical_model_config("limit-respondents")(wildcards):
        n_respondents = hierarchical_model_config("limit-respondents")(wildcards)
    else:
        n_respondents = 4000
    req_per_respondent_and_draw = 0.0064 # (MB), empirically derived
    if hierarchical_model_config("individual-covariates")(wildcards):
        req_per_respondent_and_draw = req_per_respondent_and_draw * 5
    total_memory = req_per_respondent_and_draw * n_respondents * n_draws
    requested_memory = (total_memory / threads) * 1.2
    return requested_memory if requested_memory > 16000 else 16000


def runtime_hierarchical(wildcards):
    n_tune = hierarchical_model_config("n-tune")(wildcards)
    n_draws = hierarchical_model_config("n-draws")(wildcards)
    n_iterations = n_tune + n_draws
    if hierarchical_model_config("limit-respondents")(wildcards):
        n_respondents = hierarchical_model_config("limit-respondents")(wildcards)
    else:
        n_respondents = 4000
    req_per_respondent_and_iteration = 1 / 30000 # (min), empirically derived
    if hierarchical_model_config("individual-covariates")(wildcards):
        req_per_respondent_and_iteration = req_per_respondent_and_iteration * 5
    return req_per_respondent_and_iteration * n_respondents * n_iterations * 1.2


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
        runtime = runtime_hierarchical,
        memory = memory_requirements_hierarchical
    threads: 4
    output: "build/models/hierarchical-{name}/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/analyse/bayes/hierarchical.py"


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
        rhos = "build/models/hierarchical-{name}/diagnostics/rhos.png",
        individuals = "build/models/hierarchical-{name}/diagnostics/individuals.png"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/diagnostics.py"


rule visualise_partworths:
    message: "Visualise population mean partworths of multinomial logit model."
    input: posterior = rules.multinomial_logit.output[0]
    params:
        title="{hdi_prob:.0f}% highest density interval of population-level partworths".format(
            hdi_prob=config["report"]["hdi_prob"] * 100),
        facet_by_country = True,
        variable_name = "partworths",
        hdi_prob = config["report"]["hdi_prob"],
        nice_names = config["report"]["nice-names"],
    output: "build/models/multinomial-logit/pop-means.{figure_format}"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/partworths.py"


rule visualise_population_means:
    message: "Visualise population mean partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        title="{hdi_prob:.0f}% highest density interval of population-level partworths".format(
            hdi_prob=config["report"]["hdi_prob"] * 100),
        variable_name = "alpha",
        hdi_prob = config["report"]["hdi_prob"],
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/pop-means.{figure_format}"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/partworths.py"


rule visualise_partworths_heterogeneity:
    message: "Visualise heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        title="Range of average individual-level partworths",
        variable_name = "partworths",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/individual-partworths.{figure_format}"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/partworths.py"


rule visualise_unexplained_heterogeneity:
    message: "Visualise uenxplained heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: posterior = rules.hierarchical.output[0]
    params:
        title="Unexplained heterogeneity in individual-level partworths",
        variable_name = "individuals",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-{name}/unexplained-heterogeneity.{figure_format}"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/partworths.py"


rule visualise_covariates:
    message: "Visualise the explaining power of covariates."
    input: posterior = "build/models/hierarchical-covariates/inference-data.nc"
    params:
        interval = 0.9, # show only this share of the total interval
        nice_names = config["report"]["nice-names"],
    output: "build/models/hierarchical-covariates/covariates.{figure_format}"
    resources:
        runtime = 60,
        memory = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/covariates.py"

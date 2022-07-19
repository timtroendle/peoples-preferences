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
    output: "build/models/logistic-regression.nc"
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
        runtime = 30,
        memory = 4000
    threads: 4
    output: "build/models/multinomial-logit.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/analyse/bayes/multinomial.py"


rule hierarchical:
    message: "Fit a hierarchical Bayes model using PyMC."
    input: data = rules.global_conjoint.output[0]
    params:
        n_tune = config["models"]["hierarchical"]["n-tune"],
        n_draws = config["models"]["hierarchical"]["n-draws"],
        limit_respondents = config["models"]["hierarchical"]["limit-respondents"],
        random_seed = config["models"]["hierarchical"]["random-seed"]
    resources:
        runtime = 1440
    threads: 4
    output: "build/models/hierarchical.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/analyse/bayes/hierarchical.py"


rule diagnostics:
    message: "Run diagnostics for hierarchical model."
    input:
        inference_data = rules.hierarchical.output[0]
    output:
        trace = "build/models/hierarchical/trace.png",
        pop_means = "build/models/hierarchical/pop-means.png",
        forest = "build/models/hierarchical/forest.png",
        summary = "build/models/hierarchical/summary.csv",
        rhos = "build/models/hierarchical/rhos.png"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/diagnostics.py"


rule visualise_partworths:
    message: "Visualise partworths."
    input: posterior = rules.multinomial_logit.output[0]
    output: "build/models/multinomial-logit.{figure_format}"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/bayes/partworths.py"

rule bayesm_model:
    message: "Fit hierarchical Bayes model {wildcards.model}."
    input:
        data = rules.global_conjoint.output[0]
    params:
        formula = lambda wildcards: config["models"][f"{wildcards.model}"]["formula"],
        n_iterations = lambda wildcards: config["models"][f"{wildcards.model}"]["n-iterations"]
    log: "build/logs/{model}/bayesm.log"
    output: "build/{model}/betas.feather"
    conda: "../envs/bayesm.yaml"
    script: "../scripts/analyse/bayes/bayesm.R"


rule convergence_plot:
    message: "Plot convergence of {wildcards.model} parameters."
    input:
        betas = rules.bayesm_model.output[0]
    params: burn_in = lambda wildcards: config["models"][f"{wildcards.model}"]["burn-in-length"],
    output: "build/{model}/convergence.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/bayes/convergence.py"


rule level_part_worth_utility:
    message: "Plot part-worth utility of all {wildcards.model} levels."
    input:
        betas = rules.bayesm_model.output[0]
    params: burn_in = lambda wildcards: config["models"][f"{wildcards.model}"]["burn-in-length"],
    output: "build/{model}/level-part-worths.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/bayes/levels.py"

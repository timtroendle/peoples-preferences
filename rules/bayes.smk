rule synthetic_data:
    message:
        "Create synthetic input data."
    input:
        original = rules.global_conjoint.output[0]
    params:
        n_respondents = config["parameters"]["synthetic-data"]["n-respondents"],
        n_sets_per_respondent = config["parameters"]["synthetic-data"]["n-choice-sets-per-respondent"],
        seed = config["parameters"]["synthetic-data"]["seed"]
    output:
        "build/data/synthetic.feather"
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/bayes/synthesise.py"


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
    message: "Sample {wildcards.sample} from a multinomial logit model."
    input: data = rules.global_conjoint.output[0]
    params:
        n_tune = config["models"]["multinomial"]["n-tune"],
        n_draws = config["models"]["multinomial"]["n-draws"],
        limit_respondents = config["models"]["multinomial"]["limit-respondents"],
        random_seed = config["models"]["multinomial"]["random-seed"],
    resources:
        runtime = 60,
        mem_mb_per_cpu = 4000
    threads: 4
    output: "build/results/models/multinomial-logit/{sample}/inference-data.nc"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/multinomial.py"


def hierarchical_model_config(param_name):
    def hierarchical_model_config(wildcards):
        name = f"hierarchical-{wildcards.name}"
        param = config["models"][name][param_name]
        try:
            return param[wildcards.sample]
        except   TypeError:
            return param
    return hierarchical_model_config


rule hierarchical:
    message: "Sample {wildcards.sample} from hierarchical Bayes model '{wildcards.name}' using PyMC."
    input:
        data = rules.global_conjoint_imputed.output[0]
    params:
        n_tune = hierarchical_model_config("n-tune"),
        n_draws = hierarchical_model_config("n-draws"),
        limit_respondents = hierarchical_model_config("limit-respondents"),
        random_seed = hierarchical_model_config("random-seed"),
        model_variety = hierarchical_model_config("model-variety"),
        covariances = hierarchical_model_config("covariances")
    resources:
        runtime = hierarchical_model_config("runtime"),
        mem_mb_per_cpu = lambda wildcards, threads: hierarchical_model_config("mem_mb")(wildcards) // threads
    threads: lambda wildcards: hierarchical_model_config("threads")(wildcards)
    output: "build/results/models/hierarchical-{name}/{sample}/inference-data.nc"
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/pymc.yaml"
    script: "../scripts/bayes/hierarchical.py"


rule predict:
    message:
        "Predict from hierarchical Bayes model '{wildcards.name}' using PyMC."
    input:
        in_sample = rules.global_conjoint_imputed.output[0],
        out_sample = rules.synthetic_data.output[0],
        trace = "build/results/models/hierarchical-{name}/posterior/inference-data.nc"
    params:
        limit_respondents = hierarchical_model_config("limit-respondents"),
        random_seed = hierarchical_model_config("random-seed"),
        model_variety = hierarchical_model_config("model-variety"),
        covariances = hierarchical_model_config("covariances"),
    resources:
        runtime = hierarchical_model_config("runtime"),
        mem_mb_per_cpu = lambda wildcards, threads: hierarchical_model_config("mem_mb")(wildcards) // threads
    threads: 1
    output:
        "build/results/models/hierarchical-{name}/{sample}/inference-data.nc"
    wildcard_constraints:
        sample = "prediction"
    conda:
        '../envs/pymc.yaml'
    script:
        '../scripts/bayes/hierarchical.py'


rule prediction_amces:
    message:
        "Calculate AMCEs from predicted probabilities."
    input:
        data = "build/results/models/hierarchical-{name}/prediction/inference-data.nc"
    output:
        "build/results/models/hierarchical-{name}/amce/inference-data.nc"
    conda:
        "../envs/pymc.yaml"
    script:
        "../scripts/bayes/amce.py"


rule poststratify:
    message:
        "Poststratify results from hierarchical Bayes model '{wildcards.name}'."
    input:
        conjoint = rules.global_conjoint_imputed.output[0],
        inference_data = "build/results/models/hierarchical-{name}/posterior/inference-data.nc",
        census = rules.census.output[0]
    params:
        limit_respondents = hierarchical_model_config("limit-respondents"),
        model_variety = hierarchical_model_config("model-variety"),
        covariances = hierarchical_model_config("covariances"),
    threads: 1
    output:
        "build/results/models/hierarchical-{name}/{sample}/inference-data.nc"
    wildcard_constraints:
        sample = "poststratify"
    conda:
        '../envs/pymc.yaml'
    script:
        '../scripts/bayes/hierarchical.py'


rule max_subgroup_effects:
    message:
        "Calculate the largest expected effects across subgroups."
    input:
        inference_data = "build/results/models/hierarchical-{name}/posterior/inference-data.nc"
    params:
        exclude = ["country", "respondent"]
    output:
        "build/results/models/hierarchical-{name}/{sample}/inference-data.nc"
    wildcard_constraints:
        sample = "subgroups"
    conda:
        "../envs/pymc.yaml"
    script:
        "../scripts/bayes/subgroups.py"


rule contrasts:
    message:
        "Calculate contrasts between all elements of all clusters."
    input:
        inference_data = "build/results/models/hierarchical-{name}/posterior/inference-data.nc"
    params:
        exclude = ["respondent"],
        hdi_prob = config["report"]["hdi-prob"]["default"]
    output:
        "build/results/models/hierarchical-{name}/{sample}/contrasts.feather"
    wildcard_constraints:
        sample = "prior|posterior"
    conda:
        "../envs/pymc.yaml"
    script:
        "../scripts/bayes/contrasts.py"


rule diagnostics:
    message: "Run diagnostics for hierarchical model {wildcards.name}."
    input:
        inference_data = rules.hierarchical.output[0]
    params:
        hdi_prob = 0.94,
        rho_plot_dir = lambda wildcards: f"build/results/models/hierarchical-{wildcards.name}/{wildcards.sample}/diagnostics/rho"
    output:
        trace = "build/results/models/hierarchical-{name}/{sample}/diagnostics/trace.png",
        pop_means = "build/results/models/hierarchical-{name}/{sample}/diagnostics/pop-means.png",
        forest = "build/results/models/hierarchical-{name}/{sample}/diagnostics/forest.png",
        summary = "build/results/models/hierarchical-{name}/{sample}/diagnostics/summary.feather",
        individuals = "build/results/models/hierarchical-{name}/{sample}/diagnostics/individuals.png",
        confusion = "build/results/models/hierarchical-{name}/{sample}/diagnostics/confusion-matrix.csv",
        accuracy = "build/results/models/hierarchical-{name}/{sample}/diagnostics/in-sample-prediction-accuracy.txt",
        probability = "build/results/models/hierarchical-{name}/{sample}/diagnostics/choice-probability.png",
        utility = "build/results/models/hierarchical-{name}/{sample}/diagnostics/utility.png",
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/diagnostics.py"


rule poststratification_summary:
    message: "Create a summary of poststratification."
    input:
        inference_data = rules.hierarchical.output[0]
    params:
        hdi_prob = 0.94,
    output:
        feather = "build/results/models/hierarchical-{name}/{sample}/diagnostics/summary.feather",
    wildcard_constraints:
        sample = "poststratify"
    conda:
        "../envs/analyse.yaml"
    script:
        "../scripts/bayes/summary.py"


rule brief_summary:
    message: "Create a brief sample summary for publication."
    input:
        feather = "build/results/models/hierarchical-{name}/{sample}/diagnostics/summary.feather"
    params:
        columns = config["report"]["brief-sample-summary"]["columns"],
        float_format_per_column = config["report"]["brief-sample-summary"]["float-format"]
    output:
        csv = "build/results/models/hierarchical-{name}/{sample}/diagnostics/summary-brief.csv",
    conda:
        "../envs/default.yaml"
    script:
        "../scripts/utils/feather_to_csv.py"


rule visualise_partworths:
    message: "Visualise population mean partworths of multinomial logit model."
    input: data = rules.multinomial_logit.output[0]
    params:
        facet_by_country = True,
        variable_names = "partworths",
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/multinomial-logit/{sample}/pop-means.vega.json")
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_intercept:
    message: "Visualise intercept of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "alpha",
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/intercept.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_country_differences:
    message: "Visualise country differences of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "effect_country",
        facet_by_country = True,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/country-differences.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_country_differences_sigma_respondent:
    message: "Visualise country differences of sigma of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "sigma_respondent",
        facet_by_country = True,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/varying-variation.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_country_means:
    message: "Visualise country means of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = ["alpha", "effect_country"],
        facet_by_country = True,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/country-means.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_population_means:
    message: "Visualise poststratified population means of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = ["pop_means"],
        facet_by_country = False,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/pop-means.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "poststratify"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_sample_partworths:
    message: "Visualise partworths of sample on the national level of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = ["sample_partworths"],
        facet_by_country = False,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/sample.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "poststratify"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_sample_bias:
    message: "Visualise bias on the national level of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = ["national_bias"],
        facet_by_country = False,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/bias.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "poststratify"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_max_subgroup_effects:
    message: "Visualise largest expected effect across subgroups of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = ["max_subgroup_effect"],
        facet_by_country = False,
        aggregate_individuals = False,
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/max-subgroup-effect.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "subgroups"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_respondent_heterogeneity:
    message: "Visualise heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "partworths",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/individual-partworths.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_regional_heterogeneity:
    message: "Visualise heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "regional_partworths",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/regional-partworths.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "poststratify"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule map_regional_heterogeneity:
    message: "Map heterogeneity of partworths of {wildcards.name} on {wildcards.layer} for {wildcards.level}."
    input:
        data = rules.hierarchical.output[0],
        regions = "build/data/geoboundaries/{layer}.feather",
        world = "data/automatic/geoboundaries/WLD/ADM0.geojson"
    output: temp("build/results/models/hierarchical-{name}/{sample}/maps/{layer}-{level}.vega.json")
    resources:
        runtime = 10,
        mem_mb_per_cpu = 16000
    wildcard_constraints:
        sample = "poststratify"
    conda: "../envs/geo.yaml"
    script: "../scripts/bayes/map-regions.py"


rule visualise_unexplained_heterogeneity:
    message: "Visualise unexplained heterogeneity of partworths of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        variable_names = "effect_respondent",
        aggregate_individuals = True,
        hdi_prob = None, # has no use here
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/unexplained-heterogeneity.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"


rule visualise_covariates:
    message: "Visualise the explaining power of covariates."
    input: data = "build/models/hierarchical-covariates/{sample}/inference-data.nc"
    params:
        interval = 0.9, # show only this share of the total interval
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-covariates/{sample}/covariates.vega.json")
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/covariates.py"


rule visualise_mu_left_intercept:
    message: "Density plot of left intercept."
    input:
        data = rules.hierarchical.output[0]
    params:
        variable_name = "mu_left_intercept",
        nice_variable_name = "Mean additional partworth utility of left option"
    output: temp("build/results/models/hierarchical-{name}/{sample}/left-option.vega.json")
    wildcard_constraints:
        sample = "prior|posterior"
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
    output: temp("build/results/models/hierarchical-{name}/{sample}/varying/left-intercept.vega.json")
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/varying.py"


rule visualise_varying_attribute_levels:
    message: "Visualise heterogeneity of attribute level partworth of hierarchical model {wildcards.name}."
    input: data = rules.hierarchical.output[0]
    params:
        varying_variable_name = "partworths",
        pop_means_variable_name = "alpha",
        level_choice = lambda wildcards, output: {"level": wildcards["level"].replace("___", ":").replace("__", " ")},
        narrow_hdi = config["report"]["hdi-prob"]["narrow"],
        wide_hdi = config["report"]["hdi-prob"]["default"]
    output: temp("build/results/models/hierarchical-{name}/{sample}/varying/{level}.vega.json")
    wildcard_constraints:
        sample = "prior|posterior"
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/varying.py"


rule visualise_amces:
    message: "Visualise AMCEs of hierarchical model {wildcards.name}."
    input: data = rules.prediction_amces.output[0]
    params:
        variable_names = "p",
        hdi_prob = config["report"]["hdi-prob"]["default"],
        nice_names = config["report"]["nice-names"],
    output: temp("build/results/models/hierarchical-{name}/{sample}/pop-means.vega.json")
    wildcard_constraints:
        sample = "amce"
    resources:
        runtime = 60,
        mem_mb_per_cpu = 64000
    conda: "../envs/analyse.yaml"
    script: "../scripts/bayes/partworths.py"

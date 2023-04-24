rule national_respondents:
    message: "Statistical overview over respondents in {wildcards.country_id}."
    input:
        data = rules.global_conjoint.output[0]
    output: "build/results/analysis/respondent-stats-{country_id}.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/respondents.py"


rule clustering:
    message: "Cluster respondents by '{wildcards.cluster}' using hierarchical clustering."
    input:
        data = rules.global_conjoint.output[0]
    params:
        features = lambda wildcards: config["clustering"][wildcards.cluster]["features"],
        n_cluster = lambda wildcards: config["clustering"][wildcards.cluster]["n-cluster"]
    output:
        data = "build/results/clustering/{cluster}/clustered.feather",
        tree = "build/results/clustering/{cluster}/tree.png",
        umap = "build/results/clustering/{cluster}/umap.png",
        silh = "build/results/clustering/{cluster}/silhouette.png"
    conda: "../envs/clustering.yaml"
    script: "../scripts/analyse/clustering.R"


rule cluster_analysis:
    message: "Understand importance of '{wildcards.feature}' for cluster '{wildcards.cluster}'."
    input:
        data = rules.clustering.output.data
    params:
        features = lambda wildcards: config["clustering"][wildcards.cluster]["features"]
    output: "build/results/clustering/{cluster}/umap-{feature}.png"
    conda: "../envs/clustering.yaml"
    script: "../scripts/analyse/umap.R"


rule conditional_mm_of_cluster:
    message: "Create conditional MM plot for cluster '{wildcards.cluster}'."
    input:
        data = rules.clustering.output.data
    params:
        estimate = "mm",
        by = "cluster"
    output: "build/results/clustering/{cluster}/conditional-mm.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule robustness_check_choice_set_number:
    message: "Create conditional MM plot for CHOICE_SET."
    input:
        data = rules.global_conjoint.output
    params:
        estimate = "mm",
        by = "CHOICE_SET"
    output: "build/results/robustness/conditional-mm-choice-set.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule robustness_check_label_number:
    message: "Create conditional MM plot for LABEL."
    input:
        data = rules.global_conjoint.output
    params:
        estimate = "mm",
        by = "LABEL"
    output: "build/results/robustness/conditional-mm-label.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule design_validation_plot:
    message: "Plot conditional probabilities of all levels."
    input:
        data = rules.global_conjoint.output[0]
    output: "build/results/robustness/design-validation.png"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/design_validation.py"


rule likert_plot:
    message: "Plot likert items."
    input:
        data = rules.global_conjoint.output[0]
    params:
        plot_items = config["report"]["nice-names"]["likert-items"],
        colors = config["report"]["colors"]["likert"],
        by_country = False,
        type = "likert"
    output: "build/results/analysis/likert-items.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/shares.py"


use rule likert_plot as likert_plot_by_country with:
    params:
        plot_items = config["report"]["nice-names"]["likert-items"],
        colors = config["report"]["colors"]["likert"],
        by_country = True,
        type = "likert"
    output: "build/results/analysis/likert-items-by-country.vega.json"


rule agreement_plot:
    message: "Plot agreement items."
    input:
        data = rules.global_conjoint.output[0]
    params:
        plot_items = config["report"]["nice-names"]["agreement-items"],
        by_country = False,
        type = "agreement"
    output: "build/results/analysis/agreement-items.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/shares.py"


use rule agreement_plot as agreement_plot_by_country with:
    params:
        plot_items = config["report"]["nice-names"]["agreement-items"],
        by_country = True,
        type = "agreement"
    output: "build/results/analysis/agreement-items-by-country.vega.json"


rule gender_plot:
    message: "Plot gender distribution."
    input:
        data = rules.global_conjoint.output[0]
    params:
        plot_items = {"Q3_GENDER": "Gender"},
        by_country = True,
        category_colors = config["report"]["colors"]["categories"],
        type = "demographics"
    output: "build/results/analysis/gender.vega.json"
    conda: "../envs/analyse.yaml"
    script: "../scripts/analyse/shares.py"


use rule gender_plot as area_plot with:
    params:
        plot_items = {"Q6_AREA": "Area"},
        by_country = True,
        category_colors = config["report"]["colors"]["categories"],
        type = "demographics"
    output: "build/results/analysis/area.vega.json"


use rule gender_plot as education_plot with:
    params:
        plot_items = {"Q9_EDUCATION": "Education"},
        by_country = True,
        category_colors = config["report"]["colors"]["categories"],
        type = "demographics"
    output: "build/results/analysis/education.vega.json"


use rule gender_plot as income_plot with:
    params:
        plot_items = {"Q10_INCOME": "Income"},
        by_country = True,
        category_colors = config["report"]["colors"]["categories"],
        type = "demographics"
    output: "build/results/analysis/income.vega.json"

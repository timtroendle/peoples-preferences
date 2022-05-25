rule amce_plot:
    message: "Create plot of AMCEs."
    input:
        script = "scripts/analyse/amce_plot.R",
        data = rules.global_conjoint.output[0]
    output: "build/amce.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/amce_plot.R"


rule mm_plot:
    message: "Create plot of MMs."
    input:
        script = "scripts/analyse/cregg_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
    output: "build/mm.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/cregg_plot.R"


rule national_respondents:
    message: "Statistical overview over respondents in {wildcards.country_id}."
    input:
        script = "scripts/analyse/respondents.py",
        data = rules.global_conjoint.output[0]
    output: "build/{country_id}/respondent-stats.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/respondents.py"


rule H1:
    message: "Create plot for H1."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS * PRICES + LAND + TRANSMISSION + OWNERSHIP",
    output: "build/H1.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H2:
    message: "Create plot for H2."
    input:
        script = "scripts/analyse/cregg_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "amce",
    output: "build/H2.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/cregg_plot.R"


rule H4:
    message: "Create plot for H4."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND * PRICES + TRANSMISSION + OWNERSHIP",
    output: "build/H4.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H5:
    message: "Create plot for H5."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + PRICES * OWNERSHIP + TRANSMISSION",
    output: "build/H5.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H6:
    message: "Create plot for H6."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm_diff",
        by = "Q6_AREA"
    output: "build/H6.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H7:
    message: "Create plot for H7."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS * LAND + PRICES + TRANSMISSION + OWNERSHIP",
    output: "build/H7.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H8:
    message: "Create plot for H8."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND*TECHNOLOGY + PRICES + TRANSMISSION + OWNERSHIP",
    output: "build/H8.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H9:
    message: "Create plot for H9."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS*TECHNOLOGY + LAND + PRICES + TRANSMISSION + OWNERSHIP",
    output: "build/H9.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H11:
    message: "Create plot for H11."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "RESPONDENT_COUNTRY",
    output: "build/H11.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H16:
    message: "Create plot for H16."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q3_GENDER",
    output: "build/H16.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H17:
    message: "Create plot for H17."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q4_BIRTH_YEAR_aggregated",
    output: "build/H17.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H18:
    message: "Create plot for H18."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q8_YEARS_REGION_aggregated",
    output: "build/H18.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H19:
    message: "Create plot for H19."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q10_INCOME_aggregated",
    output: "build/H19.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H20:
    message: "Create plot for H20."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q9_EDUCATION_aggregated",
    output: "build/H20.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H21:
    message: "Create plot for H21."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q11_CLIMATE_CONCERN_aggregated",
    output: "build/H21.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H22:
    message: "Create plot for H22."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q12_PARTY_aggregated",
    output: "build/H22.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H23:
    message: "Create plot for H23."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm",
        by = "Q7_RENEWABLES",
    output: "build/H23.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule H24:
    message: "Create plot for H24."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND * OWNERSHIP + PRICES + TRANSMISSION",
    output: "build/H24.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H25:
    message: "Create plot for H25."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS * TRANSMISSION + LAND + PRICES + OWNERSHIP",
    output: "build/H25.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H26:
    message: "Create plot for H26."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND * TRANSMISSION + PRICES + OWNERSHIP",
    output: "build/H26.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"



rule H27:
    message: "Create plot for H27."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + OWNERSHIP + PRICES * TRANSMISSION",
    output: "build/H27.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H28:
    message: "Create plot for H28."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY * TRANSMISSION + SHARE_IMPORTS + LAND + OWNERSHIP + PRICES",
    output: "build/H28.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H29:
    message: "Create plot for H29."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY * OWNERSHIP + SHARE_IMPORTS + LAND + PRICES + TRANSMISSION",
    output: "build/H29.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule cluster:
    message: "Cluster respondents by '{wildcards.cluster}' using hierarchical clustering."
    input:
        script = "scripts/analyse/cluster.R",
        data = rules.global_conjoint.output[0]
    params:
        features = lambda wildcards: config["cluster"][wildcards.cluster]["features"],
        n_cluster = lambda wildcards: config["cluster"][wildcards.cluster]["n-cluster"]
    output:
        data = "build/cluster/{cluster}/clustered.feather",
        tree = "build/cluster/{cluster}/tree.png",
        umap = "build/cluster/{cluster}/umap.png",
        silh = "build/cluster/{cluster}/silhouette.png"
    conda: "../envs/cluster.yaml"
    script: "../scripts/analyse/cluster.R"


rule cluster_analysis:
    message: "Understand importance of '{wildcards.feature}' for cluster '{wildcards.cluster}'."
    input:
        script = "scripts/analyse/umap.R",
        data = rules.cluster.output.data
    params:
        features = lambda wildcards: config["cluster"][wildcards.cluster]["features"]
    output: "build/cluster/{cluster}/umap-{feature}.png"
    conda: "../envs/cluster.yaml"
    script: "../scripts/analyse/umap.R"


rule conditional_mm_of_cluster:
    message: "Create conditional MM plot for cluster '{wildcards.cluster}'."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.cluster.output.data
    params:
        estimate = "mm",
        by = "cluster"
    output: "build/cluster/{cluster}/conditional-mm.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule robustness_check_choice_set_number:
    message: "Create conditional MM plot for CHOICE_SET."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output
    params:
        estimate = "mm",
        by = "CHOICE_SET"
    output: "build/robustness/conditional-mm-choice-set.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule robustness_check_label_number:
    message: "Create conditional MM plot for LABEL."
    input:
        script = "scripts/analyse/conditional_mm_plot.R",
        data = rules.global_conjoint.output
    params:
        estimate = "mm",
        by = "LABEL"
    output: "build/robustness/conditional-mm-label.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"


rule design_validation_plot:
    message: "Plot conditional probabilities of all levels."
    input:
        script = "scripts/analyse/design_validation.py",
        data = rules.global_conjoint.output[0]
    output: "build/robustness/design-validation.png"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/design_validation.py"

rule amce_plot:
    message: "Create plot of AMCEs in {wildcards.country_id}."
    input:
        script = "scripts/analyse/amce_plot.R",
        data = rules.national_conjoint.output[0]
    output: "build/{country_id}/amce.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/amce_plot.R"


rule global_amce_plot:
    message: "Create plot of AMCEs."
    input:
        script = "scripts/analyse/amce_plot.R",
        data = rules.global_conjoint.output[0]
    output: "build/amce.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/amce_plot.R"


rule mm_plot:
    message: "Create plot of MMs in {wildcards.country_id}."
    input:
        script = "scripts/analyse/mm_plot.R",
        data = rules.national_conjoint.output[0]
    output: "build/{country_id}/mm.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/mm_plot.R"


rule global_mm_plot:
    message: "Create plot of MMs."
    input:
        script = "scripts/analyse/cregg_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "mm"
    output: "build/mm.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/cregg_plot.R"


rule respondents:
    message: "Statistical overview over respondents in {wildcards.country_id}."
    input:
        script = "scripts/analyse/respondents.py",
        data = rules.national_conjoint.output[0]
    params:
        codes = {k:v for k,v in config["codes"].items() if k in ["Q3_GENDER", "Q6_AREA", "Q10_INCOME"]}
    output: "build/{country_id}/respondent-stats.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/analyse/respondents.py"


rule H1:
    message: "Create plot for H1."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS * PRICES + LAND + TRANSMISSION + OWNERSHIP"
    output: "build/H1.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H2:
    message: "Create plot for H2."
    input:
        script = "scripts/analyse/cregg_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        estimate = "amce"
    output: "build/H2.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/cregg_plot.R"


rule H4:
    message: "Create plot for H4."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND * PRICES + TRANSMISSION + OWNERSHIP"
    output: "build/H4.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H5:
    message: "Create plot for H5."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND + PRICES * OWNERSHIP + TRANSMISSION"
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
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS * LAND + PRICES + TRANSMISSION + OWNERSHIP"
    output: "build/H7.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H8:
    message: "Create plot for H8."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS + LAND*TECHNOLOGY + PRICES + TRANSMISSION + OWNERSHIP"
    output: "build/H8.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/interaction_plot.R"


rule H9:
    message: "Create plot for H9."
    input:
        script = "scripts/analyse/interaction_plot.R",
        data = rules.global_conjoint.output[0]
    params:
        formula = "CHOICE_INDICATOR ~ TECHNOLOGY + SHARE_IMPORTS*TECHNOLOGY + LAND + PRICES + TRANSMISSION + OWNERSHIP"
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
        by = "RESPONDENT_COUNTRY"
    output: "build/H11.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/analyse/conditional_mm_plot.R"

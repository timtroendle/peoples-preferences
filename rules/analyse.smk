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

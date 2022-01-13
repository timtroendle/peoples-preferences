rule amce_plot:
    message: "Create plot of AMCEs in {wildcards.country}."
    input:
        script = "scripts/amce_plot.R",
        data = "data/raw-data-{country}.csv"
    output: "build/{country}/amce.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/amce_plot.R"


rule respondents:
    message: "Statistical overview over respondents in {wildcards.country}."
    input:
        script = "scripts/respondents.py",
        data = "data/raw-data-{country}.csv"
    output: "build/{country}/respondent-stats.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/respondents.py"

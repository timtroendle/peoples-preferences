rule amce_plot:
    message: "Create plot of AMCEs in {wildcards.country}."
    input:
        script = "scripts/amce_plot.R",
        data = "data/raw-data-{country}.csv"
    output: "build/{country}/amce-plot.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/amce_plot.R"

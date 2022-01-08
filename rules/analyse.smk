rule amce_plot:
    message: "Create plot of AMCEs."
    input:
        script = "scripts/amce_plot.R",
        data = "data/rds_prod.experiment.184514.stacked.csv"
    output: "build/amce-plot.png"
    conda: "../envs/cjoint.yaml"
    script: "../scripts/amce_plot.R"

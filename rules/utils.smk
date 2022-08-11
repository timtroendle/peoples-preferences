rule feather_to_csv:
    message: "Transform {wildcards.filename}.feather to csv."
    input:
        feather = "build/{filename}.feather"
    output: "build/{filename}.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/utils/feather_to_csv.py"

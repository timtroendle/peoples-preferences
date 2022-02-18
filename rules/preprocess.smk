rule national_conjoint:
    message: "Preprocess conjoint data for country {wildcards.country_id}."
    input:
        script = "scripts/preprocess/national.py",
        conjointly = config["data-sources"]["conjointly"],
        respondi = config["data-sources"]["respondi"]
    params:
        population = lambda wildcards: config["parameters"]["population-count"][wildcards.country_id],
        pre_test_threshold = config["parameters"]["pre-test-threshold"]
    output: "build/{country_id}/conjoint.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/national.py"


rule global_conjoint:
    message: "Merge all national conjoint datasets."
    input:
        script = "scripts/preprocess/merge.py",
        datasets = expand("build/{country_id}/conjoint.csv", country_id=COUNTRY_IDS)
    output: "build/conjoint.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/merge.py"

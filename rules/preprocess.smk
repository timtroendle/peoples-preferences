rule national_conjoint_raw:
    message: "Preprocess conjoint data for country {wildcards.country_id}."
    input:
        conjointly = config["data-sources"]["conjointly"],
        respondi = config["data-sources"]["respondi"]
    params:
        population = lambda wildcards: config["parameters"]["population-count"][wildcards.country_id],
        pre_test_threshold = config["parameters"]["pre-test-threshold"],
        q12_party_base = config["parameters"]["Q12-party-base"]
    output: "build/{country_id}/raw.feather"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/national.py"


rule global_conjoint_raw:
    message: "Merge all national conjoint datasets."
    input:
        datasets = expand("build/{country_id}/raw.feather", country_id=COUNTRY_IDS)
    output: "build/raw.feather"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/merge.py"


rule global_conjoint:
    message: "Adjust data types and derive new features."
    input:
        data = rules.global_conjoint_raw.output[0]
    params:
        types = config["data-types"],
        aggregated_levels = config["new-features"]["aggregated-levels"],
        categorised_levels = config["new-features"]["categorised-levels"]
    output: "build/conjoint.feather"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/types.py"

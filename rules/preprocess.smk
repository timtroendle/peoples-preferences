import pycountry


def alpha3to2(alpha3: str) -> str:
    return pycountry.countries.lookup(alpha3).alpha_2


rule download_geonames:
    message: "Download geonames file for country {wildcards.country_id}."
    params:
        url = lambda wildcards: config["data-sources"]["geonames"].format(country_id=alpha3to2(wildcards["country_id"]))
    output:
        zip = "data/automatic/geonames/{country_id}.zip"
    shell:
        "curl -sLo {output.zip} '{params.url}'"


rule geonames:
    message: "Unzip geonames file for country {wildcards.country_id}."
    input:
        zip = rules.download_geonames.output[0]
    params:
        alpha2 = lambda wildcards: alpha3to2(wildcards.country_id)
    output:
        txt = "build/data/geonames/{country_id}.txt"
    shell: """
        unzip -o {input.zip} -d build/data/geonames
        mv build/data/geonames/{params.alpha2}.txt {output.txt}
        """


rule download_german_census:
    message: "Download German census data {wildcards.table}."
    params:
        url = lambda wildcards: config["data-sources"]["zensus"]["url"].format(
            user=os.environ["ZENSUS_USER"],
            password=os.environ["ZENSUS_PASSWORD"],
            table=wildcards["table"]
        )
    output:
        data = "data/automatic/zensus/{table}.csv"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/zensus.py"


rule german_census:
    message: "Preprocess German census data for feature {wildcards.feature}."
    input:
        data = lambda wildcards: "data/automatic/zensus/{table}.csv".format(
            table=config["data-sources"]["zensus"]["features"][wildcards.feature]
        )
    params:
        name_mapping = lambda wildcards: config["census"]["DEU"][wildcards.feature]
    output:
        "build/data/census/DEU/{feature}.feather"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/zensus.py"


rule census:
    message: "Merge all census files."
    input:
        features = expand("build/data/census/DEU/{feature}.feather", feature=config["data-sources"]["zensus"]["features"])
    output:
        "build/data/census.nc"
    conda:
        "../envs/preprocess.yaml"
    script:
        "../scripts/preprocess/census.py"


rule national_conjoint_raw:
    message: "Preprocess conjoint data for country {wildcards.country_id}."
    input:
        conjointly = config["data-sources"]["conjointly"],
        respondi = config["data-sources"]["respondi"],
        geonames = rules.geonames.output[0]
    params:
        pre_test_threshold = config["parameters"]["pre-test-threshold"],
        q12_party_base = config["parameters"]["Q12-party-base"]
    output: "build/data/{country_id}.feather"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/national.py"


rule global_conjoint_raw:
    message: "Merge all national conjoint datasets."
    input:
        datasets = expand("build/data/{country_id}.feather", country_id=COUNTRY_IDS)
    output: "build/data/raw.feather"
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
    output: "build/data/conjoint.feather"
    conda: "../envs/preprocess.yaml"
    script: "../scripts/preprocess/types.py"


rule global_conjoint_imputed:
    message: "Impute missing data."
    input:
        data = rules.global_conjoint.output[0]
    params:
        seed = config["parameters"]["impute"]["seed"],
        features = config["parameters"]["impute"]["features"]
    output:
        data = "build/data/conjoint-imputed.feather",
        error = "build/data/imputation-oob-error.feather"
    conda:
        "../envs/impute.yaml"
    script:
        "../scripts/preprocess/impute.R"

import datetime
from pathlib import Path

from snakemake.utils import min_version

PANDOC = "pandoc --filter pantable --filter pandoc-crossref --citeproc"
COUNTRY_IDS = ["DEU", "POL", "PRT", "DNK"]
NONE_BASELINE_ATTRIBUTE_LEVELS = [
    'TECHNOLOGY___Open-field__PV', 'TECHNOLOGY___Wind',
    'LAND___1%', 'LAND___2%', 'LAND___4%', 'LAND___8%',
    'TRANSMISSION___+0%__.', 'TRANSMISSION___+25%__.', 'TRANSMISSION___+50%__.', 'TRANSMISSION___+75%__.',
    'SHARE_IMPORTS___10%', 'SHARE_IMPORTS___50%', 'SHARE_IMPORTS___90%',
    'PRICES___+15%', 'PRICES___+30%', 'PRICES___+45%', 'PRICES___+60%',
    'OWNERSHIP___Community', 'OWNERSHIP___Private'
]

configfile: "config/default.yaml"
include: "rules/preprocess.smk"
include: "rules/analyse.smk"
include: "rules/bayes.smk"
include: "rules/utils.smk"
include: "./rules/sync.smk"
localrules: all, report, clean
wildcard_constraints:
    country_id = "|".join(COUNTRY_IDS),
    figure_format = "png|pdf",
    level = "|".join([a.replace("+", "\+") for a in NONE_BASELINE_ATTRIBUTE_LEVELS]),
    sample = "prior|posterior|prediction",
    layer= "ADM0|ADM1|ADM2",
min_version("7.8")
envvars:
    "ZENSUS_USER",
    "ZENSUS_PASSWORD"

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'peoples-preferences succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'peoples-preferences failed' {config[email]}")


def full_hierarchical_model_analysis(model: str, sample: str):
    return [
        f"build/results/models/hierarchical-{model}/{sample}/diagnostics/summary.csv",
        f"build/results/models/hierarchical-{model}/{sample}/intercept.png",
        f"build/results/models/hierarchical-{model}/{sample}/country-differences.png",
        f"build/results/models/hierarchical-{model}/{sample}/country-means.png",
        f"build/results/models/hierarchical-{model}/{sample}/individual-partworths.png",
        f"build/results/models/hierarchical-{model}/{sample}/unexplained-heterogeneity.png",
        f"build/results/models/hierarchical-{model}/{sample}/left-option.png",
        f"build/results/models/hierarchical-{model}/{sample}/varying/left-intercept.png"
    ] + [
        f"build/results/models/hierarchical-{model}/{sample}/varying/{level.replace(':', '_')}.png"
        for level in NONE_BASELINE_ATTRIBUTE_LEVELS
    ]

rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.html",
        "build/report.docx",
        "build/test-report.html",
        full_hierarchical_model_analysis(model="nocovariates-nocovariances", sample="prior"),
        full_hierarchical_model_analysis(model="nocovariates-nocovariances", sample="posterior"),
        full_hierarchical_model_analysis(model="covariates-nocovariances", sample="prior"),
        full_hierarchical_model_analysis(model="covariates-nocovariances", sample="posterior"),
        "build/results/models/hierarchical-nocovariates-nocovariances/poststratify/pop-means.png",
        "build/results/models/hierarchical-nocovariates-nocovariances/poststratify/diagnostics/summary.csv",
        "build/results/analysis/likert-items.png",
        "build/results/analysis/agreement-items.png",
        "build/results/analysis/gender.png",
        "build/results/analysis/area.png",
        "build/results/analysis/income.png",
        "build/results/analysis/education.png",
        "build/results/analysis/likert-items-by-country.png",
        "build/results/analysis/agreement-items-by-country.png"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--number-sections --embed-resources --standalone --mathml --to html5"
    elif suffix == "pdf":
        return "--number-sections --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        expand("build/results/analysis/respondent-stats-{country_id}.csv", country_id=COUNTRY_IDS),
        "build/results/analysis/experimental-design.png",
        "build/results/models/hierarchical-nocovariates-nocovariances/posterior/varying/left-intercept.png",
        "build/results/models/hierarchical-nocovariates-nocovariances/posterior/varying/TECHNOLOGY___Wind.png",
        "build/results/models/hierarchical-nocovariates-nocovariances/posterior/diagnostics/summary-brief.csv",
        "build/results/models/hierarchical-covariates-nocovariances/subgroups/max-subgroup-effect.png",
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


rule push:
    message: "Copy cluster build to {params.push_directory}."
    params: push_directory = Path(config["push-directory"] + "/" + datetime.date.today().isoformat()).as_posix()
    shell:
        """
        mkdir {params.push_directory}
        cp -R build/cluster/ {params.push_directory}
        cp -R build/logs/ {params.push_directory}/logs/
        """


rule dag:
     message: "Plot dependency graph of the workflow."
     conda: "envs/dag.yaml"
     shell:
         """
         snakemake --rulegraph > build/dag.dot
         dot -Tpdf -o build/dag.pdf build/dag.dot
         """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    message: "Run tests"
    input:
        test_dir = "tests",
        tests = map(str, Path("tests").glob("**/test_*.py")),
        national_conjoints = expand("build/data/{country_id}.feather", country_id=COUNTRY_IDS),
        conjoint = "build/data/conjoint.feather",
        conjoint_imputed = "build/data/conjoint-imputed.feather",
        covariate_model = "build/results/models/hierarchical-covariates-nocovariances/prior/inference-data.nc"
    params:
        config = config
    output: "build/test-report.html"
    conda: "envs/test.yaml"
    script: "tests/test_runner.py"

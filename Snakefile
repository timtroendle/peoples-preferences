import datetime
from pathlib import Path

from snakemake.utils import min_version

PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-secnos --citeproc"
COUNTRY_IDS = ["DEU", "POL", "PRT", "DNK"]
NONE_BASELINE_ATTRIBUTE_LEVELS = [
    'TECHNOLOGY:Open-field__PV', 'TECHNOLOGY:Wind',
    'LAND:1%', 'LAND:2%', 'LAND:4%', 'LAND:8%',
    'TRANSMISSION:+0%__.', 'TRANSMISSION:+25%__.', 'TRANSMISSION:+50%__.', 'TRANSMISSION:+75%__.',
    'SHARE_IMPORTS:10%', 'SHARE_IMPORTS:50%', 'SHARE_IMPORTS:90%',
    'PRICES:+15%', 'PRICES:+30%', 'PRICES:+45%', 'PRICES:+60%',
    'OWNERSHIP:Community', 'OWNERSHIP:Private'
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
    level = "|".join([a.replace("+", "\+") for a in NONE_BASELINE_ATTRIBUTE_LEVELS])
min_version("7.8")

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'peoples-preferences succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'peoples-preferences failed' {config[email]}")


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.html",
        "build/test-report.html",
        "build/results/models/multinomial-logit/pop-means.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/diagnostics/summary.csv",
        "build/results/models/hierarchical-nocovariates-nocovariances/pop-means.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/country-differences.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/country-means.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/individual-partworths.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/unexplained-heterogeneity.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/left-option.pdf",
        "build/results/models/hierarchical-nocovariates-nocovariances/varying/left-intercept.pdf",
        expand(
            "build/results/models/hierarchical-nocovariates-nocovariances/varying/{level}.pdf",
            level=NONE_BASELINE_ATTRIBUTE_LEVELS
        ),
        "build/results/likert-items.pdf",
        "build/results/agreement-items.pdf",
        "build/results/gender.pdf",
        "build/results/area.pdf",
        "build/results/income.pdf",
        "build/results/education.pdf",
        "build/results/likert-items-by-country.pdf",
        "build/results/agreement-items-by-country.pdf"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--number-sections --self-contained --to html5"
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
        "report/reset.css",
        "report/report.css",
        "report/fonts/KlinicSlabBook.otf",
        "report/fonts/KlinicSlabBookIt.otf",
        "report/fonts/KlinicSlabMedium.otf",
        "report/fonts/KlinicSlabMediumIt.otf",
        expand("build/results/respondent-stats-{country_id}.csv", country_id=COUNTRY_IDS),
        expand("build/results/clustering/4-with-country/umap-{feature}.png",
               feature=["RESPONDENT_COUNTRY", "Q3_GENDER", "Q6_AREA"]),
        "build/results/clustering/4-with-country/conditional-mm.png",
        "build/results/clustering/4-with-country/umap.png",
        "build/results/robustness/conditional-mm-choice-set.png",
        "build/results/robustness/conditional-mm-label.png",
        "build/results/robustness/design-validation.png",
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
        covariate_model = "build/results/models/hierarchical-nocovariates-nocovariances/inference-data.nc" # FIXME use covariate model
    params:
        config = config
    output: "build/test-report.html"
    conda: "envs/test.yaml"
    script: "tests/test_runner.py"

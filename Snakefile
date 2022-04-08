PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-secnos --citeproc"
COUNTRY_IDS = ["DEU", "POL", "PRT", "DNK"]

configfile: "config/default.yaml"
include: "rules/preprocess.smk"
include: "rules/analyse.smk"


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
        "build/test-report.html"


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
        expand("build/{country_id}/respondent-stats.csv", country_id=COUNTRY_IDS),
        expand("build/{country_id}/amce.png", country_id=COUNTRY_IDS),
        expand("build/{country_id}/mm.png", country_id=COUNTRY_IDS),
        "build/amce.png",
        "build/mm.png",
        "build/H1.png",
        "build/H2.png",
        "build/H4.png",
        "build/H5.png",
        "build/H6.png",
        "build/H7.png",
        "build/H8.png",
        "build/H9.png",
        "build/H11.png",
        "build/H16.png",
        "build/H17.png",
        "build/H18.png",
        "build/H19.png",
        "build/H20.png",
        "build/H21.png",
        "build/H22.png",
        "build/H23.png",
        "build/H24.png",
        "build/H25.png",
        "build/H26.png",
        "build/H27.png",
        "build/H28.png",
        "build/H29.png",
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
        national_conjoints = expand("build/{country_id}/conjoint.csv", country_id=COUNTRY_IDS)
    params:
        config = config
    output: "build/test-report.html"
    conda: "envs/test.yaml"
    script: "tests/test_runner.py"

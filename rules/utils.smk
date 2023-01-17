rule feather_to_csv:
    message: "Transform {wildcards.filename}.feather to csv."
    input:
        feather = "build/{filename}.feather"
    output: "build/{filename}.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/utils/feather_to_csv.py"


rule render_vega_lite:
    message: "Render Vega Lite spec {wildcards.filename}.json to pdf."
    input:
        json = "build/{path}/{filename}.vega.json"
    output:
        pdf = "build/{path}/{filename}.pdf"
    conda: "../envs/vega.yaml"
    # vl2pdf not usable because of https://github.com/queryverse/VegaLite.jl/issues/383
    shell: "vl2vg {input.json} | vg2pdf > {output.pdf}"

localrules: render_vega_lite_to_pdf, render_vega_lite_to_png


rule feather_to_csv:
    message: "Transform {wildcards.filename}.feather to csv."
    input:
        feather = "build/{filename}.feather"
    output: "build/{filename}.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/utils/feather_to_csv.py"


rule render_vega_lite_to_pdf:
    message: "Render Vega Lite spec {wildcards.filename}.json to pdf."
    input:
        json = "build/{path}/{filename}.vega.json"
    output:
        pdf = "build/{path}/{filename}.pdf"
    conda: "../envs/vega.yaml"
    # vl2pdf not usable because of https://github.com/queryverse/VegaLite.jl/issues/383
    shell: "vl2vg '{input.json}' | vg2pdf > '{output.pdf}'"


rule render_vega_lite_to_png:
    message: "Render Vega Lite spec {wildcards.filename}.json to png."
    input:
        json = "build/{path}/{filename}.vega.json"
    output:
        png = "build/{path}/{filename}.png"
    conda: "../envs/vega.yaml"
    # vl2png not usable because of https://github.com/queryverse/VegaLite.jl/issues/383
    shell: "vl2vg {input.json} | vg2png --scale 4 > {output.png}" # scale to >300dpi

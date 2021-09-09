rule scatter_data:
    input:
        expression = "results/quantification/{quant_pipeline}/{expression_unit}.tsv.gz",
        meta = rules.collect_metadata.output
    output:
        directory("results/scatter/{quant_pipeline}/{expression_unit}/")
    resources:
        time=240,
        mem=128000,
        cpus=2
    script:
        "../scripts/scatter_quantifications_by_te.R"

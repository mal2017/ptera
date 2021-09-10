rule bundle_expression:
    """
    Expects tsv files with filled headers (a la exporting via reader or vroom) in the
    format:

    feature\tLine1\tLine2 ...
    gene1\t100\t200 ...
    ...

    Additionally, genes only must be in FBgn ids.

    """
    input:
        expression = "results/quantification/{quant_pipeline}/{expression_unit}.tsv.gz",
        meta = rules.collect_metadata.output
    output:
        sqlite = "results/quantification/{quant_pipeline}/{expression_unit}.sqlite"
    resources:
        time=20,
        mem=20000,
        cpus=2
    script:
        "../scripts/bundle_quantifications.R"

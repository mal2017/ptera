rule bundle_expression:
    """
    Expects tsv files with filled headers (a la exporting via reader or vroom) in the
    format:

    feature\tLine1\tLine2 ...
    gene1\t100\t200 ...
    ...

    Additionally, genes must be in FBgn ids. TEs or other features can be anything.
    This allows for easy creation of separate tables in the resulting sqlite db.
    """
    input:
        male_expression = "results/quantification/{quant_pipeline}/male.{feature_level}.{expression_unit}.tsv.gz",
        female_expression = "results/quantification/{quant_pipeline}/female.{feature_level}.{expression_unit}.tsv.gz",
        meta = rules.collect_metadata.output
    output:
        sqlite = "results/quantification/{quant_pipeline}/{feature_level}.{expression_unit}.sqlite"
    resources:
        time=20,
        mem=20000,
        cpus=2
    script:
        "../scripts/bundle_quantifications.R"

rule linear_model:
    input:
        sqlite = lambda wc: expand("results/quantification/{p}/{u}.sqlite",p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"), u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS"))
    output:
        tidied = "results/linear_models/{model_id}/lm.tidy.tsv.gz",
        glanced = "results/linear_models/{model_id}/lm.glance.tsv.gz",
        augmented = "results/linear_models/{model_id}/lm.augment.tsv.gz",
        fits = "results/linear_models/{model_id}/lm.fits.rds"
    params:
        formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FORMULA"),
        split_by= lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SPLIT_BY"),
        transforms = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_TRANSFORM"),
        gene_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_GENE_FILTER"),
        te_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_TE_FILTER"),
    resources:
        time=360,
        mem=256000,
        cpus=2
    script:
        "../scripts/linear_model_v01.R"

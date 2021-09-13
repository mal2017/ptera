localrules: collect_chunked_linear_models, scatter_genes_for_lm

rule scatter_genes_for_lm:
    """
    Note: Escape brackets on wcs when using in conjunction w/ scattergather.
    """
    input:
        lambda wc: "results/quantification/{p}/{u}.tsv.gz".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"), u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS"))
    output:
        temp(expand("results/linear_models/{{model_id}}/chunk_{ch}",ch=[str(x).zfill(3) for x in range(0,config.get("LM_CHUNKS",80))]))
    params:
        chunks = config.get("LM_CHUNKS",80),
        prefix = lambda wc: "results/linear_models/{model_id}/chunk_".format(model_id=wc.model_id)
    shell:
        """
        zless {input} | cut -f 1 | tail -n+2 | grep "FBgn" | \
            split -e -d -a 3 \
            -n r/{params.chunks} - {params.prefix}
        """

rule chunked_linear_model:
    """
    Linear models were fit with R's `lm` function. See configfile and script for other details.
    """
    input:
        genes = "results/linear_models/{model_id}/chunk_{lmchunk}",
        sqlite = lambda wc: expand("results/quantification/{p}/{u}.sqlite",p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"), u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS"))
    output:
        tidied = temp("results/linear_models/{model_id}/chunk_{lmchunk}.tidy.tsv"),
        glanced = temp("results/linear_models/{model_id}/chunk_{lmchunk}.glance.tsv"),
        #augmented = temp("results/linear_models/{model_id}/chunk_{lmchunk}.augment.tsv"),
        fits = temp("results/linear_models/{model_id}/fits/chunk_{lmchunk}.fits.rds")
    params:
        formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FORMULA"),
        split_by= lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SPLIT_BY"),
        transforms = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_TRANSFORM"),
        gene_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_GENE_FILTER"),
        te_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_TE_FILTER"),
    resources:
        time=20,
        mem=12000,
        cpus=2
    script:
        "../scripts/linear_model_v01.R"

rule collect_chunked_linear_models:
    """
    Note: Escape brackets on wcs when using in conjunction w/ scattergather.
    """
    input:
        lambda wc: expand("results/linear_models/{{model_id}}/chunk_{ch}.{lmr}.tsv",ch = [str(x).zfill(3) for x in range(0,config.get("LM_CHUNKS",80))],lmr=wc.lmresult)
    output:
        "results/linear_models/{model_id}/lm.{lmresult}.tsv.gz"
    shell:
        "xsv cat rows -d '\t' {input} | gzip -c > {output}"


# rule linear_model:
#     input:
#         sqlite = lambda wc: expand("results/quantification/{p}/{u}.sqlite",p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"), u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS"))
#     output:
#         tidied = "results/linear_models/{model_id}/lm.tidy.tsv.gz",
#         glanced = "results/linear_models/{model_id}/lm.glance.tsv.gz",
#         augmented = "results/linear_models/{model_id}/lm.augment.tsv.gz",
#         fits = "results/linear_models/{model_id}/lm.fits.rds"
#     params:
#         formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FORMULA"),
#         split_by= lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SPLIT_BY"),
#         transforms = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_TRANSFORM"),
#         gene_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_GENE_FILTER"),
#         te_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_TE_FILTER"),
#     resources:
#         time=360,
#         mem=256000,
#         cpus=2
#     script:
#         "../scripts/linear_model_v01.R"

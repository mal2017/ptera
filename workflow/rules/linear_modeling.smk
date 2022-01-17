rule filter_transform_scale:
    """
    This rule takes parameters from a given linear model run section in config.yaml
    and filters, transforms, and scales accordingly. Optionally, it will correct the centered and scaled matrix
    by regressing out a specified number of prinicipal components.
    """
    input:
        sefile = lambda wc: "results/quantification/{p}/se.{fl}.rds".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"),fl=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_FEATURE_LEVEL")),
        mat =lambda wc: "results/quantification/{p}/{s}.{fl}.{c}.{u}.tsv.gz".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"),s=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SEX"),fl=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_FEATURE_LEVEL"), c=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_CNNORM"),u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS")),
    output:
        mat = "results/linear_models/{model_id}/expression.tsv.gz"
    params:
        filt = config.get("LM_FEATURE_FILT"),
        transforms = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_TRANSFORM"),
        scale = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_VARIABLE_SCALE"),
        pcs = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_CORRECT_N_PCS"),
    script:
        "../scripts/filter_transform_scale.R"


rule scatter_genes_for_lm:
    """
    Note: Escape brackets on wcs when using in conjunction w/ scattergather.

    https://www.unix.com/shell-programming-and-scripting/247089-creating-all-possible-bi-combinations-list-grep-awk.html
    """
    input:
        expression = rules.filter_transform_scale.output.mat
    output:
        temp(expand("results/linear_models/{{model_id}}/chunk_{ch}",ch=[str(x).zfill(4) for x in range(0,config.get("LM_CHUNKS"))]))
    params:
        split_invoc = "split -e -d -a 4 -n r/{ch} - results/linear_models/{{model_id}}/chunk_".format(ch= config.get("LM_CHUNKS"),)
    resources:
        time=240,
        mem=48000,
        cpus=1
    script:
        "../scripts/get_all_combos_tes_as_y.R"

rule chunked_linear_model:
    """
    Linear models were fit with R's `lm` function. See configfile and script for other details.
    """
    input:
        chunk = "results/linear_models/{model_id}/chunk_{lmchunk}",
        dat = rules.filter_transform_scale.output.mat,
        cd = rules.collect_metadata.output,
        ol = refs_wf("results/overlaps/overlaps.tsv.gz"),
        copies = wgs_wf("results/copies/copies.tsv") if config.get("INCL_COPY_ESTIMATION_IN_EXPORT") else rules.dummy_copies.output.feats
    output:
        tidy = temp("results/linear_models/{model_id}/chunk_{lmchunk}.tidy.tsv"),
    params:
        formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FORMULA"),
        alt_formula = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_ALT_FORMULA"),
        mads_filter = lambda wc: config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_MADS_FILTER"),
    conda:
        "../envs/tls_linmod_v1.yaml"
    resources:
        time=20,
        mem=10000,
        cpus=1
    script:
        "../scripts/linear_model_v03.R"

rule collect_chunked_linear_models:
    """
    Note: Escape brackets on wcs when using in conjunction w/ scattergather.
    """
    input:
        lambda wc: expand("results/linear_models/{{model_id}}/chunk_{ch}.{lmr}.tsv",ch = [str(x).zfill(4) for x in range(0,config.get("LM_CHUNKS",80))],lmr=wc.lmresult)
    output:
        "results/linear_models/{model_id}/lm.{lmresult}.tsv.gz"
    resources:
        time=60,
        mem=24000,
        cpus=1
    shell:
        "xsv cat rows -d '\t' {input} | tr ',' '\t' | gzip -c > {output}"

rule correct_tls_lm_coefs:
    input:
        coefs = "results/linear_models/{model_id}/lm.tidy.tsv.gz",
        expression = lambda wc: "results/quantification/{p}/{s}.{fl}.{c}.{u}.tsv.gz".format(p=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_PIPELINE"),s=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_SEX"),fl=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_FEATURE_LEVEL"), c=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_CNNORM"),u=config.get("LM_MODELS_TO_FIT").get(wc.model_id).get("LM_FIT_FOR_UNITS")),
    output:
        tsv = "results/linear_models/{model_id}/lm.tidy.corrected.tsv.gz"
    resources:
        time=60,
        mem=24000,
        cpus=1
    conda:
        "../envs/correct_coefs.yaml"
    script:
        "../scripts/correct_lm_coefs.R"

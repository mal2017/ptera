
rule make_r_objs:
    input:
        frydir = rules.alevin_fry_quant.output.dir
    output:
        sce = "results/downstream/raw_r_objs/{sample}.usa.sce.rds",
        #seur = "results/downstream/raw_r_objs/{sample}.usa.seur.rds",
        #h5ad = "results/downstream/initial_filtering/{sample}.usa.filt.h5ad",
    resources:
        time=5,
        mem=20000,
        cpus=2
    conda:
        "../envs/make_r_objs.yaml"
    script:
        "../scripts/import_alevin_fry.R"

rule initial_filtering:
    """
    Applies low count outlier, mitochondrial, and doublet filtering per sample.
    """
    input:
        sce = rules.make_r_objs.output.sce,
        #seur = rules.make_r_objs.output.seur,
    output:
        sce = "results/downstream/initial_filtering/{sample}.usa.filt.sce.rds",
        g_doublet = "results/downstream/figs/gg/{sample}_doublets.rds",
        g_miqc = "results/downstream/figs/gg/{sample}_mito.rds",
        dat_doublet = "results/downstream/figs/tsv/{sample}_doublets.tsv",
        dat_miqc = "results/downstream/figs/tsv/{sample}_mito.tsv",
        png_doublet ="results/downstream/figs/png/{sample}_doublets.png",
        png_miqc ="results/downstream/figs/png/{sample}_mito.png",
    resources:
        time=10,
        mem=20000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    script:
        "../scripts/initial_filtering.R"

rule sce_single_sample_dimred:
    """
    Finds HVGs, identifies number of PCs to use, and runs TSNE + UMAP per sample.
    """
    input:
        sce = rules.initial_filtering.output.sce
    output:
        sce = "results/downstream/single_sample_dimred/{sample}.usa.filt.dimred.sce.rds",
        dec = "results/downstream/single_sample_dimred/{sample}.usa.filt.dimred.dec.rds",
        g_hvg = "results/downstream/figs/gg/{sample}_hvg.rds",
        g_elbow = "results/downstream/figs/gg/{sample}_elbow.rds",
        dat_hvg = "results/downstream/figs/tsv/{sample}_hvg.tsv",
        dat_elbow = "results/downstream/figs/tsv/{sample}_elbow.tsv",
        png_hvg ="results/downstream/figs/png/{sample}_hvg.png",
        png_elbow ="results/downstream/figs/png/{sample}_elbow.png",
    resources:
        time=20,
        mem=20000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    script:
        "../scripts/sce_single_sample_dimred.R"

rule sce_integrate_by_tissue_and_sex:
    input:
        sces = lambda wc: expand("results/downstream/single_sample_dimred/{s}.usa.filt.dimred.sce.rds",s=[x.sample_name for x in pep.samples if (x.tissue == wc.tissue) and (x.sex == wc.sex) and (x.sample_name in SAMPLES)]),
        decs = lambda wc: expand("results/downstream/single_sample_dimred/{s}.usa.filt.dimred.dec.rds",s=[x.sample_name for x in pep.samples if (x.tissue == wc.tissue) and (x.sex == wc.sex) and (x.sample_name in SAMPLES)])
    output:
        sce = "results/downstream/integrated_by_tissue_and_sex/{tissue}_{sex}.sce.rds",
        ggp = "results/downstream/figs/gg/{tissue}_{sex}_batches.rds",
        png =  "results/downstream/figs/png/{tissue}_{sex}_batches.png",
        dat =  "results/downstream/figs/tsv/{tissue}_{sex}_batches.tsv",
    resources:
        time=40,
        mem=40000,
        cpus=2
    conda:
        "../envs/osca.yaml"
    script:
        "../scripts/sce_integrate_by_tissue_and_sex.R"

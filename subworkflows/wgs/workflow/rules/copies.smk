rule wgs_mosdepth:
    """
    Mosdepth is run to generate depth values in genomewide bins as well as bins
    across TEs. By defauly, reads matching sam flags 1796 are excluded from the calculation.
    Here I just allow multimappers because
    """
    input:
        cram = rules.wgs_samtools_markdup.output.cram,
        ref = rules.wgs_copy_resos_for_singularity.output.ref
    output:
        multiext('results/depth/{sample}','.mosdepth.global.dist.txt','.mosdepth.region.dist.txt', ".mosdepth.summary.txt", ".regions.bed.gz",".regions.bed.gz.csi")
    threads:
        12
    resources:
        time=20,
        mem=10000,
        cpus=12
    params:
        pfx = 'results/depth/{sample}',
        ws = config.get('MOSDEPTH_WINDOW_SIZE'),
        excl = config.get("MOSDEPTH_SAM_FLAG_EXCLUDE")
    singularity:
        "docker://quay.io/biocontainers/mosdepth:0.3.2--h01d7912_0"
    priority: 3
    shell:
        """
        mosdepth -n -t {threads} --by {params.ws} -f {input.ref} \
            --use-median {params.pfx} {input.cram} -F {params.excl}
        """

rule wgs_estimate_te_copies:
    input:
        cov = 'results/depth/{sample}.mosdepth.summary.txt',
        fl = main_wf("results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2symbol.tsv")
    output:
        csv = 'results/copies/by_strain/{sample}.csv'
    resources:
        time=30,
        mem=20000,
        cpus=2
    priority: 3
    script:
        "../scripts/copies.R"

rule wgs_collect_te_copy_estimates:
    input:
        expand('results/copies/by_strain/{s}.csv',s=SAMPLES)
    output:
        temp('results/copies/tmp.all.csv')
    shell:
        "xsv cat rows {input} > {output}"

rule wgs_add_strain_to_copy_estimates:
    input:
        "config/sample_table.csv",
        rules.wgs_collect_te_copy_estimates.output
    output:
        "results/copies/copies.tsv"
    shell:
        """
        xsv join --right sample_name {input[0]} sample_name {input[1]} | \
        xsv select '!1' | tr ',' '\t' > {output}
        """

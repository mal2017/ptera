rule fasterq_dump:
    """
    We downloaded reads from the NCBI Sequence Read Archive using the
    fasterq-dump utility.
    """
    output:
        temp(directory("results/reads/dump/{sample}/{experiment}/{run}/"))
    threads:
        22
    resources:
        time=240,
        mem=config.get("FASTERQDUMP_MEM","8000"),
        cpus=1
    priority: 1
    conda:
        "../envs/sratools.yaml"
    shell:
        """
        mkdir -p {output} &&
        fasterq-dump --mem {resources.mem}MB -s -S --include-technical -e {threads} -O {output}/ {wildcards.run}
        """

rule concat_runs:
    """
    We concatenated single-end fastqs from the same SRX accession (same library).
    """
    input:
        lambda wc: flatten([expand("results/reads/dump/{s}/{e}/{r}/",s=wc.sample,e=x,r=get_runs(wc.sample,x)) for x in get_experiments(wc.sample)])
        #lambda wc: expand("results/reads/dump/{s}/{e}/{r}/",s=wc.sample,e=wc.experiment,r=get_accession(wc, "Run"))
    output:
        temp("results/reads/concat/{sample}.fastq")
    params:
        fqs = lambda wc: flatten([expand("results/reads/dump/{s}/{e}/{r}/{r}.fastq",s=wc.sample,e=x,r=get_runs(wc.sample,x)) for x in get_experiments(wc.sample)])
    threads:
        1
    resources:
        time=120,
        mem=24000,
        cpus=1
    priority: 2
    shell:
        "cat {params.fqs} > {output}"

rule concat_pe_runs:
    """
    We concatenated single-end fastqs from the same SRX accession (same library).
    """
    input:
        lambda wc: flatten([expand("results/reads/dump/{s}/{e}/{r}/",s=wc.sample,e=x,r=get_runs(wc.sample,x)) for x in get_experiments(wc.sample)])
        #lambda wc: expand("results/reads/dump/{s}/{e}/{r}/",s=wc.sample,e=wc.experiment,r=get_accession(wc, "Run"))
    output:
        r1=temp("results/reads/concat/{sample}_r1.fastq"),
        r2=temp("results/reads/concat/{sample}_r2.fastq")
    params:
        fqs1 = lambda wc: flatten([expand("results/reads/dump/{s}/{e}/{r}/{r}_{m}.fastq",s=wc.sample,e=x,r=get_runs(wc.sample,x),m=1) for x in get_experiments(wc.sample)]),
        fqs2 = lambda wc: flatten([expand("results/reads/dump/{s}/{e}/{r}/{r}_{m}.fastq",s=wc.sample,e=x,r=get_runs(wc.sample,x),m=2) for x in get_experiments(wc.sample)])
    threads:
        1
    resources:
        time=60,
        mem=20000,
        cpus=1
    priority: 2
    shell:
        """
        cat {params.fqs1} > {output.r1} &&
        cat {params.fqs2} > {output.r2}
        """

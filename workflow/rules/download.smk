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
        time=60,
        mem=config.get("FASTERQDUMP_MEM","8000"),
        cpus=1
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
        time=60,
        mem=20000,
        cpus=1
    shell:
        "cat {params.fqs} > {output}"

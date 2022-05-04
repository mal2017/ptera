rule fastp_trim_se:
    """
    We trimmed reads from the combined fastqs for each library with fastp with
    default parameters.
    """
    input:
        r1=rules.concat_runs.output
    output:
        r1 = temp("results/reads/fastp-se/{sample}_r1.trimmed.fq.gz"),
        html = "results/reads/fastp-se/{sample}_fastp.html",
        json = "results/reads/fastp-se/{sample}_fastp.json"
    threads:
        6
    resources:
        time=60,
        mem=64000,
        cpus=4
    singularity:
        "docker://quay.io/biocontainers/fastp:0.22.0--h2e03b76_0"
    priority: 3
    shell:
        """
        fastp --in1 {input.r1} \
        --out1 {output.r1} \
        -j {output.json} -h {output.html} \
        -w {threads} -L -R {wildcards.sample}_fastp
        """

rule fastp_trim_pe:
    """
    We trimmed reads from the combined fastqs for each library with fastp with
    default parameters.
    """
    input:
        r1="results/reads/concat/{sample}_r1.fastq",
        r2="results/reads/concat/{sample}_r2.fastq"
    output:
        r1 = temp("results/reads/fastp-pe/{sample}_r1.trimmed.fq.gz"),
        r2 = temp("results/reads/fastp-pe/{sample}_r2.trimmed.fq.gz"),
        html = "results/reads/fastp-pe/{sample}_fastp.html",
        json = "results/reads/fastp-pe/{sample}_fastp.json"
    threads:
        6
    resources:
        time=60,
        mem=64000,
        cpus=4
    singularity:
        "docker://quay.io/biocontainers/fastp:0.22.0--h2e03b76_0"
    priority: 3
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
        --out1 {output.r1} --out2 {output.r2} \
        -j {output.json} -h {output.html} \
        -w {threads} -L -R {wildcards.sample}_fastp
        """

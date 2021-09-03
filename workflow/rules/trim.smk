rule fastp_trim_se:
    """
    We trimmed reads from the combined fastqs for each library with fastp with
    default parameters.
    """
    input:
        r1=rules.concat_runs.output
    output:
        r1 = temp("results/reads/fastp/{sample}_r1.trimmed.fq.gz"),
        html = "results/reads/fastp/{sample}_fastp.html",
        json = "results/reads/fastp/{sample}_fastp.json"
    threads:
        12
    resources:
        time=60,
        mem=20000,
        cpus=12
    singularity:
        "docker://quay.io/biocontainers/fastp:0.22.0--h2e03b76_0"
    priority: 40
    shell:
        """
        fastp --in1 {input.r1} \
        --out1 {output.r1} \
        -j {output.json} -h {output.html} \
        -w {threads} -L -R {wildcards.sample}_fastp
        """

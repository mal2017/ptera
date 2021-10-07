localrules: wgs_copy_resos_for_singularity

rule wgs_copy_resos_for_singularity:
    """
    This moves the resource to a suitable spot given singularity's filesystem
    mount point for the job. Additionally removes descriptions from fasta headers,
    which cause bwa-mem2 to throw unhelpful errors during alignment.

    https://www.biostars.org/p/17654/
    https://www.biostars.org/p/228886/
    """
    input:
        ref = "../../" + config.get("COMBINED_TE_GENOME_FA")
    output:
        ref = "results/mapping/index/ref.fasta.gz",
        fai = "results/mapping/index/ref.fasta.gz.fai",
        gzi = "results/mapping/index/ref.fasta.gz.gzi"
    shell:
        """
        bgzip -d -c {input.ref} | \
            awk '{{print $1}}' |
            bgzip -c > {output.ref} &&
        samtools faidx {output.ref}
        """

rule wgs_bwa_mem2_index:
    input:
        rules.wgs_copy_resos_for_singularity.output.ref
    output:
        expand("results/mapping/index/idx.{suf}",suf=["0123","amb","ann","bwt.2bit.64","pac"])
    params:
        pfx = "results/mapping/index/idx"
    threads:
        24
    resources:
        time=60,
        mem=48000,
        cpus=24
    singularity:
        "docker://quay.io/biocontainers/bwa-mem2:2.2.1--h9a82719_1"
    shell:
        """
        bwa-mem2 index -p {params.pfx} {input}
        """

rule wgs_bwa_mem2_align:
    input:
        r1 = rules.fastp_trim_pe.output.r1,
        r2 = rules.fastp_trim_pe.output.r2,
        idx = rules.wgs_bwa_mem2_index.output
    output:
        temp("results/mapping/aligned/{sample}.sam")
    threads:
        24
    params:
        idx = "results/mapping/index/idx"
    resources:
        time=480,
        mem=128000,
        cpus=24
    priority: 2
    singularity:
        "docker://quay.io/biocontainers/bwa-mem2:2.2.1--h9a82719_1"
    shell:
        """
        bwa-mem2 mem -t {threads} {params.idx} {input.r1} {input.r2} > {output}
        """

rule wgs_samtools_fixmate:
    input:
        sam = rules.wgs_bwa_mem2_align.output,
        ref = rules.wgs_copy_resos_for_singularity.output.ref
    output:
        cram = temp("results/mapping/aligned/{sample}.fixm.cram"),
    threads:
        8
    resources:
        time=20,
        mem=20000,
        cpus=8
    priority: 2
    singularity:
        "docker://quay.io/biocontainers/samtools:1.13--h8c37831_0"
    shell:
        """
        samtools fixmate -@ {threads} -m -O CRAM \
            --reference {input.ref} \
            {input.sam} {output.cram}
        """

rule wgs_samtools_sort:
    input:
        rules.wgs_samtools_fixmate.output.cram
    output:
        cram = temp("results/mapping/aligned/{sample}.srt.cram"),
        crai = temp("results/mapping/aligned/{sample}.srt.cram.crai")
    singularity:
        "docker://quay.io/biocontainers/samtools:1.13--h8c37831_0"
    resources:
        time=240,
        mem=20000,
        cpus=12
    threads:
        12
    priority: 2
    shell:
        """
        samtools sort -@ {threads} -m 1G {input} -o {output.cram} &&
        samtools index -@ {threads} {output.cram}
        """

rule wgs_samtools_markdup:
    input:
        cram = rules.wgs_samtools_sort.output.cram,
        ref = rules.wgs_copy_resos_for_singularity.output.ref
    output:
        cram = "results/mapping/aligned/{sample}.markdup.cram",
        crai = "results/mapping/aligned/{sample}.markdup.cram.crai"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.13--h8c37831_0"
    resources:
        time=240,
        mem=20000,
        cpus=12
    threads:
        12
    priority: 2
    shell:
        """
        samtools markdup -@ {threads} --reference {input.ref} -O CRAM \
            {input.cram} {output.cram} &&
        samtools index -@ {threads} {output.cram}
        """

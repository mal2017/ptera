localrules: copy_bwa_mem2_indices_to_mount

rule copy_bwa_mem2_indices_to_mount:
    """
    Copy to a dir that singularity will see.
    """
    input:
        a = refs_wf(expand("results/indices/bwa_mem2/idx.{suf}",suf=["0123","amb","ann","bwt.2bit.64","pac"])),
        ref = refs_wf("results/indices/bwa_mem2/ref.fasta.gz"),
        fai = refs_wf("results/indices/bwa_mem2/ref.fasta.gz.fai"),
        gzi = refs_wf("results/indices/bwa_mem2/ref.fasta.gz.gzi"),
    output:
        idx = expand("results/indices/bwa_mem2/idx.{suf}",suf=["0123","amb","ann","bwt.2bit.64","pac"]),
        ref = "results/indices/bwa_mem2/ref.fasta.gz",
        fai = "results/indices/bwa_mem2/ref.fasta.gz.fai",
        gzi = "results/indices/bwa_mem2/ref.fasta.gz.gzi",
    params:
        idx = "results/indices/bwa_mem2"
    shell:
        """
        mkdir -p {params.idx} &&
        cp {input.a} {params.idx}
        cp {input.ref} {output.ref}
        cp {input.fai} {output.fai}
        cp {input.gzi} {output.gzi}
        """

rule dna_bwa_mem2_align:
    input:
        reads = lambda wc: [rules.fastp_trim_pe.output.r1,rules.fastp_trim_pe.output.r2] if is_paired_end(wc.sample) else rules.fastp_trim_se.output.r1,
        idx = rules.copy_bwa_mem2_indices_to_mount.output.idx
    output:
        temp("results/mapping/dna_bwa_mem2/{sample}.sam")
    threads:
        24
    params:
        idx = "results/indices/bwa_mem2/idx"
    resources:
        time=480,
        mem=128000,
        cpus=24
    priority: 2
    singularity:
        "docker://quay.io/biocontainers/bwa-mem2:2.2.1--h9a82719_1"
    shell:
        """
        bwa-mem2 mem -t {threads} {params.idx} {input.reads} > {output}
        """

rule dna_samtools_fixmate:
    input:
        sam = rules.dna_bwa_mem2_align.output,
        ref = rules.copy_bwa_mem2_indices_to_mount.output.ref
    output:
        cram = temp("results/mapping/dna_bwa_mem2/{sample}.fixm.cram"),
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

rule dna_samtools_sort:
    input:
        rules.dna_samtools_fixmate.output.cram
    output:
        cram = temp("results/mapping/dna_bwa_mem2/{sample}.srt.cram"),
        crai = temp("results/mapping/dna_bwa_mem2/{sample}.srt.cram.crai")
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

rule dna_samtools_markdup:
    input:
        cram = rules.dna_samtools_sort.output.cram,
        ref = rules.copy_bwa_mem2_indices_to_mount.output.ref
    output:
        cram = "results/mapping/dna_bwa_mem2/{sample}.markdup.cram",
        crai = "results/mapping/dna_bwa_mem2/{sample}.markdup.cram.crai"
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

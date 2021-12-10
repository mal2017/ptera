localrules: copy_for_rpm

rule copy_for_rpm:
    input:
        repeats =  rules.get_resources_in_mount.output.CONSENSUS_TE_FASTA,
        fasta = rules.get_resources_in_mount.output.GENOME_FASTA
    output:
        fasta = temp('results/repeatmasker/genome.fasta'),
        repeats = temp('results/repeatmasker/repeats.fasta'),
    shell:
        """
        bgzip -d -c {input.repeats} | \
            awk '{{print $1}}' > {output.repeats}
        bgzip -d -c {input.fasta} | \
            awk '{{print $1}}' > {output.fasta}
        """

rule repeatmasker:
    input:
        repeats = rules.copy_for_rpm.output.repeats,
        fasta = rules.copy_for_rpm.output.fasta
    output:
        'results/repeatmasker/genome.fasta.cat.gz',
        'results/repeatmasker/genome.fasta.masked',
        'results/repeatmasker/genome.fasta.out',
        'results/repeatmasker/genome.fasta.ori.out',
        'results/repeatmasker/genome.fasta.out.gff',
        'results/repeatmasker/genome.fasta.tbl',
    threads:
        24
    params:
        dir = "results/repeatmasker/"
    singularity:
        "docker://quay.io/biocontainers/repeatmasker:4.1.2-p1"
    resources:
        mem=128000,
        cpus=24,
        time=240,
    shell:
        """
        RepeatMasker -e ncbi -pa {threads} -s \
            -lib {input.repeats} -gff \
            -no_is -nolow -dir {params.dir} \
            {input.fasta}
        """

rule bedops_parse_repeatmasker:
    input:
        rpmskr = rules.repeatmasker.output[2],
    output:
        "results/repeatmasker/reference_insertions.bed"
    singularity:
        "docker://quay.io/biocontainers/bedops:2.4.39--h7d875b9_1"
    resources:
        mem=2000,
        cpus=2,
        time=5,
    shell:
        """
        rmsk2bed < {input.rpmskr} | cut -f 1,2,3,4,5 > {output}
        """

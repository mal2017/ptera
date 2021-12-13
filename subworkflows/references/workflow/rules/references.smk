rule make_transcripts_and_consensus_tes_fasta:
    """
    We generated the combined transcriptome reference by concatenating the set of transcript
    sequences and the set of consensus TE sequences.
    """
    input:
        tes = "results/" + config.get("CONSENSUS_TE_FASTA"),
        txs = "results/" + config.get("FULL_TRANSCRIPT_FASTA"),
        trna = "results/" + config.get("TRNA_FASTA"),
        ncrna = "results/" + config.get("NCRNA_FASTA"),
        miscrna = "results/" + config.get("MISCRNA_FASTA"),
    output:
        fa = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.fasta.gz",
        dummy_decoy = touch("results/references/transcripts_and_consensus_tes/dummy_decoy.txt")
    shell:
        "cat {input.tes} {input.txs} {input.trna} {input.ncrna} {input.miscrna} > {output.fa}"

rule make_transcripts_and_consensus_tes_gtf:
    """
    We generated the combined transcriptome reference by concatenating the set of transcript
    sequences and the set of consensus TE sequences.
    """
    input:
        te_fasta = "results/" + config.get("CONSENSUS_TE_FASTA"),
        host_gtf = "results/" + config.get("TRANSCRIPTOME_GTF")
    output:
        gtf = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf",
    singularity:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.52.0--r41hd029910_0"
    script:
        "../scripts/make_transcripts_and_consensus_tes_gtf.R"

rule make_transcripts_and_consensus_tes_tx2gene:
    """
    We generated the combined transcriptome reference by concatenating the set of transcript
    sequences and the set of consensus TE sequences.
    """
    input:
        gtf = rules.make_transcripts_and_consensus_tes_gtf.output.gtf
    output:
        tx2id = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2id.tsv",
        tx2symbol = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2symbol.tsv",
        tx2txsymbol = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2txsymbol.tsv",
    singularity:
        "docker://quay.io/biocontainers/bioconductor-rtracklayer:1.52.0--r41hd029910_0"
    script:
        "../scripts/make_transcripts_and_consensus_tes_tx2gene.R"

localrules: dna_ref_prep

rule dna_ref_prep:
    """
    This moves the resource to a suitable spot given singularity's filesystem
    mount point for the job. Additionally removes descriptions from fasta headers,
    which cause bwa-mem2 to throw unhelpful errors during alignment.

    https://www.biostars.org/p/17654/
    https://www.biostars.org/p/228886/
    """
    input:
        ref = "results/" + config.get("COMBINED_TE_GENOME_FA")
    output:
        ref = "results/indices/bwa_mem2/ref.fasta.gz",
        fai = "results/indices/bwa_mem2/ref.fasta.gz.fai",
        gzi = "results/indices/bwa_mem2/ref.fasta.gz.gzi"
    shell:
        """
        bgzip -d -c {input.ref} | \
            awk '{{print $1}}' |
            bgzip -c > {output.ref} &&
        samtools faidx {output.ref}
        """

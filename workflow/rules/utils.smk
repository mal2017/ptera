# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def get_runs(sample,experiment):
    """
    Returns the run accessions for a given sample (biosample).
    """
    df = pep.subsample_table[pep.subsample_table.sample_name == sample]
    return list(set(df[df.Experiment == experiment].Run))

def get_experiments(sample):
    """
    Returns the experiment accessions for a given sample (biosample).
    """
    return list(set(pep.subsample_table[pep.subsample_table.sample_name == sample].Experiment))

flatten = lambda t: [item for sublist in t for item in sublist]

# ------------------------------------------------------------------------------
# Utility Rules
# ------------------------------------------------------------------------------

localrules: make_transcripts_and_consensus_tes_fasta, get_pipeline_info

rule make_transcripts_and_consensus_tes_fasta:
    """
    We generated the combined transcriptome reference by concatenating the set of transcript
    sequences and the set of consensus TE sequences.
    """
    input:
        tes = config.get("CONSENSUS_TE_FASTA"),
        txs = config.get("FULL_TRANSCRIPT_FASTA")
    output:
        fa = "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.fasta.gz",
        dummy_decoy = touch("results/references/transcripts_and_consensus_tes/dummy_decoy.txt")
    shell:
        "cat {input.tes} {input.txs} > {output.fa}"

rule make_transcripts_and_consensus_tes_gtf:
    """
    We generated the combined transcriptome reference by concatenating the set of transcript
    sequences and the set of consensus TE sequences.
    """
    input:
        te_fasta = config.get("CONSENSUS_TE_FASTA"),
        host_gtf = config.get("TRANSCRIPTOME_GTF")
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

rule get_pipeline_info:
    output:
        "results/meta/pipeline_meta.txt"
    shell:
        """
        git rev-parse HEAD > {output} &&
        git config --get remote.origin.url >> {output}
        """

rule se_export_txt:
    """
    Generic rule for exporting raw or normalized expression from a serialized
    SummarizedExperiment object.
    """
    input:
        se="results/quantification/{quant_pipeline}/se.{feature_level}.rds"
    output:
        txt="results/quantification/{quant_pipeline}/{sex}.{feature_level}.{expression_unit}.tsv.gz"
    resources:
        time=240,
        mem=20000,
        cpus=1
    params:
        DESEQ2_FITTYPE = config.get("DESEQ2_FITTYPE"),
    script:
        "../scripts/se_export_txt.R"

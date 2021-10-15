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

rule dummy_copies:
    """
    Rule generates a strain x feature matrix populated with 1s. This allows
    this pipeline to run when WGS is available and when WGS is not available.
    As with situations when WGS is available, samples are joined by their strain,
    so this must be accurately included in the sample_table.csv in the config
    directory.

    Format required is below, but only strain and est.copies need to be populated.
    Strain, sample_name, sequence, length, bases, median.cov, est.copies

    Note that this includes features at the tx and gene level. Similarly, it includes
    TE scaffold names and TE feature names, all set to have est.copies = 1.
    """
    input:
        samples = "config/sample_table.csv",
        feats = rules.make_transcripts_and_consensus_tes_tx2gene.output.tx2id,
    output:
        feats = temp("results/dummy-copies.tsv"),
    script:
        "../scripts/dummy-copies.R"

rule se_export_txt:
    """
    Generic rule for exporting raw or normalized expression from a serialized
    SummarizedExperiment object.
    """
    input:
        se="results/quantification/{quant_pipeline}/se.{feature_level}.rds",
        copies = wgs_wf("results/copies/copies.tsv") if config.get("INCL_COPY_ESTIMATION_IN_EXPORT") else rules.dummy_copies.output.feats
    output:
        txt="results/quantification/{quant_pipeline}/{sex}.{feature_level}.{cnnorm}.{expression_unit}.tsv.gz"
    resources:
        time=240,
        mem=20000,
        cpus=1
    params:
        DESEQ2_FITTYPE = config.get("DESEQ2_FITTYPE"),
        copy_adjustment = lambda wc: True if wc.cnnorm == "per_est_copy" else False
    script:
        "../scripts/se_export_txt.R"

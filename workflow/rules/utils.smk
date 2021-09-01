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

localrules: make_transcripts_and_consensus_tes_fasta

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

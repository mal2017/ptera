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

def is_paired_end(sample):
    """
    Returns the libary layout for a given sample. Assumes that all libraries for a given sample
    use the same layout, otherwise returns a string.
    """
    layouts_represented = set(pep.subsample_table[pep.subsample_table.sample_name == sample].LibraryLayout)
    if len(layouts_represented) > 1:
        return "ERROR! make sure each sample has only paired end or single end libraries."
    return all([x == "PAIRED" for x in layouts_represented])

def get_background_sample(sample):
    """
    Returns the 'input' or WCE sample for IP type experiments or others
    where a background sample is required for analysis. Operates on the level of samples, not
    subsamples. Fails/returns NA if the column 'input' doesn't exist in the sample_table.
    """
    x=pep.sample_table[pep.sample_table.sample_name == sample].input[0]
    return x

flatten = lambda t: [item for sublist in t for item in sublist]

# ------------------------------------------------------------------------------
# Utility Rules
# ------------------------------------------------------------------------------

localrules: get_pipeline_info


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
        feats = refs_wf("results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2id.tsv")
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
        se="results/quantification/{quant_pipeline}/se.{feature_level}.{quant_rep}.rds",
    output:
        txt="results/quantification/{quant_pipeline}/{sex}.{feature_level}.{cnnorm}.{expression_unit}.{quant_rep}.tsv.gz"
    resources:
        time=240,
        mem=20000,
        cpus=1
    params:
        DESEQ2_FITTYPE = config.get("DESEQ2_FITTYPE"),
    script:
        "../scripts/se_export_txt.R"

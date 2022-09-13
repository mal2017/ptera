library(tximeta)
library(BiocFileCache)


# -----------------------------------------------------------------------------
# make tximeta json record
# -----------------------------------------------------------------------------

#jsonFile <- "results/quantification/vanilla_salmon_tes_transcripts/tximeta.json"
jsonFile <- snakemake@output[["json"]]

#indexDir <- "results/quantification/vanilla_salmon_tes_transcripts/index/"
indexDir <- snakemake@input[["idx"]]

#flybase_release <- "FB2021_04"
flybase_release <- snakemake@params[["flybase_release"]]

#genome_version <- "r6.41"
genome_version <- snakemake@params[["genome_version"]]

#tx_fasta <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.fasta.gz"
#tx_gtf <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf"
tx_fasta <- snakemake@input[["fasta"]]
tx_gtf <- snakemake@input[["gtf"]]

makeLinkedTxome(indexDir = indexDir, source = "Flybase + repeats",
                organism = "Drosophila melanogaster",
                jsonFile = jsonFile,
                release = flybase_release,
                fasta = tx_fasta,
                gtf = tx_gtf,
                genome = genome_version)

# remove cache
if (interactive()) {
  bfcloc <- getTximetaBFC()
} else {
  bfcloc <- tempdir()
}
bfc <- BiocFileCache(bfcloc)
bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
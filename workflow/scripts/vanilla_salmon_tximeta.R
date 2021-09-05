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

# # remove cache
# if (interactive()) {
#   bfcloc <- getTximetaBFC()
# } else {
#   bfcloc <- tempdir()
# }
# bfc <- BiocFileCache(bfcloc)
# bfcremove(bfc, bfcquery(bfc, "linkedTxomeTbl")$rid)
#
# loadLinkedTxome(jsonFile)

# -----------------------------------------------------------------------------
# make tximeta object
# -----------------------------------------------------------------------------

# get sample table
#samples_fl <- "config/sample_table.csv"
samples_fl <- snakemake@input[["samples"]]
samples <- read.csv(samples_fl,header = T)

# name the quant files so reordering is possible
#files <- Sys.glob("results/quantification/vanilla_salmon_tes_transcripts/quant/*/quant.sf")
#terminus_files <- Sys.glob("results/quantification/vanilla_salmon_tes_transcripts/terminus/*/quant.sf")
files <- snakemake@input[["salmon_files"]]
#terminus_files <- snakemake@params[["terminus_files"]]
names(files) <- gsub("\\/+quant.sf","",x=gsub(".+quant\\/+","",files))
#names(terminus_files) <- gsub("\\/+quant.sf","",x=gsub(".+terminus\\/+","",terminus_files))

# useful for testing a smaller subset, also good as foolproofing
# in case some samples aren't quantified and this script is run manually
samples <- samples[samples$sample_name %in% names(files),]

# reorder to match sample table
files <- files[samples$sample_name]
#terminus_files <- terminus_files[samples$sample_name]

samples$names <- samples$sample_name
samples$files <- files
#samples_terminus <- samples
#samples_terminus$files <- samples_terminus$files <- terminus_files

## params for import
#counts_from_abundance_salmon <- "lengthScaledTPM"
#counts_from_abundance_salmon_txout <- "dtuScaledTPM"
#counts_from_abundance_terminus <- "lengthScaledTPM"

counts_from_abundance_salmon <- snakemake@params[["counts_from_abundance_salmon"]]
counts_from_abundance_salmon_txout <- snakemake@params[["counts_from_abundance_salmon_txout"]]
#counts_from_abundance_terminus <- snakemake@params[["counts_from_abundance_terminus"]]

#tx2gene_fl <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2symbol.tsv"
#tx2feature_fl <- "results/quantification/vanilla_salmon_tes_transcripts/terminus.tx2group.tsv"
tx2gene_fl <- snakemake@input[["tx2gene"]]
#tx2feature_fl <- snakemake@input[["tx2feature"]]

# get the conversion for transcripts
tx2gene <- read.table(tx2gene_fl,header = T)

# get the conversion for the renamed terminus grps w/ symbols in the name
#tx2feature <- read.table(tx2feature_fl,header = T)
#rownames(tx2feature) <- tx2feature$TXNAME
#tx2groupid_symbols <- tx2feature[,c("TXNAME","GROUPID_SYMBOLS")]

salmon_tx_se <- tximeta(samples,
                     type="salmon",
                     #customMetaInfo = file.path("..","..",basename(jsonFile)),
                     useHub = F,
                     skipMeta = F,
                     tx2gene=tx2gene,
                     txOut = T,
                     cleanDuplicateTxps = T,
                     markDuplicateTxps = T,
                     countsFromAbundance = counts_from_abundance_salmon_txout)

salmon_se <- summarizeToGene(salmon_tx_se,countsFromAbundance = counts_from_abundance_salmon)

# terminus plays weird with tximeta because the feature names don't match with "genes"
# in the gtf. So this is vanilla tximport but straight to a SummarizedExperiment with
# not as much metadata as the regular salmon import.
#terminus_se <- tximeta(samples_terminus,
#                        type="salmon",
#                        #customMetaInfo = file.path("..","..",basename(jsonFile)),
#                        useHub = F,
#                        skipMeta = T,
#                        txOut = T,
#                        cleanDuplicateTxps = T,
#                        markDuplicateTxps = T,
#                        countsFromAbundance = counts_from_abundance_terminus)

# --------------------------------------------------------------------------------------
# extra metadata
# --------------------------------------------------------------------------------------

# add pipeline hash
pipeline_meta <- readLines(snakemake@input[["pipeline_meta"]],n = 2)
pipeline_meta <-list(hash=pipeline_meta[1], url=pipeline_meta[2])

salmon_tx_se@metadata$pipeline_info <- pipeline_meta
salmon_se@metadata$pipeline_info <- pipeline_meta
#terminus_se@metadata$pipeline_info <- pipeline_meta

# id conversion for terminus
#terminus_se@elementMetadata <- DataFrame(tx2feature[terminus_se@NAMES,])

# TODO: add DGRP line metadata here directly, ie wolbachia, etc

# export
saveRDS(salmon_tx_se, snakemake@output[["salmon_tx"]])
saveRDS(salmon_se, snakemake@output[["salmon"]])
#saveRDS(terminus_se,snakemake@output[["terminus"]])

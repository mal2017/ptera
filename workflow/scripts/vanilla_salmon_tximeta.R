library(tximeta)
library(BiocFileCache)

jsonFile <- snakemake@input[["json"]]

loadLinkedTxome(jsonFile)



# -----------------------------------------------------------------------------
# make tximeta object
# -----------------------------------------------------------------------------

tx_fasta <- snakemake@input[["fasta"]]
tx_gtf <- snakemake@input[["gtf"]]

# get sample table
#samples_fl <- "config/sample_table.csv"
samples_fl <- snakemake@input[["samples"]]
samples <- read.csv(samples_fl,header = T)

# name the quant files so reordering is possible
#files <- Sys.glob("results/quantification/vanilla_salmon_tes_transcripts/quant/*/0/quant.sf")
files <- snakemake@input[["salmon_files"]]

# from before replicated salmon runs
#names(files) <- gsub("\\/+quant.sf","",x=gsub(".+quant\\/+","",files))
names(files) <- gsub("\\/\\d+\\/+quant.sf","",x=gsub(".+quant\\/+","",files))


# useful for testing a smaller subset, also good as foolproofing
# in case some samples aren't quantified and this script is run manually
samples <- samples[samples$sample_name %in% names(files),]

# reorder to match sample table
files <- files[samples$sample_name]

samples$names <- samples$sample_name
samples$files <- files

## params for import
#counts_from_abundance_salmon <- "lengthScaledTPM"
#counts_from_abundance_salmon_txout <- "dtuScaledTPM"

counts_from_abundance_salmon <- snakemake@params[["counts_from_abundance_salmon"]]
counts_from_abundance_salmon_txout <- snakemake@params[["counts_from_abundance_salmon_txout"]]

#tx2gene_fl <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2symbol.tsv"
#tx2feature_fl <- "results/quantification/vanilla_salmon_tes_transcripts/terminus.tx2group.tsv"
tx2gene_fl <- snakemake@input[["tx2gene"]]
#tx2feature_fl <- snakemake@input[["tx2feature"]]

# get the conversion for transcripts
tx2gene <- read.table(tx2gene_fl,header = T)

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

gtf <- rtracklayer::import(tx_gtf)
#summary(gtf$type)

allowed.types <- c("mRNA")

allowed.features <- unique(gtf[gtf$type %in% allowed.types]$transcript_id)

# strictly speaking, not necessary as these types aren't quantified in the current pipeline.
# assumes that TEs are type "mRNA" or at least one of the other allowed types above.
salmon_tx_se <- salmon_tx_se[rownames(salmon_tx_se) %in% allowed.features,]

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

# ------------------------------------------
# other per-sample annotation info
# ------------------------------------------
# TODO: add DGRP line metadata here directly, ie wolbachia, etc

# export
saveRDS(salmon_tx_se, snakemake@output[["salmon_tx"]])
saveRDS(salmon_se, snakemake@output[["salmon"]])
#saveRDS(terminus_se,snakemake@output[["terminus"]])

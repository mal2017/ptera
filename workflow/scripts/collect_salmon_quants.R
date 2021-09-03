library(tximport)

# tx2gene_fl <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2symbol.tsv"
# tx2feature_fl <- "results/quantification/vanilla_salmon_tes_transcripts/terminus.tx2group.tsv"
# samples_fl <- "config/sample_table.csv"

tx2gene_fl <- snakemake@input[["tx2gene"]]
tx2feature_fl <- snakemake@input[["tx2feature"]]
samples_fl <- snakemake@input[["samples"]]

counts_from_abundance_salmon <- snakemake@params[["counts_from_abundance_salmon"]]
counts_from_abundance_salmon_txout <- snakemake@params[["counts_from_abundance_salmon_txout"]]
counts_from_abundance_terminus <- snakemake@params[["counts_from_abundance_terminus"]]
counts_from_abundance_terminus_txout <- snakemake@params[["counts_from_abundance_terminus_txout"]]

# counts_from_abundance_salmon <- "lengthScaledTPM"
# counts_from_abundance_salmon_txout <- "dtuScaledTPM"
# counts_from_abundance_terminus <- "lengthScaledTPM"
# counts_from_abundance_terminus_txout <- "dtuScaledTPM"

# get the conversion for transcripts
tx2gene <- read.table(tx2gene_fl,header = T)

# get the conversion for the renamed terminus grps w/ symbols in the name
tx2feature <- read.table(tx2feature_fl,header = T)
tx2groupid_symbols <- tx2feature[,c("TXNAME","GROUPID_SYMBOLS")]

# get sample table
samples <- read.csv(samples_fl,header = T)

# get the list of files for quantification
files <- snakemake@input[["salmon_files"]]
# terminus_files <- snakemake@input[["terminus_files"]]
# files <- Sys.glob("results/quantification/vanilla_salmon_tes_transcripts/quant/*/quant.sf")
# terminus_files <- Sys.glob("results/quantification/vanilla_salmon_tes_transcripts/terminus/*/quant.sf")

names(files) <- gsub("\\/+quant.sf","",x=gsub(".+quant\\/+","",files))
# names(terminus_files) <- gsub("\\/+quant.sf","",x=gsub(".+terminus\\/+","",terminus_files))

# reorder to match sample table
# samples <- samples[samples$sample_name %in% names(files),]
files <- files[samples$sample_name]
# terminus_files <- terminus_files[samples$sample_name]

# -----------------------------------------------------------------------------
# Importing expression estimates
# -----------------------------------------------------------------------------

# now summarize using 4 different approaches
# this is using straight salmon results, but not summarizing to gene level. potentially
# useful for differential transcript analysis.
txi_tx <- tximport(files, 
                   type="salmon", 
                   tx2gene=tx2gene, 
                   txOut = T,
                   countsFromAbundance = counts_from_abundance_salmon_txout)

# this is using the salmon files again for transcript level quant, but scaling by
# the groups from terminus
txi_terminus_tx <- tximport(files, 
                   type="salmon", 
                   tx2gene=tx2groupid_symbols, 
                   txOut = T,
                   countsFromAbundance = counts_from_abundance_terminus_txout)

# this is standard salmon - summarize to genes
txi <- tximport(files, 
                  type="salmon", 
                  tx2gene=tx2gene, 
                  txOut = F,
                  countsFromAbundance = counts_from_abundance_salmon)

# this is terminus, done by reading in files from salmon and 
# summarizing to the group level w/ a tx2gene file.
txi_terminus <- tximport(files, 
                type="salmon", 
                tx2gene=tx2groupid_symbols, 
                txOut = F,
                countsFromAbundance = counts_from_abundance_terminus)

# this is terminus, done by reading in files directly from terminus and not collapsing
# txi_straight_terminus <- tximport(terminus_files, 
#                          type="salmon", 
#                          txOut =T,
#                          countsFromAbundance = counts_from_abundance_terminus)
# 
# # leave commented - this is just an example of how close the two terminus approaches
# # (via terminus directly or via salmon are)
# # Mostly this shows me that it is best to just go with pulling them in from salmon, as I can seamlessly
# # rename the the groups to be more informative
# head(txi_terminus$counts["FBtr0006151",])
# head(txi_straight_terminus$counts["FBtr0006151",])

# -----------------------------------------------------------------------------
# pull relevant data from txi objects
# -----------------------------------------------------------------------------

salmon <- txi$counts
salmon_tx <- txi_tx$counts
terminus <- txi_terminus$counts
terminus_tx <- txi_terminus_tx$counts

salmon_gz <- gzfile(snakemake@output[["salmon"]],'w')
write.table(salmon, quote=F, sep = "\t")
close(salmon_gz)

salmon_tx_gz <- gzfile(snakemake@output[["salmon_tx"]],'w')
write.table(salmon_tx, quote=F, sep = "\t")
close(salmon_tx_gz)

terminus_gz <- gzfile(snakemake@output[["terminus"]],'w')
write.table(terminus, quote=F, sep = "\t")
close(terminus_gz)

terminus_tx_gz <- gzfile(snakemake@output[["terminus_tx"]],'w')
write.table(terminus_tx, quote=F, sep = "\t")
close(terminus_tx_gz)


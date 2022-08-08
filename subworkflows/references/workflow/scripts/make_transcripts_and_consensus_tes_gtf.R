library(rtracklayer)

tefafile <- snakemake@input[["te_fasta"]]
gal4fafile <- snakemake@input[["gal4_fasta"]]
# tefafile <- "../../resources/Tidalbase_transposon_sequence.fasta.gz" # tefafile <- "resources/dmel_repbase_lib.fasta.gz"
gtffile <- snakemake@input[["host_gtf"]]
# gtffile <- "resources/dmel-all-r6.41.gtf"

# get te seqs
te_fa <-  c(import(tefafile),import(gal4fafile))


# the widths of the te seqs with be the coord of each te transcript in the ref
end <- width(te_fa)

# extract actual names - should work with either TIDAL or repbase fasta headers as of 210902
raw_te_names <- names(te_fa)

# wtf is univode or escape character between TIDAL TE names after the 3rd vert bar????
# x0 should correspond to the name salmon gives the te in the output file
x0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = raw_te_names, perl = T))
x <- gsub(".+\\|","", x=x0)
x1 <- gsub("#.+","",x = x)

clean_te_names <- gsub("gb\\|.+\\|","",x = x1)

names(end) <- clean_te_names

# create dataframe with te txs
te_gtf_df <- as.data.frame(end, row.names = names(end))

# populate other aspects of the te gtf
te_gtf_df[,"start"] <- 1
te_gtf_df[,"seqnames"] <- raw_te_names
te_gtf_df[,"strand"] <- "+"
te_gtf_df[,"transcript_id"] <- x0
te_gtf_df[,"transcript_symbol"] <- x0
te_gtf_df[,"gene_id"] <- x1
te_gtf_df[,"gene_symbol"] <- x1
te_gtf_df[,"type"] <- "mRNA"
te_gtf_df[,"source"] <- tefafile

# make into a gr
te_gtf <- GRanges(te_gtf_df)

# get host gtf
gtf <- import(gtffile)

# make final combined gr
combined <- c(gtf, te_gtf)

# write to disk
export(combined,snakemake@output[["gtf"]])

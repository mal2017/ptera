library(rtracklayer)

gtffile <-  snakemake@input[["gtf"]] # gtffile <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf"

gtf <- import(gtffile)

# salmon only takes mRNA 
# UPDATE 211213: Decided to quantify all and exclude unwated seqs better, rather than trying to pre-align things
#allowed_types <- c("mRNA")
#gtf2 <- gtf[gtf$type %in% allowed_types]
gtf2 <- gtf

df <- as.data.frame(gtf2)[,c("transcript_id","gene_id","gene_symbol","transcript_symbol")]

# remove straight gene entries (NA in tx id)
df <- df[!is.na(df$transcript_id),]

df <- unique(df)

# select 2 types of conversion
tx2id <- df[,c("transcript_id","gene_id")]
tx2symbol <- df[,c("transcript_id","gene_symbol")]
tx2txsymbol <- df[,c("transcript_id","transcript_symbol")]

# per tximport convention
colnames(tx2id) <- c("TXNAME","GENEID")
colnames(tx2symbol) <- c("TXNAME","GENEID")
colnames(tx2txsymbol) <- c("TXNAME","TXSYMBOL")

write.table(tx2id,file = snakemake@output[["tx2id"]],quote = F, row.names = F,col.names = T, sep = "\t")
write.table(tx2symbol,file = snakemake@output[["tx2symbol"]],quote = F, row.names = F,col.names = T, sep="\t")
write.table(tx2txsymbol,file = snakemake@output[["tx2txsymbol"]],quote = F, row.names = F,col.names = T, sep="\t")

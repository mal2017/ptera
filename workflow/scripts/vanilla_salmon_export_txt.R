library(DESeq2)

set.seed(1)

#dds_fl <- "results/quantification/vanilla_salmon_tes_transcripts/salmon_dds.rds"
dds_fl <- snakemake@input[["dds"]]

dds <- readRDS(dds_fl)

ot <- "raw" # "normcts" "vst" "rlog" "fpkm" "abundance"
ot <- snakemake@wildcards[["expression_unit"]] # "normcts" "vst" "rlog" "fpkm" "abundance"

if (ot == "raw") {
  dat <- counts(dds,normalized=F)
} else if (ot == "normcts") {
  dat <- counts(dds,normalized=T)
} else if (ot == "vst") {
  dat <- assay(vst(dds,blind=T,fitType = snakemake@params[["DESEQ2_VST_FITTYPE"]]))
} else if (ot == "rlog") {
  dat <- assay(rlogTransformation(dds,blind=T,fitType = snakemake@params[["DESEQ2_RLOG_FITTYPE"]]))
} else if (ot == "fpkm") {
  dat <- fpkm(dds)
} else if (ot == "abundance") {
  dat <- assay(dds,"abundance")
} else if (ot == "cpm") {
  dat <- fpm(dds)
}

feature <- rownames(dat)

dat <- cbind(as.data.frame(feature),dat)

#saveRDS(dds,snakemake@output[["dds"]])

o_gz <- gzfile(snakemake@output[["txt"]],'w')
write.table(dat,file = o_gz, quote=F, sep = "\t",row.names = F, col.names = T)
close(o_gz)

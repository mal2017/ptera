library(SummarizedExperiment)

#sex <- "both"
sex <- snakemake@wildcards[["sex"]]
sex <- unlist(ifelse(sex == "both",list(c("male","female")),sex))

#se_fl <- "results/quantification/vanilla_salmon_tes_transcripts/se.gene.rds"
se_fl <- snakemake@input[["se"]]
se <- readRDS(se_fl)
se <- se[,se$sex %in% sex]

#ot <- "ppc_qn" # "normcts" "vst" "rlog" "fpkm" "abundance"
ot <- snakemake@wildcards[["expression_unit"]] # "normcts" "vst" "rlog" "fpkm" "abundance"

stopifnot(ot %in% c("counts","normcts","vst","fpkm","abundance","ppc_qn","edaseq_uq","edaseq_qn"))

# normalization options accessible via DESeq2
if (ot %in% c("normcts","vst","rlog","fpkm","cpm")) {
  library(DESeq2)
  dds <- DESeqDataSet(se, ~ 1)

  if (ot == "normcts") {
    dat <- counts(dds,normalized=T)
  } else if (ot == "vst") {
    dat <- assay(vst(dds,blind=T,fitType = snakemake@params[["DESEQ2_FITTYPE"]]))
  } else if (ot == "fpkm") {
    dat <- fpkm(dds)
  } else if (ot == "cpm") {
    dat <- fpm(dds)
  }

}

# other normalization options
if (ot %in% c("counts","abundance","ppc_qn","edaseq_uq","edaseq_qn")) {
  message(paste("extracting data"))
  dat <- assay(se,ifelse(ot=="abundance","abundance","counts"))

  if (ot == "ppc_qn") {
    message(paste("performing",ot))
    dat <- preprocessCore::normalize.quantiles(dat,copy = T)
    dimnames(dat) <- dimnames(dat)
  } else if (ot == "edaseq_uq") {
    message(paste("performing"),ot)
    dat <- EDASeq::betweenLaneNormalization(dat, which="upper", round = F, offset=F)
  } else if (ot == "edaseq_qn") {
    message(paste("performing"),ot)
    dat <- EDASeq::betweenLaneNormalization(dat, which="full", round = F, offset=F)
  }
}

feature <- rownames(dat)

dat <- cbind(as.data.frame(feature),dat)

#saveRDS(dds,snakemake@output[["dds"]])

o_gz <- gzfile(snakemake@output[["txt"]],'w')
write.table(dat,file = o_gz, quote=F, sep = "\t",row.names = F, col.names = T)
close(o_gz)

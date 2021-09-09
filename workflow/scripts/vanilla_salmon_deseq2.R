library(DESeq2)
library(BiocParallel)

#cores <- 4
cores <- snakemake@threads

register(SnowParam(cores))

#se_fl <- "results/quantification/vanilla_salmon_tes_transcripts/salmon_se.rds"
se_fl <- snakemake@input[["se"]]

se <- readRDS(se_fl)

#form <- "~ Strain + sex"
#form <- snakemake@params[["formula"]]

#reduced <- "~ 1"
#reduced <- snakemake@params[["reduced"]]

dds <- DESeqDataSet(se, ~ 1)

#for (i in 1:100) {
#  dds@assays@data[[paste0("infRep",i)]] <- NULL
#}

#keep <- apply(counts(dds), 1, function(x){sum(x>10) > 8})
#dds <- dds[keep,]

#dds <- dds[1:1000,]

dds <- estimateSizeFactors(dds)

#dds <- DESeq(dds, test="LRT",reduced= as.formula(reduced), parallel = T,fitType = "local")

saveRDS(dds,snakemake@output[["dds"]])

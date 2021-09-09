library(DESeq2)
library(BiocParallel)

#cores <- 4
#cores <- snakemake@threads

#register(SnowParam(cores))

#se_fl <- "results/quantification/vanilla_salmon_tes_transcripts/salmon_se.rds"
se_fl <- snakemake@input[["se"]]

se <- readRDS(se_fl)

#se <- se[1:100,1:6]

#form <- "~ Strain + sex"
form <- snakemake@params[["formula"]]

#reduced <- "~ 1"
reduced <- snakemake@params[["reduced"]]

dds <- DESeqDataSet(se, as.formula(form))

dds <- DESeq(dds, test="LRT",reduced= as.formula(reduced), parallel = F)

saveRDS(dds,snakemake@output[["dds"]])

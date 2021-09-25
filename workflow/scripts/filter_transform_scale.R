library(tidyverse)
library(SummarizedExperiment)
library(rlang)

#sefile <- "results/quantification/vanilla_salmon_tes_transcripts/se.gene.rds"
sefile <- snakemake@input[["sefile"]]
se <- read_rds(sefile)

x <- assay(se,"counts")

# set up filter
#filt <- parse_expr( "(((!str_detect(rownames(x),'FBgn')) & rowSums(x > 1) > 10) | (rowSums(x > 10) > 100)) & rowSums(x == 0) < 0.3*ncol(x)" )
filt <- parse_expr( snakemake@params[["filt"]] )
features_2_use <-  rownames(x[eval(filt),])

rm(x); #rm(se)

#matfile <- "results/quantification/vanilla_salmon_tes_transcripts/male.gene.vst.tsv.gz"
matfile <- snakemake@input[["mat"]]
x <- vroom::vroom(matfile,num_threads = 1)

# ------------------------- perform filtering -------------------------
x <- column_to_rownames(x, "feature")[features_2_use,]


# --------------------- perform transforms ---------------------------
#transforms <- "log2(x+1)"
transforms <- snakemake@params[["transforms"]]
transforms <- parse_expr(transforms)

x <- eval(transforms)

# ------------------------- scaling ----------------------------------

# scale <- T
scale <- snakemake@params[["scale"]]
message(class(x))
message(scale)
if (scale) {
  x <- t(scale(t(x)))
}

# ------------------------- pcs -------------------------------------
#pcs <- 4
pcs <- snakemake@params[["pcs"]]

if (pcs > 0) {
  x <- WGCNA::removePrincipalComponents(x, pcs)
}

# ------------------------ export ----------------------------------

x <- as_tibble(x, rownames = "feature")

x <- arrange(x, -str_detect(feature,"FBgn"))
message(class(x))

write_tsv(x,file = snakemake@output[["mat"]])
#vroom::vroom(x,snakemake@output[["mat"]])

library(rtracklayer)
library(GenomicFeatures)
library(eisaR)
library(tidyverse)

#gtf <- import("../references/results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf")
gtf <- import(snakemake@input[["gtf"]])

tes <- gtf[!grepl("FBgn",gtf$gene_id)]

tes$type <- "exon"

gtf2 <- sort(c(gtf,tes))

export(gtf2,snakemake@output[["gtf"]])


# grl <- suppressWarnings(getFeatureRanges(
#   gtf = file.path("../references/results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.gtf"),
#   featureType = c("spliced", "intron"),
#   intronType = "separate",
#   flankLength = 5,
#   joinOverlappingIntrons = TRUE,
#   verbose = TRUE
# ))

#spliced_grl = grl[str_detect(names(grl), "-", negate = TRUE) | str_detect(names(grl), "FBtr", negate = TRUE) ]


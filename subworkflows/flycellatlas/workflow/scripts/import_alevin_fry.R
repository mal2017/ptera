library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(zellkonverter)

source("workflow/scripts/import_alevin_fry_func.R")

#proj <- "FlyCellAtlas"
proj <- snakemake@wildcards[["sample"]]

#frydir <- "results/alevin-fry/quant_res/FCA14_Female_head_adult_5dWT_Luo_sample5_S5/"
frydir <- paste0(snakemake@input[["frydir"]])

x <- load_fry(frydir = frydir, which_counts = c("S","A")) # load spliced, which included ambig
u <- load_fry(frydir =frydir, which_counts = "U") # load unspliced

assay(x,"unspliced") <- counts(u) # combined into 1 sce

#sds <- CreateSeuratObject(counts = counts(x),assay = "spliced",project = proj) # make initial seur obj
#sds[["unspliced"]] <- CreateAssayObject(assay(x,"unspliced")) # add unspliced

write_rds(x,snakemake@output[["sce"]])
#write_rds(sds,snakemake@output[["seur"]])
#zellkonverter::writeH5AD(sce = x,file = snakemake@output[["h5ad"]])

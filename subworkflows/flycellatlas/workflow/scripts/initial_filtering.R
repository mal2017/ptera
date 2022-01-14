library(tidyverse)
library(miQC)
library(scDblFinder)
library(SingleCellExperiment)
library(scater)
library(flexmix)
library(splines)
library(miQC)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(GenomicFeatures)
library(zellkonverter)

set.seed(2)

#sce_fl <- "results/alevin-fry/sce_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.sce.rds"
#seur_fl <- "results/alevin-fry/seur_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.seur.rds"

sce_fl <- snakemake@input[["sce"]]
seur_fl <- snakemake@input[["seur"]]

sce <- read_rds(sce_fl)
seur <- read_rds(seur_fl)

# remove cells with 0 counts
sce <- sce[,colSums(counts(sce)) > 0]

# get logcounts
sce <- logNormCounts(sce,assay.type=1)

# find doublets
set.seed(2)
scdbf <- scDblFinder(sce)

# get doublet finder data as tbl
doublet_dat <- colData(scdbf) %>% as_tibble(rownames = "cell")

# perform filt
sce <- sce[,scdbf$scDblFinder.class == "singlet"]

# get mito genes
genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
mito_genes <- genes[as.character(seqnames(genes)) == "chrM",]

# only retain those detected here
mito_genes <- mito_genes[mito_genes$gene_id %in% rownames(sce)]

sce <- addPerCellQC(sce, subsets = list(mito=mito_genes$gene_id))

model <- mixtureModel(sce)

# plot filtering scheme. see miqc for other plotting options.
# This encomposses most of the relevant info so is fine for now.
g_miqc <- plotFiltering(sce,model)

sce <- filterCells(sce, model = model)

# also filter the seurat obj.
seur <- 

write_rds(sce,snakemake@output[["sce"]])
write_rds(seur,snakemake@output[["seur"]])
zellkonverter::writeH5AD(sce = sce,file = snakemake@output[["h5ad"]])

# https://www.bioconductor.org/packages/release/bioc/vignettes/velociraptor/inst/doc/velociraptor.html
#library(scran)
#dec <- modelGeneVar(sce)
#top.hvgs <- getTopHVGs(dec,n=2000)

#library(velociraptor)

#velo.out <- scvelo(sce, subset.row=top.hvgs, assay.X="spliced")
#library(scRNAseq)
#sceX <- HermannSpermatogenesisData()
#sceX <- logNormCounts(sceX,assay.type=1)
#decX <- modelGeneVar(sceX)
#top.hvgsX <- getTopHVGs(decX, n=2000)
#velo.out <- scvelo(sceX, subset.row=top.hvgsX, assay.X="spliced",)
# 

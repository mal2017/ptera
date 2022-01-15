library(tidyverse)
library(miQC)
library(scDblFinder)
library(SingleCellExperiment)
library(scater)
library(scran)
library(flexmix)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(GenomicFeatures)
library(celda)

source("../../workflow/scripts/ggplot_theme.R")

set.seed(2)

#sce_fl <- "results/alevin-fry/r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.sce.rds"

sce_fl <- snakemake@input[["sce"]]

sce <- read_rds(sce_fl)

# --------------------
# light clustering reduces incidence of error messages later
# --------------------

# remove cells with 0 counts
sce <- sce[,colSums(counts(sce)) > 0]

# remove unexpressed genes
sce <- sce[rowMeans(counts(sce)) > 0,]

# --------------------
# remove ambient contamination
# --------------------

# run dcontX
set.seed(2)
dcx.sce <- decontX(sce, seed=12345)

counts(sce) <- round(decontXcounts(dcx.sce))

# remove cells/genes with low counts after decontx
sce <- sce[,colSums(counts(sce)) > 1]
sce <- sce[rowSums(counts(sce)) > 3 & rowMeans(counts(sce)) > 0,]

# --------------------
# doublet removal
# --------------------

# get size factors by deconv method - see http://bioconductor.org/books/3.14/OSCA.basic/normalization.html#normalization-transformation
set.seed(2)
clust.sce <- scran::quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clust.sce, min.mean=0.1)

# get logcounts
sce <- logNormCounts(sce,assay.type=1)

# find doublets
set.seed(2)
scdbf <- scDblFinder(sce)

# get doublet finder data as tbl
dat_doublet <- colData(scdbf) %>% as_tibble(rownames = "cell")

g_doublet <- dat_doublet %>%
  arrange(scDblFinder.score) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(rank,scDblFinder.score,color=scDblFinder.class)) +
  geom_point() +
  scale_x_log10() +
  scale_color_grey()

# perform filt
sce <- sce[,scdbf$scDblFinder.class == "singlet"]

# --------------------
# get qc metric
# --------------------

# get mito genes
genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
mito_genes <- genes[as.character(seqnames(genes)) == "chrM",]

#only retain those actually quantified here
mito_genes <- mito_genes[mito_genes$gene_id %in% rownames(sce)]

sce <- addPerCellQCMetrics(sce, subsets = list(mito=mito_genes$gene_id))

# --------------------
# remove outliers + low info cells
# --------------------

# next bit is low total outlier filtering, via scuttle
keep.umis <- !isOutlier(sce$sum,type = "both",log = T,nmads = 3)
keep.genes <- !isOutlier(sce$detected,type = "both",log = T, nmads = 3)

sce <- sce[,keep.umis & keep.genes]

sce <- sce[,colSums(counts(sce)) > 500 & sce$detected > 200]

# ---------------------------
#  remove mito reads
# --------------------------

# make model for mito removal
model <- mixtureModel(sce)

# plot filtering scheme. see miqc for other plotting options.
# This encomposses most of the relevant info so is fine for now.

# sometimes the above works fine
no_mito_outliers <- length(unique(model@cluster)) == 1

g_miqc0 <- plotMetrics(sce,model)

dat_miqc <- g_miqc0$data %>% as_tibble(rownames = "feature")

g_miqc <- ggplot(dat_miqc,aes(sum,subsets_mito_percent)) +
  geom_point()

if (!no_mito_outliers) {
  sce <- filterCells(sce, model = model)  
}


# --------------
# more qc metrics
# -----------------

# add some extra qc info before shipping out
sce <- addPerFeatureQCMetrics(sce)

# --------------
# last bits and export
# -----------------

# also filter the seurat obj.
# seurat replaces underscores with dashes
#seur <- seur[str_replace_all(rownames(sce),pattern = "_",replacement = "-"),colnames(sce)]

write_rds(sce,snakemake@output[["sce"]])
#write_rds(seur,snakemake@output[["seur"]])

write_tsv(dat_doublet,snakemake@output[["dat_doublet"]])
saveRDS(g_doublet,snakemake@output[["g_doublet"]])
ggsave(snakemake@output[["png_doublet"]],g_doublet, width=4, height = 4)

write_tsv(dat_miqc,snakemake@output[["dat_miqc"]])
saveRDS(g_miqc,snakemake@output[["g_miqc"]])
ggsave(snakemake@output[["png_miqc"]],g_miqc, width=4, height = 4)

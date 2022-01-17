library(scater)
library(scran)
library(bluster)
library(tidyverse)

source("../../workflow/scripts/ggplot_theme.R")

# see http://bioconductor.org/books/3.14/OSCA.advanced/

#sce_fl <- "results/downstream/single_sample_dimred/FCA7_Female_Adult_Head1_S4.usa.filt.dimred.sce.rds"

sce_fl <- snakemake@input[["sce"]]

sce <- read_rds(sce_fl)

set.seed(2)
sce <- runTSNE(sce,dimred="PCA.elbow")
set.seed(2)
sce <- runUMAP(sce,dimred="PCA.elbow")

# ------------------------
# Clustering options
# ------------------------

#clusts <- clusterCells(sce, use.dimred = "PCA.elbow")

clusts <- clusterCells(sce, use.dimred="PCA.elbow", BLUSPARAM=NNGraphParam(cluster.fun = "walktrap",type="rank",k=15))

colLabels(sce) <- clusts

# -----------------
# Get sil data
# -----------------

sil.approx <- approxSilhouette(reducedDim(sce, "PCA.elbow"), clusters=colLabels(sce))

sil.data <- as.data.frame(sil.approx)

# -----------------
# smooth the cluster assignments
# -----------------

sil.data$closest <- factor(ifelse(sil.data$width > 0, colLabels(sce), sil.data$other))

sil.data$cluster <- colLabels(sce)

stopifnot(all(colnames(sce) == rownames(sil.data)))

sce$label <- sil.data$closest

# -----------------
# plots
# -----------------

dat_sil <- sil.data %>% as_tibble()

set.seed(2)
g_sil <- ggplot(dat_sil,aes(cluster,width)) +
  geom_jitter(size=rel(0.7))

g_tsne <- plotTSNE(sce,colour_by="label")
g_umap <- plotUMAP(sce,colour_by="label")

dat_tsne <- g_tsne$data
dat_umap <- g_umap$data

# -----------------
# export
# -----------------

write_rds(sce,snakemake@output[["sce"]])

write_tsv(dat_sil,snakemake@output[["dat_sil"]])
saveRDS(g_sil,snakemake@output[["g_sil"]])
ggsave(snakemake@output[["png_sil"]],g_sil, width=4, height = 4)

write_tsv(dat_tsne,snakemake@output[["dat_tsne"]])
saveRDS(g_tsne,snakemake@output[["g_tsne"]])
ggsave(snakemake@output[["png_tsne"]],g_tsne, width=4, height = 4)

write_tsv(dat_umap,snakemake@output[["dat_umap"]])
saveRDS(g_umap,snakemake@output[["g_umap"]])
ggsave(snakemake@output[["png_umap"]],g_umap, width=4, height = 4)


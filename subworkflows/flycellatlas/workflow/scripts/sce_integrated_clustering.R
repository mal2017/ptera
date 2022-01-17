library(scater)
library(scran)
library(tidyverse)
library(bluster)


source("../../workflow/scripts/ggplot_theme.R")



# ---------------------------
# get clusts
# ---------------------------
#sce_fl <- "results/downstream/integrated_by_tissue_and_sex/head_female.sce.rds"

sce_fl <- snakemake@input[["sce"]]

sce <- readRDS(sce_fl)

# just for testing
#sce <- sce[, sample(seq_len(ncol(sce)),size = 5e3,replace = F)]

clusts <- clusterCells(sce, use.dimred="corrected", BLUSPARAM=NNGraphParam(shared = F,k=5))

colLabels(sce) <- clusts

# -----------------
# Get sil data
# -----------------

sil.approx <- approxSilhouette(reducedDim(sce, "corrected")[colnames(sce),], clusters=colLabels(sce))

sil.data <- as.data.frame(sil.approx)[colnames(sce),]

# -----------------
# smooth the cluster assignments
# -----------------

sil.data$closest <- factor(ifelse(sil.data$width > 0, colLabels(sce), sil.data$other))

sil.data$cluster <- colLabels(sce)

stopifnot(all(colnames(sce) == rownames(sil.data)))

sce$label <- sil.data$closest


# ---------------------------
# plot em
# ---------------------------

dat_sil <- sil.data %>% as_tibble()

set.seed(2)
g_sil <- ggplot(dat_sil,aes(cluster,width)) +
  geom_jitter(size=rel(0.7))

dat_tsne <- reducedDim(sce,"TSNE") %>% magrittr::set_colnames(c("TSNE1","TSNE2")) %>%
  as_tibble(rownames = "cell") %>%
  left_join(as_tibble(colData(sce),rownames="cell"))
  #pivot_longer(-c("cell","TSNE1","TSNE2"),names_to = "color_type",values_to = "color_by") %>%
  
g_tsne <- ggplot(dat_tsne, aes(TSNE1,TSNE2,color=label)) +
  geom_point() +
  coord_fixed()

# ---------------------------
# export
# ---------------------------
write_rds(sce,snakemake@output[["sce"]])

write_tsv(dat_sil,snakemake@output[["dat_sil"]])
saveRDS(g_sil,snakemake@output[["g_sil"]])
ggsave(snakemake@output[["png_sil"]],g_sil, width=8, height = 8)

write_tsv(dat_tsne,snakemake@output[["dat_tsne"]])
saveRDS(g_tsne,snakemake@output[["g_tsne"]])
ggsave(snakemake@output[["png_tsne"]],g_tsne, width=8, height = 8)

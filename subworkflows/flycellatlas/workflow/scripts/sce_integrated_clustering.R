library(scater)
library(scran)
library(tidyverse)


source("../../workflow/scripts/ggplot_theme.R")

sce <- readRDS("results/downstream/integrated_by_tissue_and_sex/head_female.sce.rds")



clusts <- clusterCells(sce, use.dimred = "corrected")

colLabels(sce) <- clusts

g_labels <- plotTSNE(sce,colour_by="label")


g

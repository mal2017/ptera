library(scater)
library(scran)
library(tidyverse)
library(monocle3)

source("../../workflow/scripts/ggplot_theme.R")

#sce_fl <- "results/downstream/integrated_by_tissue_and_sex_clustering/head_female.sce.rds"

sce_fl <- snakemake@input[["sce"]]

sce <- read_rds(sce_fl)

# ---------------------------
# convert to cds
# ---------------------------

X <- assay(sce,i = 1) %>% as.matrix()


cds <- new_cell_data_set(X,
                         cell_metadata = colData(sce),
                         gene_metadata = rowData(sce))


# ---------------------------
# convert to cds
# ---------------------------
cds

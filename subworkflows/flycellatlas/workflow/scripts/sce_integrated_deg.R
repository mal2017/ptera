library(scater)
library(scran)
library(tidyverse)
library(bluster)


source("../../workflow/scripts/ggplot_theme.R")


# ---------------------------
# get clusts
# ---------------------------
#merged_fl <- sce_fl <- "results/downstream/integrated_by_tissue_and_sex_clustering/ovary_female.sce.rds"
merged_fl <- snakemake@input[["merged"]]

merged <- readRDS(merged_fl)

#sce_fls <- Sys.glob("results/downstream/single_sample_clustering/*Female*ovary*usa.filt.dimred.clust.sce.rds")
sce_fls <- snakemake@input[["sces"]]

names(sce_fls) <- str_extract(sce_fls,"(?<=clustering\\/).+(?=\\.usa)")

sces <- sce_fls %>% map(read_rds)

# ---------------------------
# workflow for diffs from
# http://bioconductor.org/books/3.14/OSCA.multisample/using-corrected-values.html#after-blocking-on-the-batch
# ---------------------------

universe <- Reduce(intersect, lapply(sces, rownames))
all.sce2 <- lapply(sces, "[", i=universe,)
#all.dec2 <- lapply(decs, "[", i=universe,)

all.sce2 <- lapply(all.sce2, function(x) {
  reducedDims(x) <- reducedDims(x)["PCA"]
  reducedDims(x) <- NULL
  x
})


all.sce2 <- lapply(all.sce2, function(x) {
  rowData(x) <- rowData(all.sce2[[1]])
  x
})

all.sce2 <- lapply(all.sce2, function(x) {
  colData(x) <- colData(x)[,1,drop=F]
  x
})

combined <- do.call(cbind, all.sce2)

# combined$batch <- rep(names(all.sce2), vapply(all.sce2, ncol, 0L))

clusters.mnn <- colLabels(merged)

rm(sces)
rm(all.sce2)
rm(merged)

# ---------------------------------
# see combineMarkers()
#
# The pval.type="some" setting serves as a compromise between "all" and "any".
# A combined p-value is calculated by taking the middlemost value of the Holm-corrected p-values for each gene.
# (By default, this the median for odd numbers of contrasts and one-after-the-median for even numbers,
# but the exact proportion can be changed by setting min.prop - see ?combineParallelPValues.)
# Here, the null hypothesis is that the gene is not DE in at least half of the contrasts.
# Genes are then ranked by the combined p-value.
# The aim is to provide a more focused marker set without being overly stringent,
# though obviously it loses the theoretical guarantees of the more extreme settings.
# For example, there is no guarantee that the top set contains genes that can distinguish a cluster from any other cluster,
# which would have been possible with pval.type="any".
# For each gene and cluster, the summary effect size is defined as
# the effect size from the pairwise comparison with the min.prop-smallest p-value.
# This mirrors the p-value calculation but, again, is reported only for the benefit of the user.

# this means, by my unserstanding, that the effect size is from the smallesst p value
# where min prop is satisfied. For pval.type="some", min prop will be 0.5
# ------------------------------

m.out <- findMarkers(combined, clusters.mnn, block=combined$batch,
                     row.data=rowData(combined)[,c(1,2)],
                     pval.type="any",
                     min.prop=0.25,
                     test.type="t")


diffs_tbl_list <- m.out %>% as.list() %>%
  map(as_tibble,rownames="feature")

diffs_tbl <- diffs_tbl_list %>% map_df(select,feature:summary.logFC,.id = "label")

write_tsv(diffs_tbl, snakemake@output[["tsv"]])

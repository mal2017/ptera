library(tidyverse)
library(rtracklayer)

# ----- process ins to grl ------
#ins_fl <- "results/overlaps/collected_insertions/insertions_by_strain.tsv.gz"
ins_fl <- snakemake@input[["ins"]]

ins <- read_tsv(ins_fl)

#genes_fl <- "resources/dmel-all-r6.41.gtf"
genes_fl <- snakemake@input[["genes"]]

genes <- rtracklayer::import(genes_fl)

genes <- genes[genes$type == "gene"]

genes$name <- genes$gene_id

# ---- --------
get_overlaps_per_strain <- function(ins_obj, genes, feature_df) {
  ins_gr <- GRanges(ins_obj)

  gr <- c(ins_gr, genes)

  feature_df <- as_tibble(gr) %>% mutate(queryHits=row_number())

  # should we really be counting overlaps here?
  findOverlaps(gr,gr,ignore.strand=TRUE, maxgap = snakemake@params[["maxgap"]]) %>%
    as_tibble() %>%
    filter(queryHits != subjectHits) %>%
    right_join(feature_df,., by="queryHits") %>%
    left_join(feature_df, by=c(subjectHits="queryHits")) %>%
    dplyr::select(name.x, gene_symbol.x,name.y, gene_symbol.y) %>%
    filter(name.x!=name.y) %>%
    distinct() %>%
    mutate(gene_symbol.x = ifelse(is.na(gene_symbol.x),name.x, gene_symbol.x)) %>%
    mutate(gene_symbol.y = ifelse(is.na(gene_symbol.y),name.y, gene_symbol.y))
}

res <- ins %>%
  group_by(strain) %>%
  do(get_overlaps_per_strain(.,genes, feature_df))

res <- res %>%
  ungroup() %>%
  dplyr::rename(Strain=strain) %>%
  mutate(overlap=T)


write_tsv(res,snakemake@output[["tsv"]])

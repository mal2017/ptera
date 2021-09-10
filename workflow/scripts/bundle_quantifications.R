library(tidyverse)

#expression <- "results/quantification/vanilla_salmon_tes_transcripts/cpm.tsv.gz"
expression <- snakemake@input[["expression"]]

#sqlite <- "results/quantification/vanilla_salmon_tes_transcripts/cpm.sqlite"
sqlite <- snakemake@output[["sqlite"]]

#metadata <- "results/meta/metadata.csv"
metadata <- snakemake@input[["meta"]]

df <- vroom::vroom(expression)

gene_df <- df %>% filter(grepl("FBgn",feature)) %>%
  pivot_longer(cols = -"feature", names_to = "sample", values_to = "x") %>%
  dplyr::select(sample,host_gene = "feature",x)

transposon_df <- df %>% filter(!grepl("FBgn",feature)) %>%
  pivot_longer(cols = -"feature", names_to = "sample", values_to = "y") %>%
  dplyr::select(sample,transposon = "feature",y)

metadata_df <- vroom::vroom(metadata)  %>%
  dplyr::rename(sample=sample_name)#%>% mutate(wolbachia = as.numeric(wolbachia))

con <- DBI::dbConnect(RSQLite::SQLite(), sqlite)

copy_to(con, metadata_df, name = "sample_info", temporary=F)

copy_to(con, gene_df, name = "genes", temporary=F)

copy_to(con, transposon_df, name = "transposons", temporary=F)

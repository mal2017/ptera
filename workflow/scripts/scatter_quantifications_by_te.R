library(tidyverse)
library(vroom)

#meta_fl <- "results/meta/metadata.csv"
meta_fl <- snakemake@input[["meta"]]
metadata <- read_csv(meta_fl)

#expression_fl <- "~/work/ptera/results/quantification/vanilla_salmon_tes_transcripts/raw.tsv.gz"
expression_fl <- snakemake@input[["expression"]]

mat <- read_tsv(expression_fl)

df <- pivot_longer(mat, cols = -"feature", names_to = "sample", values_to = "x")

#df <- df %>% filter(!str_detect(feature,"FBgn") | feature %in% c("FBgn0085432","FBgn0004872"))

df <- pivot_wider(df, names_from = "feature",values_from = "x") %>%
  pivot_longer(contains("FBgn"),names_to = "host_gene",values_to = "x") %>%
  pivot_longer(cols = -c(sample,host_gene,x),names_to = "transposon", values_to="y")

df <- left_join(metadata,df, by=c("sample_name"="sample")) %>%
  relocate(sample_name,Strain,sex, wolbachia, host_gene,x, transposon, y)

#path <- "~/work/manual_analysis/210908_lm_ptera_v1/matrices"
path <- snakemake@output[[1]]
dir.create(path)

df %>%
  mutate(batch = transposon) %>%
  group_by(batch) %>%
  group_walk(~ vroom::vroom_write(.x, sprintf("%s/%s.tsv.gz", path, .y$batch)))

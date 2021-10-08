library(tidyverse)

feat_f <-  snakemake@input[["feats"]]
#feat_f <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2id.tsv"

feat <- read_tsv(feat_f)

cd_f <- snakemake@input[["samples"]]
#cd_f <- "config/sample_table.csv"

cd <- read_csv(cd_f) %>% dplyr::select(Strain) %>% distinct()

genes <- feat %>% dplyr::select(sequence = GENEID) %>% distinct()
transcripts <- feat %>% dplyr::select(sequence = TXNAME) %>% distinct()

genes <- crossing(cd, genes) %>% mutate(est.copies = 1) %>%
  mutate(sample_name = NA, length = NA, bases=NA, median.cov=NA) %>%
  dplyr::select(Strain, sample_name, sequence, length, bases, median.cov, est.copies)

transcripts <- crossing(cd, transcripts) %>% mutate(est.copies = 1) %>%
  mutate(sample_name = NA, length = NA, bases=NA, median.cov=NA) %>%
  dplyr::select(Strain, sample_name, sequence, length, bases, median.cov, est.copies)

out <- bind_rows(genes, transcripts) %>% distinct()

write_tsv(out, snakemake@output[["feats"]])


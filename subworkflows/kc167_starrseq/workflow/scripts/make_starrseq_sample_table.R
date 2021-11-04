library(tidyverse)

df <- read_csv("../../resources/kc167_starrseq_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(Sample.Name,source_name, cell_type, genotype = "genotype.variation", treatment = "treated_with", Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, Experiment, Run, BioSample) %>%
  mutate(sample_name = paste(source_name,Run, sep="_"))

df2 <- df2 %>%
  dplyr::relocate(sample_name, genotype)

input <- unique(df2 %>% filter(is.na(genotype) | is.na(treatment)) %>% pull(sample_name))

df2 <- mutate(df2,background = input)

df2 %>% dplyr::select(sample_name, genotype, treatment, background) %>%
  distinct() %>% write_csv("config/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

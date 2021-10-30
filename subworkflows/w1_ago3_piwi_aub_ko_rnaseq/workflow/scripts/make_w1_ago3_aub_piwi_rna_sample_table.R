library(tidyverse)

df <- read_csv("../../resources/w1_ago3_aub_piwi_ko_rnaseq_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(Sample.Name,Genotype, Strain, Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, Experiment, Run, BioSample) %>%
  mutate(Sample.Name = str_replace_all(Sample.Name,"\\/|;","_")) %>%
  mutate(Sample.Name = str_remove_all(Sample.Name,"[[:space:]]")) %>%
  mutate(sample_name = paste(Sample.Name,Run, sep="_"))

df2 <- df2 %>%
  dplyr::relocate(sample_name, Genotype)

df2 %>% dplyr::select(sample_name, Genotype) %>%
  distinct() %>% write_csv("config/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

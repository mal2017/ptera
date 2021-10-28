library(tidyverse)

df <- read_csv("../../resources/dpan_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument,  Strain=strain, Experiment, Run, BioSample, Sample.Name) %>%
  distinct() %>%
  mutate(Strain = 'Kc167') %>%
  mutate(sample_name = Sample.Name)

df2 <- df2 %>%
  mutate(genotype = str_extract(sample_name,"^.+(?=-CM|-WCM)")) %>%
  mutate(treatment = str_extract(sample_name,"WCM|CM")) %>%
  dplyr::relocate(sample_name,genotype,treatment)

df2 %>% dplyr::select(sample_name, genotype, treatment) %>%
  distinct() %>% write_csv("config/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

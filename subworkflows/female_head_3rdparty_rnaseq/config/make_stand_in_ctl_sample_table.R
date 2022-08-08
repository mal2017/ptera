library(tidyverse)

df <- read_csv("SraRunTable.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, Genotype, sex, tissue_type, Experiment, Run, BioSample, Sample.Name) %>%
  distinct() %>%
  mutate(sample_name = Sample.Name)

df3 <- df2 %>% filter(tissue_type == "Head" & sex == "female") %>%
  dplyr::relocate(sample_name,Genotype,sex, tissue_type)


df3 %>% dplyr::select(sample_name,Genotype,sex, tissue_type) %>%
  distinct() %>% write_csv("sample_table.csv")

df3 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("subsample_table.csv")

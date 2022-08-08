library(tidyverse)

df <- read_csv("SraRunTable.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument,  Strain=strain, Genotype, Age, cultivation_temperature, Experiment, Run, BioSample, Sample.Name) %>%
  distinct() %>%
  mutate(sample_name = Sample.Name)

df2 <- df2 %>%
  dplyr::relocate(sample_name,Genotype,Age, cultivation_temperature)

df2 %>% dplyr::select(sample_name, Genotype,Age, cultivation_temperature) %>%
  distinct() %>% write_csv("sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("subsample_table.csv")

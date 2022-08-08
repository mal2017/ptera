library(tidyverse)

df <- read_csv("data/SraRunTable.txt") %>%
  set_tidy_names(syntactic = T)

df1 <-  df %>% filter(Assay.Type == "RNA-Seq")

df2 <- df1 %>%
  dplyr::select(Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, Genotype, Experiment, Run, BioSample, Sample.Name) %>%
  distinct() %>%
  mutate(sample_name = Sample.Name)

df2 <- df2 %>%
  dplyr::relocate(sample_name,Genotype)

df2 %>% dplyr::select(sample_name, Genotype) %>%
  distinct() %>% write_csv("data/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("data/subsample_table.csv")

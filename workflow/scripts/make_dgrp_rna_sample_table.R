library(tidyverse)

df <- read_csv("resources/dgrp_expression_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(source_name, Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, sex, age, medium, Strain, Tissue, Experiment, Run, BioSample) %>%
  distinct() %>%
  mutate(sample_name = paste(str_replace_all(source_name," ","_"),BioSample, sep="_"))

df2 %>% dplyr::select(sample_name, source_name, Strain, sex, age) %>%
  distinct() %>% write_csv("config/sample_table.csv")


df2 %>% dplyr::select(-source_name, Strain) %>% distinct() %>%
  dplyr::relocate(sample_name, BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

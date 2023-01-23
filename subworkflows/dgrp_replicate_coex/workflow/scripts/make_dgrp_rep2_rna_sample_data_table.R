library(tidyverse)

# PRJNA773066

df <- read_csv("dgrp_rep2_expression_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df <- df %>% 
  mutate(Tissue= source_name) %>%
  mutate(sample_name = paste("DGRP_Line",LINE,stringr::str_to_title(sex),BioSample,sep="_")) %>%
  mutate(source_name = str_replace_all(str_remove(sample_name,"_SAMN.+"),"_"," ")) %>%
  mutate(Strain = paste0("DGRP_",LINE)) %>%
  mutate(age = "adult") %>%
  mutate(medium=exposure)

df2 <- df %>%
  dplyr::select(sample_name, source_name, Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, sex, age, Strain, medium, Tissue, , Experiment, Run, BioSample) %>%
  distinct()

df2 %>% dplyr::select(sample_name, source_name, Strain, sex, age) %>%
  distinct() %>% write_csv("sample_table.csv")

df2 %>% dplyr::select(-source_name, Strain) %>% distinct() %>%
  dplyr::relocate(sample_name, BioSample) %>%
  distinct() %>%
  write_csv("subsample_table.csv")

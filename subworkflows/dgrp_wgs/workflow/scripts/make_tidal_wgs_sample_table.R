library(tidyverse)

df <- read_csv("../../resources/tidal_wgs_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument,  Strain=strain, Experiment, Run, BioSample) %>%
  distinct() %>%
  mutate(Strain = str_replace_all(Strain,"-","_")) %>%
  mutate(sample_name = paste(Strain,BioSample, sep="_"))

df2 <- df2 %>% filter(LibraryLayout == "PAIRED") %>%
  distinct()

df2 %>% dplyr::select(sample_name, Strain) %>%
  distinct() %>% write_csv("config/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

library(tidyverse)

df <- read_csv("../../resources/tfrnai_rnaseq_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(rnai=RNAi_target_gene_name, Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Strain = Cell_Line, Instrument, Experiment, Run, BioSample, Sample.Name) %>%
  distinct() %>%
  mutate(sample_name = paste(rnai,Sample.Name, sep="_"))

df2 <- df2 %>%
  dplyr::relocate(sample_name, rnai)

df2 %>% dplyr::select(sample_name, rnai) %>%
  distinct() %>% write_csv("config/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

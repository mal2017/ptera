library(tidyverse)

df <- read_csv("SraRunTable.txt") %>%
  set_tidy_names(syntactic = T)

df1 <-  df %>% filter(Assay.Type == "RNA-Seq")

df2 <- df1 %>%
  dplyr::select(Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, Genotype, Age, Developmental_Stage, fertility_phenotype, Experiment, Run, BioSample, Sample.Name) %>%
  distinct() %>%
  mutate(sample_name = Sample.Name)


df2 <- df2 %>%
  dplyr::relocate(sample_name,Genotype,Age,Developmental_Stage, fertility_phenotype)

df3 <- df2 %>% filter(str_detect(Genotype, "TRiP|VALIUM"))

df3 <- df3 %>% mutate(Genotype = str_remove(Genotype,"Genotype: "))

df3 %>% dplyr::select(sample_name,Genotype,Age,Developmental_Stage, fertility_phenotype) %>%
  distinct() %>% write_csv("sample_table.csv")

df3 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("subsample_table.csv")

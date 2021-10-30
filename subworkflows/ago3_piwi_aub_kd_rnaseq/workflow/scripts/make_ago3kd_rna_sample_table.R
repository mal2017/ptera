library(tidyverse)

df <- read_csv("../../resources/ago3piwiaub_kd_rnaseq_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>%
  dplyr::select(kd=germline_knockdown_and_other_transgenes, Assay.Type,  LibraryLayout, LibrarySource, LibrarySelection, Instrument, Experiment, Run, BioSample, Sample.Name) %>%
  mutate(kd=tidy_names(kd,syntactic = T)) %>%
  mutate(kd = str_remove_all(kd,"sh|sh1|attP2|attp40|gfp")) %>%
  mutate(kd = str_remove_all(kd,"\\.|1|2|control")) %>%
  mutate(kd = ifelse(kd == "aubago3","aub_ago3",kd)) %>%
  distinct() %>%
  mutate(sample_name = paste(kd,Sample.Name, sep="_"))

df2 <- df2 %>%
  dplyr::relocate(sample_name, kd)

df2 %>% dplyr::select(sample_name, kd) %>%
  distinct() %>% write_csv("config/sample_table.csv")

df2 %>% distinct() %>%
  dplyr::relocate(sample_name,BioSample) %>%
  distinct() %>%
  write_csv("config/subsample_table.csv")

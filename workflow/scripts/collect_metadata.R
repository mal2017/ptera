library(tidyverse)
library(readxl)

sample_df <- read_csv("config/sample_table.csv")

wolbachia <- read_xlsx("resources/wolbachia.xlsx") %>%
  dplyr::rename(line = `DGRP Line`, wolbachia = `Infection Status`) %>%
  mutate(Strain = paste0("DGRP_",str_extract(line, "(?<=__).+"))) %>%
  mutate(wolbachia = wolbachia == "y") %>%
  dplyr::select(Strain, wolbachia)

left_join(sample_df, wolbachia) %>%
  dplyr::select(-source_name) %>%
  write_csv(snakemake@output[[1]])

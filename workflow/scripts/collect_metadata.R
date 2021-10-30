library(tidyverse)
library(readxl)

sample_df <- read_csv(snakemake@input[["sample_table"]])

wolbachia <- read_xlsx(snakemake@input[["wolbachia"]]) %>%
  dplyr::rename(line = `DGRP Line`, wolbachia = `Infection Status`) %>%
  mutate(Strain = paste0("DGRP_",str_extract(line, "(?<=__).+"))) %>%
  mutate(wolbachia = wolbachia == "y") %>%
  dplyr::select(Strain, wolbachia)

left_join(sample_df, wolbachia) %>%
  dplyr::select(-source_name) %>%
  write_csv(snakemake@output[[1]])

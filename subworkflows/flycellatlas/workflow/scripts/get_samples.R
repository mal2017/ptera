library(tidyverse)

fl <- read_tsv("../../resources/E-MTAB-10519.sdrf.txt")

x <- fl %>%
  dplyr::select(sample_name = `Source Name`,
                individual = `Characteristics[individual]`,
                `age`= `Characteristics[age]`,
                time_unit = `Unit[time unit]`,
                developmental_stage = `Characteristics[developmental stage]`,
                sex = `Characteristics[sex]`,
                tissue = `Characteristics[organism part]`,
                cell_type = `Characteristics[cell type]`,
                fq1_uri = `Comment[FASTQ_URI]...56`,
                fq2_uri = `Comment[FASTQ_URI]...58`,
                fq3_uri = `Comment[FASTQ_URI]...60`,
                library_type = `Comment[library construction]`,
                cdna_read = `Comment[cdna read]`,
                umi_read = `Comment[umi barcode read]`,
                barcode_read = `Comment[sample barcode read]`)

x2 <- x %>% 
  filter(tissue %in% c("testis","ovary","head"))

sample_table <- x2 %>% dplyr::select(-fq2_uri,-fq1_uri,, -fq3_uri, -library_type,
                     -cdna_read,-umi_read,-barcode_read) %>%
  distinct()

subsample_table <- x2 %>% select(sample_name,fq1_uri,fq2_uri,fq3_uri,
              library_type,cdna_read,umi_read,barcode_read) %>%
  distinct()

write_csv(sample_table,"config/sample_table.csv")
write_csv(subsample_table,"config/subsample_table.csv")

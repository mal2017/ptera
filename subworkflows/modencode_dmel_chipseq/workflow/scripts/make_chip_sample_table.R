library(tidyverse)

df <- read_csv("../../../../resources/modencode_dmel_chipseq_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>% filter(str_detect(source_name,regex("tcf|gro|h3|h4|myc|Dll|dTBP",ignore_case = T))) %>%
  dplyr::select(source_name,Run,Experiment,Sample.Name,LibraryLayout,sex,transgene,Instrument,Run,Sample.Name, strain, tissue, Genotype, cell_line, 
                Library.Name, antibody_name, Development_stage, Developmental_Stage) %>%
  mutate(Developmental_Stage = ifelse(is.na(Developmental_Stage),Development_stage,Developmental_Stage)) %>%
  mutate(sample_name = tidy_names(paste(source_name,Run, sep="_"),syntactic = T)) %>%
  dplyr::relocate(sample_name,Developmental_Stage)

# get only samples we want for now
df2 <- df2 %>% filter(!str_detect(source_name,"BG3|nonISO1")) %>% filter(is.na(transgene))

# add inputs. If adding samples, make sure this approach still yields logical names....
df2 <- df2 %>% mutate(input_source_name = str_replace(source_name,regex("(?<=seq_)ChIP"),"Input"))

# a special case on further inspection
special_kc<-df2 %>%
  filter(str_detect(sample_name,"dMyc") & str_detect(sample_name,"Kc167") & !str_detect(sample_name,"W3L")) %>%
  pull(sample_name)

special_s2<-df2 %>%
  filter(str_detect(sample_name,"dMyc") & str_detect(sample_name,"S2") & !str_detect(sample_name,"W3L")) %>%
  pull(sample_name)

df2 <- df2 %>% dplyr::select(input=sample_name,input_source_name=source_name) %>% 
  left_join(df2,.) %>%
  dplyr::relocate(input,.after = sample_name)

df2 <- df2 %>%
  mutate(input = ifelse(sample_name %in% special_kc,"dMyc.Kc167.late.embryo.Input_SRR1199483",input)) %>%
  mutate(input = ifelse(sample_name %in% special_s2,"dMyc.S2.late.embryo.Input_SRR1199486",input))

stopifnot(all(df2$input %in% df2$sample_name))


 df2 %>% dplyr::select(sample_name, source_name, Developmental_Stage, input) %>%
  distinct() %>% write_csv("../../config/sample_table.csv")

df2 %>%
  distinct() %>%
  write_csv("../../config/subsample_table.csv")

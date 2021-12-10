library(tidyverse)

df <- read_csv("../../../../resources/gro_cellline_paper_srarunselector.txt") %>%
  set_tidy_names(syntactic = T)

df2 <- df %>% dplyr::select(cell_line, LibraryLayout,Library.Name,Genotype,genetic_modification,Assay.Type,Run,Experiment,Instrument) %>%
  mutate(Library = Library.Name) %>%
  mutate(sample_name = tidy_names(paste(Library.Name,Run, sep="_"),syntactic = T)) %>%
  dplyr::relocate(sample_name)

# get only samples we want for now
df3 <- df2 %>% filter(str_detect(Assay.Type,"ChIP"))

# add inputs. If adding samples, make sure this approach still yields logical names....
df4 <- df3 %>% 
  mutate(input_source_name = str_replace(Library.Name,regex("ChIP"),"Input")) %>%
  mutate(input_source_name = str_remove(input_source_name,"-\\d{1}$"))

# special cases on further inspection

df5 <- df4 %>% mutate(input_source_name = ifelse(input_source_name %in% c("Input-H3ac-RNAi","Input-H4ac-RNAi","Input-H4ac","Input-H3ac","Input-polII-RNAi","Input-polII"),"Input-H3ac-H4ac-pol",input_source_name))
df5 <- df5 %>% mutate(input_source_name = ifelse(str_detect(input_source_name,"^Input-Rpb3"),"Input-Rpb3",input_source_name))
df5 <- df5 %>% mutate(input_source_name = ifelse(input_source_name == "Input-gro-RNAi","Input-gro",input_source_name))
df5 <- df5 %>% mutate(input_source_name = ifelse(input_source_name %in% c("Input-groGFP","Input-mutGFP"),"Input-GFP",input_source_name))
df5 <- df5 %>% mutate(input_source_name = ifelse(input_source_name %in% c("S2-Input-groGFP"),"S2-Input-gro",input_source_name))

stopifnot(nrow(df5[!df5$input_source_name %in% df4$Library.Name,c("sample_name","input_source_name","Library.Name")])==0)

inputs <- df5 %>% filter(str_detect(Library.Name,"Input")) %>%
  dplyr::rename(input=sample_name, library_layout_input = LibraryLayout)

df6 <- df5 %>% left_join(inputs[c("input","Library.Name","library_layout_input")], by=c(input_source_name="Library.Name"))

stopifnot(all(df6$input %in% df6$sample_name))

df6 %>% dplyr::select(sample_name, cell_line, Library.Name, LibraryLayout, Genotype, genetic_modification, input) %>%
  distinct() %>% write_csv("../../config/sample_table.csv")

df6 %>%
  distinct() %>%
  write_csv("../../config/subsample_table.csv")

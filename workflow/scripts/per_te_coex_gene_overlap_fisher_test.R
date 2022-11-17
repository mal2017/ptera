library(tidyverse)

#current_model <- "male_model_01"
current_model <- snakemake@params[["model_id"]]

#edges_fl <- "subworkflows/dgrp_coex/results/linear_models/male_model_01.gene_2_gene_edges.rds"
edges_fl <- snakemake@input[["edges"]]
#nudges_fl <- "subworkflows/dgrp_coex/results/linear_models/male_model_01.gene_2_gene_nudges.rds"
nudges_fl <- snakemake@input[["nudges"]]

#lm_fl <- "subworkflows/dgrp_coex/results/linear_models/male_model_01.collected-info.tsv.gz"
lm_fl <- snakemake@input[["lms"]]

#ol_fl <-"subworkflows/references/results/overlaps/overlaps.tsv.gz"
ol_fl <- snakemake@input[["ol"]]


# get significant models
lms <- read_tsv(lm_fl) %>% 
  filter(model == current_model) %>%
  filter(significant_x & significant_model) %>%
  dplyr::select(x=feature.x,coex.te=feature.y) %>% 
  distinct()

#tfs <- read_tsv("resources/Drosophila_melanogaster_TF.txt") %>% pull(Ensembl)

edges <- read_rds(edges_fl) %>% filter(x %in% lms$x) #%>% filter(x %in% tfs)
nudges <- read_rds(nudges_fl) %>% filter(x %in% edges$x)


# get overlapping features
# this is a symmetric matrix, elongated, ie all elements xy are also present as yx
ol <- read_tsv(ol_fl) %>% 
  filter(str_detect(Strain,"DGRP")) %>%
  dplyr::select(x=name.x,y=name.y,overlap) %>% 
  distinct()

# for each gene in lms that has a coex gene, annote by coex genes with
# overlap of the TE queried by each entry in lm
get_coex_set_overlaps <- . %>% 
  dplyr::select(x,coex.gene=y) %>%
  right_join(lms,.,by = c("x")) %>%
  filter(!is.na(coex.te)) %>% # don't consider pairs where the TE is not associated with the gene
  left_join(ol,by=c(coex.te="y",coex.gene="x")) %>%
  mutate(overlap=replace_na(overlap,F)) %>%
  filter(!is.na(coex.te))

edge_ols <- get_coex_set_overlaps(edges)
nudge_ols <- get_coex_set_overlaps(nudges)

# the number of coexpressed gene (coexpressed with putative TE regulator)
# a given TE overlaps
get_olap_cts <- . %>%
  mutate(overlap=if_else(overlap,"overlap","no.overlap")) %>%
  count(x,coex.te,overlap)

edge_olap_cts <- get_olap_cts(edge_ols)
nudge_olap_cts <- get_olap_cts(nudge_ols)

stopifnot((nudge_olap_cts$x %>% unique() %>% length())==(edge_olap_cts$x %>% unique() %>% length()))

fisher_tbls <- bind_rows(coex=edge_olap_cts,
                         other=nudge_olap_cts,.id="grp") %>%
  pivot_wider(names_from = grp, values_from = n)

fisher_tbls <- right_join(fisher_tbls,expand(fisher_tbls,x,coex.te,overlap)) %>% 
  group_by(x,coex.te) %>%
  filter(!(all(is.na(coex)) & all(is.na(other)))) %>%
  ungroup() %>%
  mutate(overlap=fct_relevel(overlap,"overlap")) %>%
  arrange(x,overlap)

fisher_tbls <- fisher_tbls %>%
  mutate(across(c("coex","other"),replace_na,0)) %>%
  nest(data=c(overlap,coex,other)) %>%
  mutate(data = map(data,column_to_rownames,"overlap"))

fisher_tbls <- fisher_tbls %>%
  mutate(fish = map(data,~broom::tidy(fisher.test(.x))))

fisher_tbls <- fisher_tbls %>%
  unnest(fish)

write_rds(fisher_tbls,snakemake@output[["rds"]])
library(tidyverse)

# raw_fl <- "results/quantification/vanilla_salmon_tes_transcripts/terminus.tx2group.raw.csv"
# tx2id_fl <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2id.tsv"
# tx2symbol_fl <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2symbol.tsv"
# tx2txsymbol_fl <- "results/references/transcripts_and_consensus_tes/transcripts_and_consensus_tes.tx2txsymbol.tsv"

raw_fl <- snakemake@input[["raw"]]
tx2id_fl <- snakemake@input[["tx2id"]]
tx2symbol_fl <- snakemake@input[["tx2symbol"]]
tx2txsymbol_fl <- snakemake@input[["tx2txsymbol"]]

raw <- read_csv(raw_fl, col_names = c("TXNAME","GROUPID"))
tx2id <- read_tsv(tx2id_fl) %>% set_names(c("TXNAME","GENEID"))
tx2symbol <- read_tsv(tx2symbol_fl) %>% set_names(c("TXNAME","GENESYMBOL"))
tx2txsymbol <- read_tsv(tx2txsymbol_fl) %>% set_names(c("TXNAME","TXSYMBOL"))

x <- left_join(raw, tx2id, by=c("TXNAME"))

x <- left_join(x, tx2symbol, by=c("TXNAME"))

x <- left_join(x, tx2txsymbol, by=c("TXNAME"))

# GROUPIDs for TEs are the long weird fasta header thing. If my clean TE name is in the GROUPID, I replace this.
x <- mutate(x, GROUPID2 = if_else(pmap_lgl(list(GENEID, GROUPID, TXNAME),
                                           .f=function(x,y,z) {grepl(x,y) & y==z}),
                                  GENEID,GROUPID))


# how many genes have transcripts split across single tx groups?
x <- x %>% group_by(GENEID, GROUPID) %>%
  add_count(name = "n_of_gene_in_group") %>%
  group_by(GENEID) %>%
  add_count(name = "n_of_gene_in_data") %>%
  group_by(GROUPID) %>%
  mutate(n_genes_in_group = length(unique(GENEID))) %>%
  ungroup()

# 1. groups containing only 1 gene that is entirely represented by the group
# naming this is easy
simple_groups <- group_by(x,GROUPID) %>%
  filter(all(n_of_gene_in_group==n_of_gene_in_data)) %>%
  filter(n_genes_in_group==1) %>%
  mutate(GROUPID_IDS = unique(GENEID), GROUPID_SYMBOLS = unique(GENESYMBOL)) %>%
  ungroup()

# 2. groups containing multiple genes where each gene is entirely represented by the group
# also pretty easy
complex_groups <- group_by(x,GROUPID) %>%
  filter(all(n_of_gene_in_group==n_of_gene_in_data)) %>%
  filter(n_genes_in_group>1) %>%
  mutate(GROUPID_IDS = paste(GROUPID2,paste(sort(unique(GENEID)),collapse = "_"),sep="_")) %>%
  mutate(GROUPID_SYMBOLS = paste(GROUPID2,paste(sort(unique(GENESYMBOL)),collapse = "_"),sep="_")) %>%
  ungroup()

# 3. groups containing single genes where each gene is represented by multiple groups, single-gene or otherwise
# gets harder, 
simple_group_multi <- group_by(x,GROUPID) %>%
  filter(n_genes_in_group==1) %>%
  filter(!GROUPID %in% simple_groups$GROUPID) %>%
  mutate(GROUPID_IDS = ifelse(n() == 1, GROUPID2, ifelse(n() <= 3,paste(GROUPID2,paste(sort(unique(TXNAME)),sep="_"),sep="_"),paste(GROUPID2,sort(unique(GENEID)),sep="_")))) %>%
  mutate(GROUPID_SYMBOLS = ifelse(n() == 1, GROUPID2, ifelse(n() <= 3,paste(GROUPID2,paste(sort(unique(TXSYMBOL)),sep="_"),sep="_"),paste(GROUPID2,sort(unique(GENESYMBOL)),sep="_")))) %>%
  ungroup()

# 4. groups containing multiple genes where each genes all representated by multiple groups
complex_group_multi <- group_by(x,GROUPID) %>%
  filter(n_genes_in_group>1) %>%
  filter(!GROUPID %in% complex_groups$GROUPID) %>%
  mutate(GROUPID_IDS = ifelse(n() <= 3,paste(GROUPID2,paste(sort(unique(TXNAME)),sep="_"),sep="_"),paste(GROUPID2,sort(unique(GENEID)),sep="_"))) %>%
  mutate(GROUPID_SYMBOLS =  ifelse(n() <= 3,paste(GROUPID2,paste(sort(unique(TXSYMBOL)),sep="_"),sep="_"),paste(GROUPID2,sort(unique(GENESYMBOL)),sep="_"))) %>%
  ungroup()

# sanity
cluster_list <- list(simple_groups, complex_groups, simple_group_multi, complex_group_multi)
new_clusters <- bind_rows(cluster_list)

stopifnot(nrow(new_clusters)==nrow(x))

stopifnot(all(new_clusters$GROUPID %in% x$GROUPID))

stopifnot(length(unique(new_clusters$GROUPID2)) == length(unique(x$GROUPID)))

stopifnot(length(unique(simple_groups$GROUPID_SYMBOLS)) == length(unique(simple_groups$GROUPID)))

stopifnot(length(unique(complex_groups$GROUPID_SYMBOLS)) == length(unique(complex_groups$GROUPID)))

stopifnot(length(unique(simple_group_multi$GROUPID_SYMBOLS)) == length(unique(simple_group_multi$GROUPID)))

stopifnot(length(unique(complex_group_multi$GROUPID_SYMBOLS)) == length(unique(complex_group_multi$GROUPID)))

stopifnot(length(unique(new_clusters$GROUPID_IDS)) == length(unique(x$GROUPID)))

stopifnot(length(unique(new_clusters$GROUPID_SYMBOLS)) == length(unique(x$GROUPID)))

# test if each old GROUPID has only a single equivalent ID in each of the new ID sets
new_clusters %>% 
  group_by(GROUPID) %>% 
  summarize(n_unique_new_groupid2 = length(unique(GROUPID2)),
            n_unique_new_groupid_ids = length(unique(GROUPID_IDS)),
            n_unique_new_groupid_symbols = length(unique(GROUPID_SYMBOLS))) %>%
  summarise(across(where(is.numeric), ~{all(.x==1)})) %>%
  gather(metric,value) %>%
  pull(value) %>%
  all() %>%
  stopifnot()

new_clusters %>%
  dplyr::select(TXNAME,GROUPID,GROUPID2,GROUPID_IDS,GROUPID_SYMBOLS) %>%
  write_tsv(snakemake@output[["tx2group"]])



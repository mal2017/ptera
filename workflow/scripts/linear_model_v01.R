library(tidyverse)
library(dbplyr)
library(rlang)

# params

# formula <- as.formula("y ~ x + wolbachia")
formula <- as.formula(snakemake@params[["formula"]])

# split_by <- c("sex")
split_by <- snakemake@params[["split_by"]]

# transforms <- list(y = "~log2(.x+10)")
transforms <- snakemake@params[["transforms"]]
message(paste("transforms:",transforms))

# filter_tes <- "sum(y > 0) == 800"
filter_tes <- snakemake@params[["te_filter"]]
filter_tes <- paste0("filter(",filter_tes,")")

# filter_genes<- "sum(x > 0) > 600 & sum(x > 7) > 100"
filter_genes <- snakemake@params[["gene_filter"]]
filter_genes <- paste0("filter(",filter_genes,")")

# db_file <- "results/quantification/vanilla_salmon_tes_transcripts/vst.sqlite"
db_file <- snakemake@input[["sqlite"]]

# open db connection
con <- DBI::dbConnect(RSQLite::SQLite(), db_file)

# just get features that pass the filter in the config
expressed_tes <- eval_tidy(expr(tbl(con,"transposons") %>% 
                                  group_by(transposon) %>% 
                                  !!parse_expr(filter_tes))) %>% 
  pull(transposon) %>% 
  unique()

expressed_genes <- eval_tidy(expr(tbl(con,"genes") %>% 
                                    group_by(host_gene) %>% 
                                    !!parse_expr(filter_genes))) %>% 
  pull(host_gene) %>% 
  unique()

sample_info <- tbl(con,"sample_info")

transposons <- tbl(con,"transposons") %>% filter(transposon %in% expressed_tes)

# last bit for testing only
genes <- tbl(con,"genes") %>% filter(host_gene %in% expressed_genes) #%>% filter(host_gene == "FBgn0085432") 

res <- left_join(sample_info, genes) %>%
  left_join(transposons) %>%
  collect()

# if transformations are in the params, loop through them and update the res.
if(!is.na(transforms)) {
  for (i in names(transforms)) {
    res <- dplyr::mutate(res,across(.cols = all_of(i),.fns=eval(parse(text=transforms[[i]]))))
  }
}

res <- res %>% 
  nest(data=-all_of(c("host_gene","transposon",split_by))) %>%
  mutate(fit = map(data, ~lm(formula, data = .))) %>%
  dplyr::select(-data) %>%
  mutate(tidied = map(fit,broom::tidy),
         glanced = map(fit,broom::glance),
         augmented = map(fit, broom::augment))

res %>%
  dplyr::select(sex, host_gene, transposon, glanced) %>%
  unnest(glanced) %>%
  vroom::vroom_write(snakemake@output[["glanced"]])

res %>%
  dplyr::select(sex, host_gene, transposon, tidied) %>%
  unnest(tidied) %>%
  vroom::vroom_write(snakemake@output[["tidied"]])

res %>%
  dplyr::select(sex, host_gene, transposon,augmented) %>%
  unnest(augmented) %>%
  vroom::vroom_write(snakemake@output[["augmented"]])

saveRDS(res %>% dplyr::select(sex,host_gene,transposon,fit) ,snakemake@output[["fits"]])


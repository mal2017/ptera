library(tidyverse)
library(dbplyr)
library(rlang)

# params

#genes_to_use <- "results/linear_models/y_~_x_+_wolbachia_by_sex_salmon_edaseq_qn/chunk_000" %>% read_tsv(col_names = "host_gene") %>% pull(host_gene)
genes_to_use <- snakemake@input[["genes"]] %>% read_tsv(col_names = "host_gene") %>% pull(host_gene)

# formula <- as.formula("y ~ x + wolbachia")
formula <- as.formula(snakemake@params[["formula"]])

# split_by <- c("sex")
split_by <- snakemake@params[["split_by"]]

# transforms <- ~{log2(.x+1)}
transforms <- snakemake@params[["transforms"]]
message(paste("transforms:",transforms))

# filter_tes <- "sum(y > 0) == 800"
filter_tes <- snakemake@params[["te_filter"]]
filter_tes <- paste0("filter(",filter_tes,")")

# filter_genes<- "sum(x > 0) > 600 & sum(x > 7) > 100"
filter_genes <- snakemake@params[["gene_filter"]]
filter_genes <- paste0("filter(",filter_genes,")")

#db_file <- "results/quantification/vanilla_salmon_tes_transcripts/gene.edaseq_qn.sqlite"
db_file <- snakemake@input[["sqlite"]]

message(paste("Data:",db_file))
message(paste("TE filter expression:",filter_tes))
message(paste("Gene filter expression:",filter_genes))
message(paste("Fitting separate models for:",split_by))
message(paste("Using formula:",format(formula)))
message(paste("Using transformations:",transforms))
message(paste("Applying to N genes:",length(genes_to_use)))

# diagnostics
#saveRDS(split_by, "~/split_by.rds")
#saveRDS(transforms, "~/transforms.rds")

# open db connection
con <- DBI::dbConnect(RSQLite::SQLite(), db_file)

# just get features that pass the filter in the config
expressed_tes <- eval_tidy(expr(tbl(con,"transposons") %>%
                                  group_by(transposon) %>%
                                  !!parse_expr(filter_tes))) %>%
  pull(transposon) %>%
  unique()

expressed_genes <- eval_tidy(expr(tbl(con,"genes") %>%
                                    filter(host_gene %in% genes_to_use) %>%
                                    group_by(host_gene) %>%
                                    !!parse_expr(filter_genes)
                             )) %>%
  pull(host_gene) %>%
  unique()

sample_info <- tbl(con,"sample_info")

transposons <- tbl(con,"transposons") %>% filter(transposon %in% expressed_tes)

# last bit for testing only
genes <- tbl(con,"genes") %>% filter(host_gene %in% expressed_genes) #%>% filter(host_gene == "FBgn0085432")

message("Pulling specified data into memory.")
res <- left_join(sample_info, genes) %>%
  left_join(transposons) %>%
  collect()

# if transformations are in the params, loop through them and update the res.
if(!is.na(transforms) & transforms!="NA") {
  for (i in names(transforms)) {
    message(paste("Transforming variable '",i,".'"))
    res <- dplyr::mutate(res,across(.cols = all_of(i),.fns=eval(parse(text=transforms[[i]]))))
  }
}

message("Running lm.")
if (!is.na(split_by) & split_by!="NA") {
  split_by <- c("host_gene","transposon",split_by)
} else {
  split_by <-  c("host_gene","transposon")
}
message(paste("Nesting by:",split_by))

# diagnostic
#save.image(file = "~/lm.RData")

res <- res %>%
  nest(data=-all_of(split_by)) %>%
  mutate(fit = map(data, ~lm(formula, data = .))) %>%
  dplyr::select(-data) %>%
  mutate(tidied = map(fit,broom::tidy),
         glanced = map(fit,broom::glance))
         #augmented = map(fit, broom::augment))

message("Exporting results.")
res %>%
  dplyr::select(split_by, glanced) %>%
  unnest(glanced) %>%
  vroom::vroom_write(snakemake@output[["glanced"]])

res %>%
  dplyr::select(split_by, tidied) %>%
  unnest(tidied) %>%
  vroom::vroom_write(snakemake@output[["tidied"]])

#res %>%
#  dplyr::select(split_by,augmented) %>%
#  unnest(augmented) %>%
#  vroom::vroom_write(snakemake@output[["augmented"]])

saveRDS(res %>% dplyr::select(split_by, fit) ,snakemake@output[["fits"]])

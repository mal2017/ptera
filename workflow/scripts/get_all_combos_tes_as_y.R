library(tidyverse)
library(vroom)

#dat <- vroom("subworkflows/dgrp_coex/results/linear_models/male_vst_control_copy_and_overlap/expression.tsv.gz", col_select = "feature")
dat <- vroom(snakemake@input[["expression"]])

#split_invoc <- "split -e -d -a 4 -n r/8 - chunk_"
split_invoc <- snakemake@params[["split_invoc"]]

tes <- dat %>% filter(!str_detect(feature,"FBgn")) %>% pull(feature)
genes <- dat %>% filter(str_detect(feature,"FBgn")) %>% pull(feature)

res <- expand_grid(tes = tes, genes=genes)

vroom_write(res, pipe(split_invoc), col_names = F)


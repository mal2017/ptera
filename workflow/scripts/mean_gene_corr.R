library(DescTools)
library(tidyverse)
library(corrr)

#fls <- c("subworkflows/dgrp_coex/results/linear_models/female/0/gene_x_gene.corrr.rds", "subworkflows/dgrp_coex/results/linear_models/female/1/gene_x_gene.corrr.rds","subworkflows/dgrp_coex/results/linear_models/female/2/gene_x_gene.corrr.rds")

fls <- snakemake@input

df <- map(fls, read_rds)

#df <- map(df,~{.x[1:10000,1:10001]}) %>% map_df(stretch, .id = "rep")

df <- map_df(df,stretch, .id = "rep")

res <- df %>% mutate(fishz = FisherZ(r))

res <- res %>%
  group_by(x, y) %>%
  summarise(mean_fishz = mean(fishz), .groups = "drop")

res <-  mutate(res,mean_r=FisherZInv(mean_fishz)) %>%
  dplyr::select(x, y, r = mean_r)

res <- retract(res,x,y,r) 

res <- as_cordf(res)

res <- shave(res)

write_rds(res, snakemake@output[["rds"]])
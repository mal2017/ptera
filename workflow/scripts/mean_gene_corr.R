library(DescTools)
library(tidyverse)
library(corrr)

fls <- snakemake@input

df <- map(fls, read_rds) %>%
  map_df(stretch, .id = "rep")

res <- df %>%
  mutate(fishz = FisherZ(r) ) %>%
  group_by(x, y) %>%
  summarise(mean_fishz = mean(fishz), .groups = "drop") %>%
  mutate(mean_r=FisherZInv(mean_fishz)) %>%
  dplyr::select(x, y, r = mean_r)

res <- retract(res,x,y,r) 

res <- as_cordf(res)

res <- shave(res)

write_rds(res, snakemake@output[["rds"]])
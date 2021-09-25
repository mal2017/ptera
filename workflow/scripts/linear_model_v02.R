library(tidyverse)
library(vroom)

#pairs_fl <- "results/linear_models/male_edaseq_qn/chunk_0000"
pairs_fl <- snakemake@input[["chunk"]]

pairs <- vroom::vroom(pairs_fl,col_names = c("feature.y","feature.x")) #%>% head(10)

#dat_fl <- "results/linear_models/male_edaseq_qn/expression.tsv.gz"
dat_fl <- snakemake@input[["dat"]]
dat <- vroom::vroom(dat_fl,num_threads = 1)

#cd_fl <- "results/meta/metadata.csv"
cd_fl <-  snakemake@input[["cd"]]
cd <- read_csv(cd_fl)

dat <- pivot_longer(dat,-feature,names_to = "sample",values_to = "score")

df <- left_join(pairs,dat, by=c(feature.y="feature"))

df <- left_join(df,dat, by=c(feature.x="feature",sample="sample"), suffix=c(".y",".x"))

df <- left_join(df,cd, by=c(sample="sample_name"))

## TODO: import overlaps for feature pairs

## TODO: import copy number for y (regressand) in each strain

df <- df %>% dplyr::rename(x="score.x", y="score.y")

# formula <- as.formula("y ~ x + wolbachia")
formula <- as.formula(snakemake@params[["formula"]])

res <- df %>%
  nest(-c(feature.x,feature.y)) %>%
  mutate(fit=map(data,~lm(formula, data = .))) %>%
  dplyr::select(-data) %>%
  mutate(tidy = map(fit,broom::tidy),
         glance = map(fit,broom::glance))

rm(df); rm(dat)

# ------------------- export ------------------------------------------
#tidy_fl <- "results/linear_models/male_edaseq_qn/chunk_0000.tidy.tsv"
#glance_fl <- "results/linear_models/male_edaseq_qn/chunk_0000.glance.tsv"

tidy_fl <- snakemake@output[["tidy"]]
glance_fl <- snakemake@output[["glance"]]
fits_fl <- snakemake@output[["fits"]]

res %>% dplyr::select(feature.x, feature.y, tidy) %>%
  unnest(tidy) %>%
  filter(term=="x") %>%
  vroom::vroom_write(tidy_fl)

res %>% dplyr::select(feature.x, feature.y, glance) %>%
  unnest(glance) %>%
  vroom::vroom_write(glance_fl)

fits <- res %>% dplyr::select(feature.x, feature.y, fit) %>%
  unite(relationship,feature.y, feature.x, sep="~") %>%
  deframe() %>%
  saveRDS(fits_fl)
  

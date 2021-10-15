library(tidyverse)
library(vroom)
library(modelr)


# -------- imort pairs for the chunk of lms to be evaluated --------
#pairs_fl <- "results/linear_models/male_vst/chunk_0000"
pairs_fl <- snakemake@input[["chunk"]]

pairs <- vroom::vroom(pairs_fl,col_names = c("feature.y","feature.x")) #%>% head(10)

# -------- prep data for the chunk of lms to be evaluated --------

#dat_fl <- "results/linear_models/male_vst/expression.tsv.gz"
dat_fl <- snakemake@input[["dat"]]
dat <- vroom::vroom(dat_fl,num_threads = 1)

#cd_fl <- "results/meta/metadata.csv"
cd_fl <-  snakemake@input[["cd"]]
cd <- read_csv(cd_fl)

dat <- pivot_longer(dat,-feature,names_to = "sample",values_to = "score")

df <- left_join(pairs,dat, by=c(feature.y="feature"))

df <- left_join(df,dat, by=c(feature.x="feature",sample="sample"), suffix=c(".y",".x"))

df <- left_join(df,cd, by=c(sample="sample_name"))

# -------------- add the overlap annotation to the dataframe ---------------

ol_fl <- "results/overlaps/overlaps.tsv.gz"
ol_fl <- snakemake@input[["ol"]]

ol <- read_tsv(ol_fl) %>%
  dplyr::select(Strain,feature.x=name.x,feature.y=name.y,overlap)

df <- df %>% left_join(ol, by = c("feature.y", "feature.x", "Strain")) %>%
  mutate(overlap = ifelse(is.na(overlap),F,overlap))

# --------begin to evaluate the lms --------

# rename these variables so the config can be more succinct
df <- df %>% dplyr::rename(x="score.x", y="score.y")

# formula <- as.formula("y ~ x + wolbachia + overlap")
formula <- as.formula(snakemake@params[["formula"]])

res <- df %>%
  nest(-c(feature.x,feature.y)) %>%
  mutate(fit=map(data,~lm(formula, data = .))) %>%
  #dplyr::select(-data) %>%
  mutate(tidy = map(fit,broom::tidy),
         glance = map(fit,broom::glance),
         aug = map2(fit, data, broom::augment_columns))

rm(df); rm(dat)

# ------------------- export ------------------------------------------
#tidy_fl <- "results/linear_models/male_edaseq_qn/chunk_0000.tidy.tsv"
#glance_fl <- "results/linear_models/male_edaseq_qn/chunk_0000.glance.tsv"

tidy_fl <- snakemake@output[["tidy"]]
glance_fl <- snakemake@output[["glance"]]
fits_fl <- snakemake@output[["fits"]]
aug_fl <- snakemake@output[["aug"]]

res %>% dplyr::select(feature.x, feature.y, tidy) %>%
  unnest(tidy) %>%
  filter(term!="(Intercept)") %>%
  vroom::vroom_write(tidy_fl)

res %>% dplyr::select(feature.x, feature.y, glance) %>%
  unnest(glance) %>%
  vroom::vroom_write(glance_fl)

fits <- res %>% dplyr::select(feature.x, feature.y, fit) %>%
  #unite(relationship,feature.y, feature.x, sep="~") %>%
  #deframe() %>%
  saveRDS(fits_fl)
  
# https://people.duke.edu/~rnau/testing.htm
# aug results will be usefull as well
aug <- res %>% dplyr::select(feature.y,feature.x,aug) %>%
  unnest(aug)

vroom::vroom_write(aug, aug_fl)

#aug %>% 
#  ggplot(aes(sample=.resid)) +
#  geom_qq() +
#  geom_qq_line() +
#  facet_wrap(feature.x ~ feature.y)
#  geom_point(aes(group=feature.y)) +
#  #geom_line(aes(group = feature.y), alpha = 1 / 3) + 
#  geom_smooth(se = FALSE) +
#  geom_rug() +
#  facet_wrap(~feature.y)

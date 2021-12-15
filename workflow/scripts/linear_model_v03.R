library(tidyverse)
library(tls)

# -------- imort pairs for the chunk of lms to be evaluated --------
#pairs_fl <- "subworkflows/dgrp_coex/results/linear_models/female_model_01/chunk_0080"
pairs_fl <- snakemake@input[["chunk"]]

pairs <- vroom::vroom(pairs_fl,col_names = c("feature.y","feature.x")) #%>% head(10)

# -------- prep data for the chunk of lms to be evaluated --------

#dat_fl <- "subworkflows/dgrp_coex/results/linear_models/female_model_01/expression.tsv.gz"
dat_fl <- snakemake@input[["dat"]]
dat <- vroom::vroom(dat_fl,num_threads = 1)

#cd_fl <- "subworkflows/dgrp_coex/results/meta/metadata.csv"
cd_fl <-  snakemake@input[["cd"]]
cd <- read_csv(cd_fl)

dat <- pivot_longer(dat,-feature,names_to = "sample",values_to = "score")

# ----------- join with metadata and pairs to create working dataframe -----------
df <- left_join(pairs,dat, by=c(feature.y="feature"))

df <- left_join(df,dat, by=c(feature.x="feature",sample="sample"), suffix=c(".y",".x"))

df <- left_join(df,cd, by=c(sample="sample_name"))

# -------------- add the overlap annotation to the dataframe ---------------

#ol_fl <- "results/overlaps/overlaps.tsv.gz"
ol_fl <- snakemake@input[["ol"]]

ol <- read_tsv(ol_fl) %>%
  dplyr::select(Strain,feature.x=name.x,feature.y=name.y,overlap)

df <- df %>% left_join(ol, by = c("feature.y", "feature.x", "Strain")) %>%
  mutate(overlap = ifelse(is.na(overlap),F,overlap))

# ----------- add the copies annotation to the dataframe

#TEST: copies_fl <- "subworkflows/dgrp_wgs/results/copies/copies.tsv"
copies_fl <- snakemake@input[["copies"]]

copies <- read_tsv(copies_fl) %>%
  group_by(sequence) %>%
  mutate(scaled.copies = scale(est.copies)[,1]) %>%
  ungroup()

# join estimated copies and set to 1 where no info is available, (0 for scaled)
df <- df %>%
  left_join(copies, by=c(Strain="Strain",feature.x="sequence")) %>%
  left_join(copies, by=c(Strain="Strain",feature.y="sequence")) %>%
  mutate_at(c("est.copies.y","est.copies.x"), replace_na, 1) %>%
  mutate_at(c("scaled.copies.x","scaled.copies.y"),replace_na,0)

# ----------------remove outliers per feature by mads ----------------------
#mads_filter <- 3
mads_filter <- snakemake@params[["mads_filter"]]

dat <- dat %>% 
  group_by(feature) %>% 
  mutate(mads = abs(score- median(score))/mad(score)) %>%
  ungroup() %>%
  filter(mads < mads_filter)

# --------begin to evaluate the lms --------

# rename these variables so the config can be more succinct
df <- df %>% dplyr::rename(x="score.x", y="score.y")

# main formula
# TESTING: formula <- "y ~ 0 + x + wolbachia + scaled.copies.y + overlap"
formula <- snakemake@params[["formula"]]

# alternate formula; mainly relevant when using the overlap term, which is frequently all F for any
# given TE/gene pair.
# TESTING: alt_formula <- "y ~ 0 + x + wolbachia + scaled.copies.y"
alt_formula <- snakemake@params[["alt_formula"]]

# custom function for tidying tls::tls fit objs.
tidy2 <-  function(x) {
  full_join(enframe(x$coefficient, name = "term", value = "estimate"),
            as_tibble(x$confidence.interval,rownames = "term"),by="term") %>%
    full_join(enframe(x$sd.est, name = "term", value = "sd.est"),
              by="term")
}

# apply tls and tidy the results. alternate formula is possible because overlaps are often all False, leading to
# tls error. This allows us to pull out all situations where overlap strongly correlates to paired expression.
res <- df %>%
  nest(-c(feature.x,feature.y)) %>%
  mutate(fit=map(data,~tls::tls(as.formula(if_else(length(unique(.$overlap))==1,alt_formula, formula)), data = ., method = "normal"))) %>%
  mutate(tidy = map(fit,tidy2))

# ------------------- export ------------------------------------------
tidy_fl <- snakemake@output[["tidy"]]

res %>% dplyr::select(feature.x, feature.y, tidy) %>%
  unnest(tidy) %>%
  vroom::vroom_write(tidy_fl)

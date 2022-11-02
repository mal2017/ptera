library(tidyverse)
library(lmtest)

# -------- imort pairs for the chunk of lms to be evaluated --------
#pairs_fl <- "subworkflows/dgrp_coex/results/linear_models/male_model_01/2/chunk_0271"
pairs_fl <- snakemake@input[["chunk"]]

pairs <- vroom::vroom(pairs_fl,col_names = c("feature.y","feature.x")) #%>% head(10)

# -------- prep data for the chunk of lms to be evaluated --------

#dat_fl <- "subworkflows/dgrp_coex/results/linear_models/male_model_01/2/expression.tsv.gz"
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

#ol_fl <- "subworkflows/references/results/overlaps/overlaps.tsv.gz"
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

#unique(copies$sequence)[!unique(copies$sequence) %in% df$feature.y]
copies <- copies %>% mutate(sequence = str_replace(sequence,"\\\\","_")) #%>% filter(str_detect(sequence,"Ddip"))

stopifnot(all(df$feature.y %in% copies$sequence))

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
# TESTING: alt_formula <- "y ~ 0 + x + wolbachia"
alt_formula <- snakemake@params[["alt_formula"]]


# execute the modeling and tests of signficance -------------------------------
possibly_tidy2 <- possibly(broom::tidy,NULL)
possibly_glance2 <- possibly(broom::glance,NULL)
possibly_anova2 <- possibly(~broom::tidy(drop1(.,test="F")),NULL)

# custom func for running lm
run_lm <- function(.x,.y) {
  # alternate formula is possible because overlaps are often all False, leading to
  # error. This allows us to pull out all situations where overlap strongly 
  # correlates to paired expression.
  message(.y);
  frm <- as.formula(if_else(length(unique(.x$overlap))==1,alt_formula, formula))
  lm(frm, data = .x)
}

possibly_run_lm <- possibly(run_lm,NULL)

res <- df %>% 
  nest(-c(feature.x,feature.y)) %>% 
  mutate(fit=map2(data,feature.x,possibly_run_lm)) %>%
  mutate(tidy = map(fit,possibly_tidy2),
         glance = map(fit,possibly_glance2),
         anova = map(fit,possibly_anova2))

# tests of LM assumptions -----------------------------------------------------

possibly_bp2 <- possibly(~broom::tidy(bptest(.)),NULL) # heteroskedasticity
# possibly_dw2 <- possibly(~broom::tidy(dwtest(.)),NULL) # dependence of errors, but not meaningful for non-ordered data
possibly_rain2 <- possibly(function(dat,y) broom::tidy(raintest(y,order.by = ~x, data=dat)),NULL) # linearity


res <- res %>%
  mutate(breusch_pagan = map(fit,possibly_bp2),
         rainbow = map2(data,fit,possibly_rain2))


# ------------------- export ------------------------------------------
tidy_fl <- snakemake@output[["tidy"]]
glance_fl <- snakemake@output[["glance"]]
anova_fl <- snakemake@output[["anova"]]
rainbow_fl <- snakemake@output[["rainbow"]]
breusch_pagan_fl <- snakemake@output[["breusch_pagan"]]

res %>% dplyr::select(feature.x, feature.y, tidy) %>%
  unnest(tidy) %>%
  vroom::vroom_write(tidy_fl)

res %>% dplyr::select(feature.x, feature.y, glance) %>%
  unnest(glance) %>%
  vroom::vroom_write(glance_fl)

res %>% dplyr::select(feature.x, feature.y, anova) %>%
  unnest(anova) %>%
  vroom::vroom_write(anova_fl)

res %>% dplyr::select(feature.x, feature.y, rainbow) %>%
  unnest(rainbow) %>%
  vroom::vroom_write(rainbow_fl)

res %>% dplyr::select(feature.x, feature.y, breusch_pagan) %>%
  unnest(breusch_pagan) %>%
  vroom::vroom_write(breusch_pagan_fl)
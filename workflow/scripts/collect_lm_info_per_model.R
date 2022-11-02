library(tidyverse)

#info_fls <- Sys.glob("data/linear_models/female_model_01/*/lm.collected-info.tsv.gz")
info_fls <- snakemake@input[["info"]]

info <- info_fls %>% map_df(read_tsv)

# annotate as valid if it passes breusch pagan and rainbow tests 
# and the overall fit is significant
info <- info %>%
  mutate(valid = p.value_breuschpagan > 0.05 & p.value_rainbow > 0.05 & p.value_ftest_r2 < 0.05)
  
# note as reproducible if all reps suggest nonsignificance or all reps are
# significant and have same signs for gene expression coefs
info <- info %>%
  group_by(model,feature.x,feature.y) %>%
  mutate(reproducible = all(valid) & 
           (
             all(p.value_anova_x > 0.05) | 
               (all(p.value_anova_x < 0.05) & length(unique(sign(estimate)))==1)
             )
         )

# for each gene/te pair, take the least significant and lowest effect size
# model rep as the representative, to be conservative.
info <- info %>%
  group_by(model,feature.x,feature.y) %>%
  slice_max(p.value_anova_x, n=1) %>%
  slice_min(abs(estimate.qnorm),n=1,with_ties = F)
  
# now adjust p values once we have 1 model per te/gene pair
info <- info %>%  
  group_by(model) %>%
  mutate(across(starts_with("p.value_anova"),
                ~p.adjust(.x,method="BH"),
                .names= "adj_{.col}")) %>%
  ungroup()
    
info <- info %>%  
  mutate(significant = adj_p.value_anova_x < 0.1 & reproducible & valid)

gene_universe <- read_tsv("http://ftp.flybase.net/releases/FB2022_04/precomputed_files/genes/fbgn_fbtr_fbpp_expanded_fb_2022_04.tsv.gz", skip = 4)

lkup <- gene_universe %>% dplyr::select(gene_ID, gene_symbol) %>% distinct()

info <- info %>% 
  group_by(model) %>%
  mutate(coef.quantile = cume_dist(abs(estimate.qnorm))) %>%
  ungroup() %>%
  left_join(lkup, by=c(feature.x = "gene_ID")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol),feature.x,gene_symbol)) %>%
  relocate(gene_symbol, .after = "feature.x")


write_tsv(info,snakemake@output[["tsv"]])
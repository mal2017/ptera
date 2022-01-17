library(scater)
library(scran)
library(tidyverse)
library(scuttle)
library(PCAtools)

source("../../workflow/scripts/ggplot_theme.R")

# see http://bioconductor.org/books/3.14/OSCA.advanced/

#sce_fl <- "results/downstream/initial_filtering/FCA62_Female_ovary_adult_5dWT_Nystul_All_Nuclei_S57.usa.filt.sce.rds"

sce_fl <- snakemake@input[["sce"]]

sce <- read_rds(sce_fl)

# initial filtering is previously done to ensure the same starting cells
# regardless of choice of downstream processing tools.

# -------------------------------------------------------------------------------------
# Highly Variable Gene Selection
# -------------------------------------------------------------------------------------

# now will model variance to get hivar genes
# all the same sample, no experimental conditions, so no need to add design

counts(sce) <- as.matrix(counts(sce))
logcounts(sce) <- as.matrix(logcounts(sce))

dec <- modelGeneCV2(sce, max.iter=200,min.mean=0.01)

fit <- metadata(dec)

# this would be the place to exclude TEs, but I think it is
# actually important to leave them in - this may cause other
# clusters to appear that are defined in part by TE expression,
# and these would be worth digging into to see if there is
# real biology there or artifactual stuff.

# here I get HVGs the staandard variance approach
hvg <- getTopHVGs(dec,n=2000,var.field = "ratio")

dat_var <- tibble(feature = names(fit$cv2),
                  mean = fit$mean,
                  metric = fit$cv2,
                  trend = fit$trend(mean))

g_var <- ggplot(dat_var, aes(mean,metric)) +
  geom_point() +
  geom_point(data = . %>% filter(feature %in% hvg),color="blue") +
  scale_y_log10() +
  scale_x_log10() +
  geom_line(aes(mean,trend),color="red") +
  #coord_cartesian(ylim=c(0.1,max(dat_var$var,na.rm = T)*1.1)) +
  theme(aspect.ratio = 0.5)

# add this info to the feature df.
rowData(sce)$is_hvg_in_sample <- rownames(rowData(sce)) %in% hvg

# -------------------------------------------------------------------------------------
# choosing N PCs
# -------------------------------------------------------------------------------------

sce <- runPCA(sce, subset_row = hvg, scale=T)

percent.var <- attr(reducedDim(sce), "percentVar")

chosen.elbow <- PCAtools::findElbowPoint(percent.var)

dat_elbow <- tibble(PC=1:length(percent.var),`Variance explained (%)`=percent.var/100)

g_elbow <- ggplot(dat_elbow, aes(PC,`Variance explained (%)`)) +
  geom_point() +
  geom_vline(xintercept =chosen.elbow,color="red",linetype="dashed") +
  xlab("PC") +
  ylab("Variance explained") +
  scale_y_continuous(labels=scales::percent)

# retain the full set, but also subset
reducedDim(sce, "PCA.elbow") <- reducedDim(sce)[,1:chosen.elbow]

plotPCASCE(sce)

# -------------------------------------------------------------------------------------
# More dimensionality reduction for checking
# -------------------------------------------------------------------------------------

#set.seed(2)
#sce <- runTSNE(sce,dimred="PCA.elbow")
#plotReducedDim(sce,"TSNE")

#set.seed(2)
#sce <- runUMAP(sce,dimred="PCA.elbow")
#plotReducedDim(sce,"UMAP")

# -------------------------------------------------------------------------------------
# Export results
# -------------------------------------------------------------------------------------

write_rds(sce,snakemake@output[["sce"]])
write_rds(dec,snakemake@output[["dec"]])

write_tsv(dat_var,snakemake@output[["dat_hvg"]])
saveRDS(g_var,snakemake@output[["g_hvg"]])
ggsave(snakemake@output[["png_hvg"]],g_var, width=4, height = 4)

write_tsv(dat_elbow,snakemake@output[["dat_elbow"]])
saveRDS(g_elbow,snakemake@output[["g_elbow"]])
ggsave(snakemake@output[["png_elbow"]],g_elbow, width=4, height = 4)

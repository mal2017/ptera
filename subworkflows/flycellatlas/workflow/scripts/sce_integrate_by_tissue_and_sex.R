library(scater)
library(scran)
library(tidyverse)
library(batchelor)

source("../../workflow/scripts/ggplot_theme.R")

# http://bioconductor.org/books/3.14/OSCA.multisample/integrating-datasets.html#quick-start

# -----------------------------
#  Get previously calculated HVG results
# -----------------------------

#dec_fls <- Sys.glob("results/downstream/single_sample_dimred/*Male*testis*.dec.rds")
dec_fls <- snakemake@input[["decs"]]

names(dec_fls) <- str_extract(dec_fls,"(?<=dimred\\/).+(?=\\.usa)")

decs <- dec_fls %>% map(read_rds)

# -----------------------------
#  Get the SCE opjects for each batch
# -----------------------------

#sce_fls <- Sys.glob("results/downstream/single_sample_clustering/*Male*testis*.sce.rds")
sce_fls <- snakemake@input[["sces"]]

names(sce_fls) <- str_extract(sce_fls,"(?<=clustering\\/).+(?=\\.usa)")

sces <- sce_fls %>% map(read_rds)

# -----------------------------
#  Get HVGs from each batch
# -----------------------------

#hvgs <- sces %>%
#  map(~{subset(rowData(.x),is_hvg_in_sample)}) %>%
#  map(rownames) %>%
#  unlist() %>%
#  unique()

#all_features <- sces %>% map(rownames) %>% unlist() %>% unique()

#hvgs <- hvgs[hvgs %in% all_features]

# -----------------------------
#  Apply mnn via batchelor
# -----------------------------

# currently, correction is
corrected <- quickCorrect(sces,
                          precomputed = decs,
                          hvg.args = list(var.field="ratio"), # for use with cv2
                          PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam()))

sce <- corrected$corrected

colnames(sce) <- make.unique(colnames(sce))

sce <- runTSNE(sce,dimred="corrected")
sce <- runUMAP(sce,dimred="corrected")

g_batches <- plotTSNE(sce,colour_by="batch")


saveRDS(sce,snakemake@output[["sce"]])
saveRDS(g_batches,snakemake@output[["ggp"]])
write_tsv(g_batches$data,snakemake@output[["dat"]])
ggsave(snakemake@output[["png"]],g_batches)

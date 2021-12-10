library(rtracklayer)
library(tidyverse)

#ins_fl <- "results/overlaps/polymorphic_insertions/dl/DGRP_flies/RAL-101_result/RAL-101_Inserts_Annotated.txt"
#deps_fl <- "results/overlaps/polymorphic_insertions/dl/DGRP_flies/RAL-101_result/RAL-101_Depletion_Annotated.txt"
#ref_fl <- "results/repeatmasker/reference_insertions.bed"

ins_fl <- snakemake@params[["poly_ins_tsv"]]
deps_fl <- snakemake@params[["poly_deps_tsv"]]
ref_fl <- snakemake@input[["ref_copies"]]

ins <- read_tsv(ins_fl)
deps <- read_tsv(deps_fl)

# ----------------------------------- process ref ------
# mostly about getting the names to make sense
ref <- rtracklayer::import(ref_fl)

raw_te_names <- ref$name
# wtf is univode or escape character between TIDAL TE names after the 3rd vert bar????
# x0 should correspond to the name salmon gives the te in the output file
x0 <- gsub("[\t\n\r\v\f\a\b ]","",x=gsub("(?<=[\t\n\r\v\f\a\b ]).+","",x = raw_te_names, perl = T))
x <- gsub(".+\\|","", x=x0)
x1 <- gsub("#.+","",x = x)

ref$name <- gsub("gb\\|.+\\|","",x = x1)

ref$source <- "reference"

# -------------------------------- process ins and deps

# http://www.bio.brandeis.edu/laulab/Tidal_Fly/UserGuide_TIDAL_outputs.html

ins <- ins %>% 
  mutate(strand=ifelse(TE_coord_start > TE_coord_end,"-","+")) %>%
  dplyr::select(chr=Chr,start=Chr_coord_5p,end=Chr_coord_3p,name=TE, strand) %>%
  mutate(start = start + 1) %>%
  mutate(chr=str_remove(chr,"chr"))%>%
  mutate(source = "Tidal") %>%
  GRanges()

deps <- deps %>% 
  filter(Chr_5p...2 == Chr_3p...24) %>%
  filter(Chr_5p...2 == Chr_mid) %>%
  dplyr::select(chr=Chr_5p...2,start=Chr_coord_5p_end,end=Chr_coord_3p_start,name=repName) %>%
  mutate(start_old = start, end_old = end, start=ifelse(start_old > end_old,end,start),end=ifelse(start_old > end_old,start_old,end)) %>%
  mutate(start = start + 1) %>%
  mutate(chr=str_remove(chr,"chr"))%>%
  GRanges()

# these match with the regions from the tidal beds, but are importable to R,
# spot checked in IGV for RAL-101
#export(ins,"~/Downloads/test_ral101_ins.bed")
#export(deps,"~/Downloads/test_ral101_deps.bed")

# ------------------ remove depletions from ref set ----------------------
# https://support.bioconductor.org/p/87081/

retained <- subsetByOverlaps(ref,deps,minoverlap = 1, type = "any", invert = T)

#export(retained,"~/Downloads/retained.bed")

retained_and_gained <- c(retained, ins)

retained_and_gained$score <- 1

retained_and_gained$strain <- snakemake@wildcards[["tidal_strain"]]
retained_and_gained$tidal_group <- snakemake@wildcards[["tidal_group"]]

# only matters if DGRP
retained_and_gained$strain <- str_replace(retained_and_gained$strain,"RAL-","DGRP_")

#
#export(retained_and_gained,"~/Downloads/retained_and_gained.bed")

write_tsv(as_tibble(retained_and_gained),snakemake@output[["tsv"]])

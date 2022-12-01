library(tidyverse)

#fls <- Sys.glob("subworkflows/references/results/overlaps/polymorphic_insertions/DGRP_flies/*tsv")
fls <- snakemake@input

x <- map_df(fls,read_tsv)

#important because TIDAL includes multiple libraries for the same strain.
if (str_detect(snakemake@wildcards[["tidal_group"]],"DGRP")) {
    x <- mutate(x,strain  = str_extract(strain,"DGRP_\\d+"))
}

x <- x %>% distinct()

write_tsv(x,snakemake@output[["tsv"]])
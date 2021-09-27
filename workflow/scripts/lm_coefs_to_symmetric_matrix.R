library(tidyverse)
library(vroom)

tidy_fl <- snakemake@input[["tidy"]]

res.tidy <- vroom(tidy_fl, num_threads = 1)

features.x <- unique(res.tidy$feature.x)
features.y <- unique(res.tidy$feature.y)

unique.features <- sort(unique(c(features.x, features.y)))

mat <- matrix(nrow=length(unique.features), ncol=length(unique.features),data = numeric(0))

colnames(mat) <- unique.features; rownames(mat) <- unique.features

x <- res.tidy$feature.x
y <- res.tidy$feature.y
z <- res.tidy$estimate

for (i in 1:nrow(res_vst.tidy)) {
  mat[x[i], y[i]] <- z[i]
  mat[y[i], x[i]] <- z[i]
}

stopifnot(isSymmetric(mat))

mat_df <- as_tibble(mat,rownames="feature")

write_tsv(mat_df,snakemake@output[["tsv"]])

#!/bin/Rscript
# izpiše meje domen za dani protein
source(here::here("scripts", "utils.r"))

d <- read.csv(paths$domains)

p <- {
    arg <- commandArgs(trailingOnly = TRUE)[1]
    if (is.na(arg)) stop("./domains_inds.r <protein-name>\n")
    arg
}

out <- d[d$protein == p, -1] |> unlist()
cat(out[1], out[2], out[3], out[4], "\n")

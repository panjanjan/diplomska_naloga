#!/bin/Rscript
# izpiše meje domen za dani protein
setwd(Sys.getenv("ROOT"))
source("./scripts/utils.r")

d <- load_data("domains")
p <- {
    arg <- commandArgs(trailingOnly = TRUE)[1]
    if (is.na(arg)) stop("manjka ime proteina: ./domains_inds.r <protein>\n")
    arg
}

out <- d[d$protein == p, -1] |> unlist()
cat(out[1], out[2], out[3], out[4], "\n")

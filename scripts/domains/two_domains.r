#!/bin/Rscript
# naredi csv, ki hrani informacije o domenah
#
# protein1 start1 end1 start2 end2
# protein2 start1 end1 start2 end2
# protein3 start1 end1 start2 end2
# ...
library(magrittr)
library(here)
source("scripts/utils.r")

out <- file.path(root, "two_domains.csv")
data <- read.csv(file.path(root, "sword_results_clean.csv"))

# izloči tiste, ki imajo samo 2 domeni
# ohrani ime proteina in meje domen
data <- data$protein %>%
    table() %>%
    {which(. < 3)} %>%
    names() %>%
    {data[data$protein %in% ., c("protein", "start", "end")]}

# združi vrstice, da bo en protein na vrstico
d2 <- data.frame()

for (i in seq(1, nrow(data), 2)) {
    dsub  <- data[i:(i+1), ]
    start <- dsub$start
    end   <- dsub$end
    dnew  <- data.frame(
        protein = dsub$protein[1],
        start1  = start[1],
        end1    = end[1],
        start2  = start[2],
        end2    = end[2]
    )
    d2 <- rbind(d2, dnew)
}

write.csv(d2, out, quote = FALSE, row.names = FALSE)

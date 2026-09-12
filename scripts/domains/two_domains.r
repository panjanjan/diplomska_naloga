#!/bin/Rscript
# naredi csv, ki hrani informacije o domenah
#
# protein1 start1 end1 start2 end2
# protein2 start1 end1 start2 end2
# protein3 start1 end1 start2 end2
# ...
# -----------------------------------------------------------------
library(magrittr)

source(here::here("scripts", "utils.r"))

# -----------------------------------------------------------------
cat("reading\n")
data <- read.csv(paths$sword_clean)

# ohrani proteine z 2 domenama. število vrstic v CSVju določa št.
# domen. glej ./build_decompositions_csv.py
# ohrani ime proteina in meje domen
cat("trimming\n")
data <- data$protein %>%
    table() %>%
    {which(. < 3)} %>%
    names() %>%
    {data[data$protein %in% ., c("protein", "start", "end")]}

# združi vrstice, da bo en protein na vrstico
cat("building\n")
d2 <- data.frame()

# ker vemo da bosta samo 2 vrstici na protein
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

cat("writing\n")
write.csv(d2, paths$domains, quote = FALSE, row.names = FALSE)



# NOTE: Well..... lahko bi še bolj zmanjšal redundanco in shranil
# samo end1 ali start2. Vedno se začne z 1, vedno bo end2 na koncu
# proteina. Ampak I CBA.

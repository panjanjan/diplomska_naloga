#!/bin/Rscript
# dodaj podatek o verigi v PDB datoteke iz ATLAS
# pomembno za delovanje SWORD2
library(stringr)
library(bio3d)
setwd(Sys.getenv("ROOT"))

old_dir <- "atlas_db/PDB"
new_dir <- "atlas_db/PDB_chained"
if (!dir.exists(new_dir)) dir.create(new_dir)

pdbs <- load_data("pdb")

for (file in pdbs) {
     cat(file, ": ")

     # preberi datoteko in dobi njeno ime iz absolutne poti
     pdb <- read.pdb(file)
     fname <- basename(file)

     # iz imena povleči podatek o verigi
     chain <- str_replace(fname, "\\w+_(.*).pdb", "\\1")

     # vstavi verigo v PDB
     pdb$atom$chain <- rep(chain, nrow(pdb$atom))

     # izpiši v novo datoteko
     new_fname <- file.path(new_dir, fname)
     write.pdb(pdb, new_fname)

     cat("done\n")
}

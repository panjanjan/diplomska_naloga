#!/bin/Rscript
# WARN: poženi da vidiš če dela pravilno

# dodaj podatek o verigi v PDB datoteke iz ATLAS
# pomembno za delovanje SWORD2
library(stringr)
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

# -----------------------------------------------------------------
n_cores <- min(detectCores() - 1, 10)

old_dir <- paths$pdb
new_dir <- here::here("atlas_db", "PDB_chained")
if (!dir.exists(new_dir)) dir.create(new_dir)

pdbs <- list.files(paths$pdb, full.names = TRUE)

mclapply(pdbs, \(file) {
     # preberi datoteko in dobi njeno ime iz absolutne poti
     pdb <- read.pdb(file)
     fname <- basename(file)

     # iz imena povleči podatek o verigi
     chain <- str_replace(fname, "\\w+_(.*).pdb", "\\1")

     # vstavi verigo v PDB
     pdb$atom$chain <- rep(chain, nrow(pdb$atom))

     # izpiši v novo datoteko
     new_fname <- here::here(new_dir, fname)
     # write.pdb(pdb, new_fname)
     print(new_fname)
}, mc.cores = n_cores)

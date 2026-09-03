#!/bin/Rscript
# doda podatek o verigi v PDB datoteke iz ATLAS
# pomembno za delovanje SWORD2
library(stringr)
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

# -----------------------------------------------------------------
n_cores <- min(detectCores() - 1, 10)

# trenutne PDB datoteke bivajo tukaj
old_dir <- paths$pdb

# nove PDB datoteke bodo bivale tukaj
new_dir <- here::here("atlas_db", "PDB_chained")
if (!dir.exists(new_dir)) dir.create(new_dir)

# datoteke
pdbs <- list.files(paths$pdb, full.names = TRUE)

print(commandArgs(trailingOnly = TRUE)[1])

# naredi na par datotekah če testiraš
arg <- commandArgs(trailingOnly = TRUE)[1]
if (!is.na(arg) && arg == "-t") {
  new_dir <- "TEST_CHAINS"
  pdbs <- head(pdbs)
}

if (!dir.exists(new_dir)) dir.create(new_dir)

invisible(mclapply(pdbs, \(file) {
    # preberi PDB datoteko in iz imena dobi podatek o verigi
    pdb <- read.pdb(file)
    fname <- basename(file)
    chain <- str_replace(fname, "\\w+_(.*).pdb", "\\1")

    # vstavi v PDB
    pdb$atom$chain <- rep(chain, nrow(pdb$atom))

    # izpiši v novo datoteko
    new_fname <- here::here(new_dir, fname)
    print(new_fname)
    write.pdb(pdb, new_fname)
}, mc.cores = n_cores))

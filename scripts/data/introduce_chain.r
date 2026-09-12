#!/bin/Rscript
# doda podatek o verigi v PDB datoteke iz ATLAS
# pomembno za delovanje SWORD2
library(stringr)
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

# -----------------------------------------------------------------
n_cores <- min(detectCores() - 1, 10)

# trenutne PDB datoteke
pdbs <- list.files(paths$pdb_old, full.names = TRUE)

# nove PDB datoteke
if (!dir.exists(paths$pdb)) dir.create(paths$pdb)

invisible(mclapply(pdbs, \(file) {
    # preberi PDB datoteko in iz imena dobi podatek o verigi
    pdb <- read.pdb(file)
    fname <- basename(file)
    chain <- str_replace(fname, "\\w+_(.*).pdb", "\\1")

    # vstavi v PDB
    pdb$atom$chain <- rep(chain, nrow(pdb$atom))

    # izpiši v novo datoteko
    new_fname <- here::here(paths$pdb, fname)
    cat(new_fname, "\n")
    write.pdb(pdb, new_fname)
}, mc.cores = n_cores))

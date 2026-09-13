#!/bin/Rscript
# Izračuna RMSF vrednosti za vsak protein in njegove replikate
# ter jih zapiše v CSV datoteke.
#
# Za izračune uporabi 3 selekcije atomov:
# vsi atomi brez vodikov, atomi verige (backbone) in samo C-alfa
# atomi.
#
# Za vsak protein ustvari 3 CSVje, saj imajo različne selekcije
# različno število vrnjenih RMSF vrednosti.
# --------------------------------------------------------------
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

# --------------------------------------------------------------
# število jeder za paralelizacijo
n_cores <- min(detectCores() - 1, 10)

# PDB datotetke, DCD trajektorije, proteini
data  <- list(
    pdb     = list.files(paths$pdb,  pattern = ".pdb", full.names = TRUE),
    traj    = list.files(paths$traj, pattern = ".dcd", full.names = TRUE),
    domains = read.csv(paths$domains)
)

n_all <- nrow(data$domains)
n_replicates <- 3

# --------------------------------------------------------------
run <- function(i) {
    protein <- data$domains$protein[i]
    # cat("[", i, "/", n_all, "] ", protein, " ... ", sep = "")

    dcdfiles <- grep(protein, data$traj, value = TRUE)
    pdbfile  <- grep(protein, data$pdb, value = TRUE)
    print(paste(dcdfiles, pdbfile))
    assertthat::are_equal(length(dcdfiles), n_replicates)
    assertthat::are_equal(length(pdbfile), 1)

    pdb <- read.pdb(pdbfile, verbose = FALSE)

    # ustvari 3 selekcije
    sel_list <- list(
        all = atom.select(pdb, "noh"),
        bb  = atom.select(pdb, "backbone"),
        ca  = atom.select(pdb, "calpha")
    )

    # izračunaj RMSFje za vse replikate
    r1 <- run_replicate(dcdfiles[1], pdb, sel_list)
    r2 <- run_replicate(dcdfiles[2], pdb, sel_list)
    r3 <- run_replicate(dcdfiles[3], pdb, sel_list)

    # združi in shrani vrednosti po selekcijah
    out <- file.path(paths$rmsf, paste0(protein, "_rmsf_all.csv"))
    d <- data.frame(
        R1 = r1$all,
        R2 = r2$all,
        R3 = r3$all
    )
    write.csv(d, out, quote = FALSE, row.names = FALSE)

    out <- file.path(paths$rmsf, paste0(protein, "_rmsf_bb.csv"))
    d <- data.frame(
        R1  = r1$bb,
        R2  = r2$bb,
        R3  = r3$bb
    )
    write.csv(d, out, quote = FALSE, row.names = FALSE)

    out <- file.path(paths$rmsf, paste0(protein, "_rmsf_ca.csv"))
    d <- data.frame(
        R1  = r1$ca,
        R2  = r2$ca,
        R3  = r3$ca
    )
    write.csv(d, out, quote = FALSE, row.names = FALSE)

    # cat("done\n")
}

# vrne seznam RMSF vrednosti za vse tri selekcije
run_replicate <- function(dcdfile, pdb, sel_list) {
    dcd <- read.dcd(dcdfile, verbose = FALSE)
    list(
         all = rmsf(dcd[, sel_list$all$xyz]),
         bb  = rmsf(dcd[, sel_list$bb$xyz]),
         ca  = rmsf(dcd[, sel_list$ca$xyz])
    )
}

# --------------------------------------------------------------
cat("using", n_cores, "cores\n")
nothing <- mclapply(1:n_all, run, mc.cores = n_cores)

# for (i in 1:n_all) run(i)

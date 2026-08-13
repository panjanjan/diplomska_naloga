#!/bin/Rscript
# Kljub statistično značilne razlike med RMSF domen, ne pomeni, da
# v njih prihaja do meddomenskega gibanja. Na primer ena domena je
# skoraj toga, druga precej fleksibilna in ni meddomenskega gibanja.
# Ker je razlika v RMSF visoka, bo test statistično značilen.

# Ta skripta je namenjena razreševanju pristranskosti notranje (ne)fleksibilnosti.
# Po zgledu <članek iz late 1990s, ne spomnim se točno> primerja notranjo
# in zunanjo fleksibilnost domen. Notranja fleksibilnost v tem primeru
# pomeni fleksibilnost domene "neodvisno" (ali je res?) od druge domene,
# analogno zunanja fleksibilnost (prisotnost druge domene). Če notranja
# prevlada nad zunanjo, je opisan zgornji primer in ni meddomenskega
# gibanja. Če je notranja dovolj manjša od zunanje, lahko govorimo o
# meddomenskem gibanju.
# ------------------------------------------------------------------------------
library(bio3d)
library(parallel)
library(magrittr)

setwd(Sys.getenv("ROOT"))
source("./scripts/utils.r")

results_target    <- "rmsf_ratios_results.csv"
replicates_target <- "rmsf_ratios_replicates.txt"
proteins_target   <- "rmsf_ratios_proteins.txt"

# meja za p-vrednosti
cutoff <- 0.05

# število jeder za paralelizacijo
n_cores <- min(detectCores() - 1, 10)

data  <- load_data(c("pdb", "traj", "domains"))
n_all <- nrow(data$domains)

# ------------------------------------------------------------------------------
run <- function(i) {
    protein <- data$domains$protein[i]

    dcdfiles <- grep(protein, data$traj, value = TRUE)
    pdbfile  <- grep(protein, data$pdb, value = TRUE)
    assertthat::are_equal(length(dcdfiles), 3)
    assertthat::are_equal(length(pdbfile), 1)

    pdb <- read.pdb(pdbfile, verbose = FALSE)

    domain_bounds <- data$domains[i, -1] |> unlist()
    inds <- atom.select(pdb, "noh", resno = domain_bounds[1]:domain_bounds[2])

    # rmsf-ji
    rmsf1 <- run_replicate(dcdfiles[1], pdb, inds)
    rmsf2 <- run_replicate(dcdfiles[2], pdb, inds)
    rmsf3 <- run_replicate(dcdfiles[3], pdb, inds)

    rm(pdb)

    # izvedi testa
    test1 <- run_test(rmsf1, domain_bounds)
    test2 <- run_test(rmsf2, domain_bounds)
    test3 <- run_test(rmsf3, domain_bounds)

    protein_names <- sub(".dcd", "", basename(dcdfiles))

    # združi rezultate
    df <- rbind(test1, test2, test3)
    df <- cbind(protein_names, df)

    df
}

# gg
run_replicate <- function(dcdfile, pdb, inds) {
}

# izvede statistični test nad rmsfji domen
run_test <- function(rmsf, domain_bounds) {
    rmsf_a <- rmsf[domain_bounds[1]:domain_bounds[2]]
    rmsf_b <- rmsf[domain_bounds[3]:domain_bounds[4]]

    ks_res <- ks.test(rmsf_a, rmsf_b, alternative = "two.sided")

    # logaritmiraj za t-test
    rmsf_a <- log(rmsf_a)
    rmsf_b <- log(rmsf_b)

    t_res <- t.test(rmsf_a, rmsf_b, alternative = "two.sided")

    # enovrstični dataframe
    data.frame(
        ks_stat   = ks_res$statistic,
        ks_pval   = ks_res$p.value,
        ks_pass   = ks_res$p.value < cutoff,
        t_stat    = t_res$statistic,
        t_pval    = t_res$p.value,
        t_pass    = t_res$p.value < cutoff,
        both_pass = (t_res$p.value < cutoff) & (ks_res$p.value < cutoff)
    )
}

# ------------------------------------------------------------------------------
cat("using", n_cores, "cores\n")
cat("using", cutoff, "as cutoff for p-values\n")

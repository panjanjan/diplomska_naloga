#!/bin/Rscript
# WARN: poženi, da vidiš če dela pravilno

# Nad proteini izvede statistična testa t in ks.
#
# Iz trajektorij izračuna RMSF vrednosti in jih razdeli na dva dela glede na domeno.
# Trajektorije so bile poravnane na prvi frame pred depozicijo v ATLAS bazo.
#
# S statističnim testom primerja ali se RMSF vrednosti med domenama razlikujejo.
# Ob statistični značilnosti "lahko" pričakujemo meddomensko gibanje.

# TODO: 3 selekcije: vse brez H, backbone, C-alpha pri izračunu RMSFjev. Glej inds, run_replicate

# ------------------------------------------------------------------------------
library(bio3d)
library(parallel)
library(magrittr)

source(here::here("scripts", "utils.r")

# Rezultate testov shrani v 'results_target'.
# Replikate proteinov, ki imajo statistično značilne razlike shrani v 'replicates_target'.
# Imena proteinov, katerih vsi trije replikati passajo, shrani v 'proteins_target'
out <- list(
    results    = here::here("outputs", "rmsf_test_results.csv")
    replicates = here::here("outputs", "rmsf_test_replicates.txt")
    proteins   = here::here("outputs", "rmsf_test_proteins.txt")
)

# meja za p-vrednosti
cutoff <- 0.05

# število jeder za paralelizacijo
n_cores <- min(detectCores() - 1, 10)

# PDB datotetke, DCD trajektorije, podatki o domenah
data  <- list(
    pdb = list.files(paths$pdb, pattern = ".pdb"),
    traj = list.files(paths$traj, pattern = ".dcd"),
    domains = read.csv(paths$domains)
)

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

# vrne vektor rmsf vrednosti za replikat
run_replicate <- function(dcdfile, pdb, inds) {
    cat(dcdfile, "\n")
    dcd <- read.dcd(dcdfile, verbose = FALSE)
    rmsf(dcd)
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

# --------------------------------------------------------------------
cat("using", n_cores, "cores\n")
cat("using", cutoff, "as cutoff for p-values\n")

# rezultati testov
results <- mclapply(1:n_all, run, mc.cores = n_cores)
results <- do.call(rbind, results)

# imena replikatov proteinov ki passajo OBA testa
passed_replicates <- results[which(results$both_pass), "protein"]

# imena proteinov ki passajo z vsakim replikatom
passed_proteins <- passed_replicates %>%
    sub("_R.", "", .) %>%
    table() %>%
    grep(3, ., value = TRUE) %>%
    names(.)

write.csv(results,            out$results, quote = FALSE, row.names = FALSE)
writeLines(passed_replicates, out$replicates)
writeLines(passed_proteins,   out$proteins)

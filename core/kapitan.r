#!/bin/Rscript
library(bio3d)
library(magrittr)
setwd(Sys.getenv("ROOT"))

# CSV datoteke z izračuni bodo shranjene v "core/*", kjer "*" označuje poddirektorij.
# Glej spodnji seznam. Datoteke so urejene so po proteinih. Glede na pripono datoteke je
# razvidno kateri izračuni so shranjeni notri. Vsaka vsebuje stolpce R1-R3, ki predstavljajo
# izračune po replikatu.
#
# - core/rmsf/*_rmsf.csv      : RMSF vrednosti
# - core/distances/*_dist.csv : razdalje med masnimi centri domen
# - core/psd/*_psd.csv        : "power spectral density" razdalj masnih centrov
# - core/autocorr/*_ac.csv    : rezultati avtokorelacije razdalj masnih centrov
#
data_root <- "core"

paths <- list(
    rmsf      <- file.path(data_root, "rmsf"),
    distances <- file.path(data_root, "distances"),
    psd       <- file.path(data_root, "psd"),
    autocorr  <- file.path(data_root, "autocorr")
)

for (path in paths) {
    if (!dir.exists(path)) dir.create(path)
}

csv_suffix <- list(
    rmsf     = "rmsf.csv",
    dist     = "dist.csv",
    psd      = "psd.csv",
    autocorr = "ac.csv"
)

# wrapper
csv_writer <- \(out_data, out_path) write.csv(out_data, out_path, row.names = FALSE, quote = FALSE)

# ------ GLEDE RMSF -------------------------------------------------------------------------------
#
# Rezultati statističnih testov na RMSF bodo shranjeni v "core/rmsf". Datoteke so
# razdeljene glede na selekcijo v "atom_selections" (glej for-loop dol), kar označujejo
# pripone datotek (all, "bb" kot "backbone", "ca" kot "C-alpha").
#
# - core/rmsf/stat_test_*.csv  : rezultati statističnih testov
# - core/rmsf/replicates_*.csv : replikati, ki passajo teste
# - core/rmsf/proteins_*.csv   : proteini, katerih vsi replikati passajo teste
#
stat_file_names <- list(
    all = list(
        results    = file.path(paths$rmsf, "stat_test_all.csv")
        replicates = file.path(paths$rmsf, "replicates_all.txt")
        proteins   = file.path(paths$rmsf, "proteins_all.txt")
    ),
    backbone = list(
        results    = file.path(paths$rmsf, "stat_test_bb.csv")
        replicates = file.path(paths$rmsf, "replicates_bb.txt")
        proteins   = file.path(paths$rmsf, "proteins_bb.txt")
    ),
    calpha = list(
        results    = file.path(paths$rmsf, "stat_test_ca.csv")
        replicates = file.path(paths$rmsf, "replicates_ca.txt")
        proteins   = file.path(paths$rmsf, "proteins_ca.txt")
    ),
)

# ------ GLEDE MASNIH CENTROV ---------------------------------------------------------------------
# TODO: ...

# ------ GLEDE POWER SPECTRAL DENSITY -------------------------------------------------------------
# TODO: ...

# ------ GLEDE AVTOKORELACIJ ----------------------------------------------------------------------
# TODO: ...

# ------ GLEDE VHODNIH PODATKOV -------------------------------------------------------------------
#
# Inputs za izračune so:
#
# - PDB datoteke proteinov
# - trajektorije v DCD formatu
# - domene, ki jih je določil SWORD2
#
# Seznam proteinov, ki so bili izbrani glede na kakovost rezultatov SWORD2 so
# shranjeni v "two_domains.csv". Vsak protein ima 3 trajektorije (3 replikati),
# 1 PDB datoteko in 1 vrstico v two_domains.
#
data <- list(
    pdb = list.files("atlas_db/PDB", ".pdb", full.names = TRUE),
    dcd = list.files("atlas_db/TRAJ", ".dcd", full.names = TRUE),
    domains = read.csv("two_domains.csv")
)

for (protein in data$domains$protein) {
    # najdi datoteke za protein
    pdbfile <- grep(protein, data$pdb, value = TRUE)
    dcdfiles <- grep(protein, data$traj, value = TRUE)
    assertthat::are_equal(length(pdbfile), 1)
    assertthat::are_equal(length(dcdfiles), 3)

    pdb <- read.pdb(pdbfile, verbose = FALSE)

    # najdi meje domen
    domain_bounds <- data$domains[i, -1] |> unlist()

    # ustvari 3 selekcije
    make_sel <- \(pattern, start, end)

    atom_selections <- list(
        all      = atom.select(pdb, "noh"),
        backbone = atom.select(pdb, "backbone"),
        calpha   = atom.select(pdb, "calpha")
    )

    # TODO: poženi izračune
    run_rmsf(protein, domain_bounds, dcdfiles, atom_selections)
    run_distances()
    run_fft()
    run_autocorr()
}

run_rmsf <- function(domain_bounds, dcdfiles, atom_selections) {
    # meja za p-vrednosti
    cutoff <- 0.05

    # TODO: za vsak selection moram naredit nekaj
    # ...
    for (selection in names(atom_selections)) {
        sel <- atom_selections[[selection]]

        # rmsf glede na izbiro atomov
        rmsf1 <- rmsf(dcdfiles[1][, sel$xyz])
        rmsf2 <- rmsf(dcdfiles[2][, sel$xyz])
        rmsf3 <- rmsf(dcdfiles[3][, sel$xyz])

        # naredi CSV z RMSFji
        out_name <- paste(protein, selection, prefix$rmsf, sep = "_")
        csv_writer(
            out_data = data.frame(R1 = rmsf1, R2 = rmsf2, R3 = rmsf3),
            out_path = file.path(data_target, out_name)
        )

        # statistična analiza
        run_stat_tests()
        test1 <- calc_stat_tests(rmsf1)
        test2 <- calc_stat_tests(rmsf2)
        test3 <- calc_stat_tests(rmsf3)
    }
}

run_stat_tests()

# izvedi t-test in ks-test
stat_tests <- \(rmsf) {
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

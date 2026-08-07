#!/bin/Rscript
library(bio3d)
setwd(Sys.getenv("ROOT"))

source("scripts/utils.r")
source("scripts/core/rmsf.r")
source("scripts/core/distances.r")
source("scripts/core/fft.r")
source("scripts/core/autocorr.r")

data <- load_data(c("pdb", "traj")) # glorified 'list.files'
domains <- load_data("domains") # prebere csv z domenami

for (protein in domains$protein) {
    pdbfile <- grep(protein, data$pdb, value = TRUE)
    dcdfiles <- grep(protein, data$traj, value = TRUE)
    assertthat::are_equal(length(pdbfile), 1)
    assertthat::are_equal(length(dcdfiles), 3)

    pdb <- read.pdb(pdbfile, verbose = FALSE)
    domain_bounds <- data$domains[i, -1] |> unlist()

    # ustvari 3 selekcije
    make_sel <- \(pattern, start, end) atom.select(pdb, pattern, resno = domain_bounds[start]:domain_bounds[end])

    atom_selections <- list(
        all      = list(a = make_sel("noh", 1, 2),      b = make_sel("noh", 3, 4)),
        backbone = list(a = make_sel("backbone", 1, 2), b = make_sel("backbone", 3, 4)),
        calpha   = list(a = make_sel("calpha", 1, 2),   b = make_sel("calpha", 3, 4)),
    )

    # NOTE: mogoče mase za distances, ipd. za ostale, ni nujno da vsak
    # potrebuje iste argumente

    # poženi enega za drugim
    run_rmsf(pdb, atom_selections)
    run_distances(pdb, atom_selections)
    run_fft(pdb, atom_selections)
    run_autocorr(pdb, atom_selections)
}

#!/bin/Rscript
# Proteinsko domeno lahko predstavimo kot masni center atomov, ki jo sestavljajo.
# Meddomensko gibanje lahko opazimo ob relativno velikih premikih teh masnih
# centrov skozi celotno trajektorijo.
# ------------------------------------------------------------------------------
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

# ------------------------------------------------------------------------------
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

target <- here::here("atlas_db", "COM")
if (!dir.exists(target)) dir.create(target, recursive = TRUE)

# ------------------------------------------------------------------------------
# * vrne indekse aminokislin, ki sestavljajo domeno glede na start in end
#   stolpca
# TODO: implement this
dom_inds <- function(dom_df) {
  inds <- c()
  for
}

# * ustvari CSVje z imeni {protein}_dist.csv
# * stolpci R1-R3
# * vsaka vrstica vsebuje razdaljo med masnima centroma domen trenutnega framea
#
# frame R1 R2 R3
# 1     x  x  x
# 2     x  x  x
# 3     x  x  x
# ...
run <- function(i) {
    protein <- data$domains$protein[i]

    dcdfiles <- grep(protein, data$traj, value = TRUE)
    stopifnot(length(dcdfiles) == 3)

    pdbfile <- grep(protein, data$pdb, value = TRUE)
    stopifnot(length(pdbfile) == 1)

    pdb <- read.pdb(pdbfile, verbose = FALSE)

    # BUG: posodobi logiko ker je lahko več vrstic na domeno, dom_inds()
    dom1_inds <- data$domains[i, -1] |> unlist()
    # BUG: posodobi logiko ker je lahko več vrstic na domeno

    # določi kje sta domeni, ignoriraj vodike pri selekciji
    inds_a <- atom.select(pdb, "noh", resno = domain_bounds[1]:domain_bounds[2])
    inds_b <- atom.select(pdb, "noh", resno = domain_bounds[3]:domain_bounds[4])

    # najde mase atomov za izračun masnega centra
    mass_a <- atom2mass(pdb$atom[inds_a$atom, "elety"])
    mass_b <- atom2mass(pdb$atom[inds_b$atom, "elety"])

    # razdalje za vsak replikat
    r1 <- run_replicate(dcdfiles[1], pdb, inds_a, inds_b, mass_a, mass_b)
    r2 <- run_replicate(dcdfiles[2], pdb, inds_a, inds_b, mass_a, mass_b)
    r3 <- run_replicate(dcdfiles[3], pdb, inds_a, inds_b, mass_a, mass_b)

    n_frames <- max(nrow(r1), nrow(r2), nrow(r3))

    # nimajo vsi enako število frame-ov for some reason
    # na koncu skopira zadnjo vrstico da se zapolni do željene velikosti
    pad_replicate <- function(replicate, target_len) {
        cur_len <- nrow(replicate)
        if (cur_len < target_len) {
            # če je začetna dolžina 3, željena pa 5 bo idx = [1,2,3,3,3]
            idx <- c(1:cur_len, rep(cur_len, target_len - cur_len))
            replicate <- replicate[idx, ]
            replicate$frame <- 1:target_len
        }
        replicate
    }

    r1 <- pad_replicate(r1, n_frames)
    r2 <- pad_replicate(r2, n_frames)
    r3 <- pad_replicate(r3, n_frames)

    df <- data.frame(
        "frame" = 1:n_frames,
        "R1" = r1[ns],
        "R2" = r2[ns],
        "R3" = r3[ns]
    )

    out <- paste0(protein, "_dist.csv")
    csv <- file.path(target, out)
    write.csv(df, csv, quote = FALSE, row.names = FALSE)
}

# * izračuna razdalje med masnimi centri domen
# * vrne vektor števil (razdalj)
run_replicate <- function(dcdfile, pdb, inds_a, inds_b, mass_a, mass_b) {
    cat(dcdfile, "\n")
    dcd <- read.dcd(dcdfile, verbose = FALSE)

    # poravnava na prvo domeno
    aligned <- fit.xyz(
        fixed = pdb$xyz,
        mobile = dcd,
        fixed.inds = inds_a$xyz,
        mobile.inds = inds_a$xyz
    )

    # razdeli koordinate v trajektoriji glede na domene
    coords_a <- aligned[, inds_a$xyz]
    coords_b <- aligned[, inds_b$xyz]

    # preko koordinat in mas izračuna masne centre za vsak frame
    com_a <- com.xyz(coords_a, mass = mass_a)
    com_b <- com.xyz(coords_b, mass = mass_b)

    # com sta matrike oblike n×3 (x,y,z)
    # vrne evklidske razdalje med koordinatami
    (com_a - com_b)**2 |>
        rowSums() |>
        sqrt()
}

# ------------------------------------------------------------------------------
cat("using", n_cores, "cores\n")
invisible(mclapply(1:n_all, run, mc.cores = n_cores))

#!/bin/Rscript
# Proteinsko domeno lahko predstavimo kot masni center atomov, ki jo sestavljajo.
# Meddomensko gibanje lahko opazimo ob relativno velikih premikih teh masnih
# centrov skozi celotno trajektorijo.
#
# Selekcijo atomov sprejme kot parameter: Rscript calc_dist.r "<sel>"
# ------------------------------------------------------------------------------
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

sels <- c("noh", "backbone", "calpha")

sel <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(sel %in% sels)

# ------------------------------------------------------------------------------
# število jeder za paralelizacijo
n_cores <- min(detectCores() - 1, 8)

# PDB datotetke, DCD trajektorije proteini
# dist da preskoči obstoječe datoteke
data  <- list(
    pdb     = list.files(paths$pdb,  pattern = ".pdb", full.names = TRUE),
    traj    = list.files(paths$traj, pattern = ".dcd", full.names = TRUE),
    domains = read.csv(paths$domains),
    dist    = list.files(paths$dist)
)

proteins <- unique(data$domains$protein)

target <- here::here("outputs", "COM")
if (!dir.exists(target)) dir.create(target, recursive = TRUE)

# ------------------------------------------------------------------------------
# * vrne indekse aminokislin, ki sestavljajo i-to domeno
dom_inds <- function(protein, i) {
  d <- data$domains[data$domains == protein, ]
  dom <- d[d$domain == i, ]
  bounds <- list(start = dom$start, end = dom$end)
  lapply(1:nrow(dom), \(i) bounds$start[i]:bounds$end[i]) |>
    unlist()
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
run <- function(protein) {
  # preskoči če že obstajajo vsi trije CSVji tega proteina
  expected <- paste(protein, sels, "dist.csv", sep = "_")
  if (all(expected %in% data$dist)) {
    cat(protein, "done\n")
    return(invisible(NULL))
  }

  dcdfiles <- grep(protein, data$traj, value = TRUE)
  stopifnot(length(dcdfiles) == 3)

  pdbfile <- grep(protein, data$pdb, value = TRUE)
  stopifnot(length(pdbfile) == 1)

  pdb <- read.pdb(pdbfile, verbose = FALSE)

  # indeksi od obeh domen
  domain_bounds <- list(
    dom1 = dom_inds(protein, 1),
    dom2 = dom_inds(protein, 2)
  )

  # določi kje sta domeni
  inds_d1 <- atom.select(pdb, sel, resno = domain_bounds$dom1)
  inds_d2 <- atom.select(pdb, sel, resno = domain_bounds$dom2)

  # najde mase atomov za izračun masnega centra
  mass_d1 <- atom2mass(pdb$atom[inds_d1$atom, "elety"])
  mass_d2 <- atom2mass(pdb$atom[inds_d2$atom, "elety"])

  # razdalje za vsak replikat
  #: ~1001x vec
  r1 <- run_replicate(dcdfiles[1], pdb, inds_d1, inds_d2, mass_d1, mass_d2)
  r2 <- run_replicate(dcdfiles[2], pdb, inds_d1, inds_d2, mass_d1, mass_d2)
  r3 <- run_replicate(dcdfiles[3], pdb, inds_d1, inds_d2, mass_d1, mass_d2)

  rm(pdb)

  n_frames <- max(length(r1), length(r2), length(r3))

  # nimajo vsi enako število frame-ov for some reason
  # na koncu skopira zadnjo vrstico da se zapolni do željene velikosti
  pad_replicate <- function(replicate, target_len) {
    cur_len <- length(replicate)
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
    "R1" = r1,
    "R2" = r2,
    "R3" = r3
  )

  out <- paste(protein, sel, "dist.csv", sep = "_")
  csv <- file.path(target, out)
  write.csv(df, csv, quote = FALSE, row.names = FALSE)
}

# * izračuna razdalje med masnimi centri domen
# * vrne vektor števil (razdalj)
run_replicate <- function(dcdfile, pdb, inds_d1, inds_d2, mass_d1, mass_d2) {
  cat(dcdfile, "\n")
  dcd <- read.dcd(dcdfile, verbose = FALSE)

  # poravnava na prvo domeno
  aligned <- fit.xyz(
    fixed = pdb$xyz,
    mobile = dcd,
    fixed.inds = inds_d1$xyz,
    mobile.inds = inds_d1$xyz
  )

  rm(dcd)

  # razdeli koordinate v trajektoriji glede na domene
  coords_d1 <- aligned[, inds_d1$xyz]
  coords_d2 <- aligned[, inds_d2$xyz]

  # preko koordinat in mas izračuna masne centre za vsak frame
  com_d1 <- com.xyz(coords_d1, mass = mass_d1)
  com_d2 <- com.xyz(coords_d2, mass = mass_d2)

  # com sta matrike oblike n×3 (x,y,z)
  # vrne evklidske razdalje med koordinatami
  (com_d1 - com_d2)**2 |>
    rowSums() |>
    sqrt()
}

# ------------------------------------------------------------------------------
cat("using", n_cores, "cores\n")

mclapply(proteins, run, mc.cores = n_cores)

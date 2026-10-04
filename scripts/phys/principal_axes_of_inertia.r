#!/bin/Rscript

# izračuna kote med vztrajnostnimi osemi domen
# ------------------------------------------------------------------------------
library(bio3d)
library(parallel)

source(here::here("scripts", "utils.r"))

n_cores <- min(detectCores() - 1, 10)

data <- list(
  pdb = list.files(paths$pdb, ".pdb", full.names = TRUE),
  traj = list.files(paths$traj, ".dcd", full.names = TRUE),
  domains = read.csv(paths$domains)
)
n_all <- nrow(data$domains)

target <- file.path("outputs", "PAI")
if (!dir.exists(target)) dir.create(target)

# ------------------------------------------------------------------------------
# * bio3d objekti imajo koordinate v obliki [x1,y2,z1,x2,...], kar je nadležno
# * pretvori koordinate v n×3 matriko, n×(x,y,z), tako kot vrne `bio3d::com.xyz`
matrix_coords <- function(coords) {
  n <- length(coords)
  x_inds <- seq(1, n, 3)
  y_inds <- seq(2, n, 3)
  z_inds <- seq(3, n, 3)
  x <- coords[x_inds]
  y <- coords[y_inds]
  z <- coords[z_inds]
  matrix(
    c(x, y, z),
    ncol = 3,
    byrow = FALSE,
    dimnames = list(NULL, c("x", "y", "z"))
  )
}

# * izračuna inertia tensor
# * vrne 3x3 matriko
inertia_tensor <- function(coords, masses) {
  x <- coords[, "x"]
  y <- coords[, "y"]
  z <- coords[, "z"]

  i_xx <- sum(masses * (y^2 + z^2))
  i_yy <- sum(masses * (x^2 + z^2))
  i_zz <- sum(masses * (x^2 + y^2))
  i_xy <- -sum(masses * x * y)
  i_xz <- -sum(masses * x * z)
  i_yz <- -sum(masses * y * z)

  matrix(
    c(
      i_xx, i_xy, i_xz,
      i_xy, i_yy, i_yz,
      i_xz, i_yz, i_zz
    ),
    nrow = 3, byrow = TRUE
  )
}

# * csv z imenom {protein}_angles.csv
# * stolpci R1-R3
# * vsaka vrstica vsebuje kot med prvo vztrajnostno osjo ene domene in prvo,
# drugo, tretjo vztrajnostjo osjo druge domene
#   domen v trenutnem frame-u
#
# frame R1 R2 R3
# 1     x  x  x
# 2     x  x  x
# 3     x  x  x
# ...
run <- function(i) {
  protein <- data$domains$protein[i]

  pdbfile <- grep(protein, data$pdb, value = TRUE)
  dcdfiles <- grep(protein, data$traj, value = TRUE)
  stopifnot(length(pdbfile) == 1)
  stopifnot(length(dcdfiles) == 3)

  # najdi meje domen
  domain_bounds <- data$domains[i, -1] |> unlist()

  # izberi domeni
  domain_bounds <- list(
    dom1 = dom_inds(protein, 1),
    dom2 = dom_inds(protein, 2)
  )
  pdb <- read.pdb(pdbfile, verbose = FALSE)
  inds_d1 <- atom.select(pdb, "noh", resno = domain_bounds$dom1)
  inds_d2 <- atom.select(pdb, "noh", resno = domain_bounds$dom2)

  # najde mase atomov za izračun masnega centra
  mass_d1 <- atom2mass(pdb$atom[inds_d1$atom, "elety"])
  mass_d2 <- atom2mass(pdb$atom[inds_d2$atom, "elety"])

  # določimo kote za vsak frame za vsak replikat
  # v enem prehodu izračunamo vse tri glavne osi
  r1 <- run_replicate(dcdfiles[1], pdb, inds_d1, inds_d2, mass_d1, mass_d2)
  r2 <- run_replicate(dcdfiles[2], pdb, inds_d1, inds_d2, mass_d1, mass_d2)
  r3 <- run_replicate(dcdfiles[3], pdb, inds_d1, inds_d2, mass_d1, mass_d2)

  n_frames <- max(nrow(r1), nrow(r2), nrow(r3))

  # nimajo vsi enako število frame-ov
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

  angles <- data.frame(
    frame = 1:n_frames,
    R1_p1 = r1$P1,
    R1_p2 = r1$P2,
    R1_p3 = r1$P3,
    R2_p1 = r2$P1,
    R2_p2 = r2$P2,
    R2_p3 = r2$P3,
    R3_p1 = r3$P1,
    R3_p2 = r3$P2,
    R3_p3 = r3$P3
  )

  out <- paste0(protein, "_angles.csv")
  csv <- file.path(target, out)
  write.csv(angles, csv, quote = FALSE, row.names = FALSE)
}

# * izračuna kot med vsemi tremi vztrajnostnimi osmi domen za vsak frame
#   v trajektoriji
# * vrne data.frame z stolpci R1, R2, R3
run_replicate <- function(dcdfile, pdb, inds_d1, inds_d2, mass_d1, mass_d2) {
  cat(dcdfile, "\n")
  dcd <- read.dcd(dcdfile, verbose = FALSE)
  n_frames <- nrow(dcd)

  # razdeli koordinate v trajektoriji glede na domene
  coords_d1 <- dcd[, inds_d1$xyz]
  coords_d2 <- dcd[, inds_d2$xyz]

  # preko koordinat in mas izračuna masne centre za vsak frame
  com_d1 <- com.xyz(coords_d1, mass = mass_d1)
  com_d2 <- com.xyz(coords_d2, mass = mass_d2)

  # izračuna kot med vsakimi od treh osi v frame-u
  # vektor dolžine 3
  angles <- lapply(1:n_frames, \(i) {
    crds_d1 <- matrix_coords(coords_d1[i, ])
    crds_d2 <- matrix_coords(coords_d2[i, ])

    # centrira glede na masni center
    centered_d1 <- scale(crds_d1, center = com_d1[i, ], scale = FALSE)
    centered_d2 <- scale(crds_d2, center = com_d2[i, ], scale = FALSE)

    inertia_d1 <- inertia_tensor(centered_d1, mass_d1)
    inertia_d2 <- inertia_tensor(centered_d2, mass_d2)

    # lastni vektorji
    axes_d1 <- eigen(inertia_d1)$vectors
    axes_d2 <- eigen(inertia_d2)$vectors

    vapply(1:3, \(axis) {
      paxis_d1 <- axes_d1[, axis]
      paxis_d2 <- axes_d2[, axis]

      sum(paxis_d1 * paxis_d2) |>
        abs() |>
        acos()
    }, numeric(1))
  })

  # vse elemente seznama združi v matriko
  angles_mat <- do.call(rbind, angles)

  cat("done\n")

  data.frame(
    frame = seq_len(n_frames),
    P1 = angles_mat[, 1],
    P2 = angles_mat[, 2],
    P3 = angles_mat[, 3]
  )
}

# ------------------------------------------------------------------------------
cat("using", n_cores, "cores\n")
mclapply(1:n_all, run, mc.cores = n_cores)

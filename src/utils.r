paths <- list(
    traj        = here::here("atlas_db", "trajectories"),
    pdb_old     = here::here("atlas_db", "PDB"),
    pdb         = here::here("atlas_db", "PDB_chained"),
    sword       = here::here("outputs", "sword_results.csv"),
    sword_clean = here::here("outputs", "sword_results_clean.csv"),
    domains     = here::here("outputs", "two_domains.csv"),
    t_test_res  = here::here("outputs", "rmsf_ttest_results.csv"),
    ks_test_res = here::here("outputs", "rmsf_kstest_results.csv"),
    dist        = here::here("outputs", "COM"),
    angles      = here::here("outputs", "PAI"),
    atlas_rmsf  = here::here("atlas_db", "RMSF"),
    rmsf        = here::here("outputs", "RMSF"),
    psd         = here::here("outputs", "PSD"),
    acf         = here::here("outputs", "ACF")
)

# ustvari inverted window za plottanje
# no more flashbangs
dark_plot <- function() {
    par(
        bg = "black", # Barva ozadja
        fg = "white", # Osnovna barva (okvirji in črte)
        col.axis = "white", # Oznake na oseh (številke)
        col.lab = "white", # Imeta osi
        col.main = "white" # Glavni naslov
    )
    plot(1)
}

light_plot <- function() {
    par(
        bg = "white", # Barva ozadja
        fg = "black", # Osnovna barva (okvirji in črte)
        col.axis = "black", # Oznake na oseh (številke)
        col.lab = "black", # Imeta osi
        col.main = "black" # Glavni naslov
    )
    plot(1)
}

mplt <- function(s, title) {
  or <- "orange"
  plot(
    s,
    type = "l",
    col = "gray",
    plot.first = grid(),
    ...
  )
  points(
    s,
    pch = 19,
    col = adjustcolor(or, alpha.f = 0.6)
  )
  abline(h = mean(s), col = "red", lwd = 1)
}

# ustvari PDB datoteko iz trajektorije proteina
# traj: abs path do DCD trajektorije
# pdb: abs path do PDB datoteke
pdb_from_trajectory <- function(traj, pdb) {
    cat("this may take a while\n")
    t <- bio3d::read.dcd(traj, verbose = FALSE)
    p <- bio3d::read.pdb(pdb, verbose = FALSE)
    o <- here::here(sub("dcd", "pdb", basename(dcdfile)))
    bio3d::write.pdb(
        pdb = p,
        xyz = t,
        file = o
    )
    for (i in 1:10) cat("ne pozabi dodat na .gitignore!!!!!!!!!!!!!!!!!!\n")
}

load_protein <- function(protein, replicate) {
    pdb <- paths$pdb |>
        list.files(pattern = protein, full.names = TRUE) |>
        bio3d::read.pdb(verbose = FALSE)
    dcd <- paths$traj |>
        list.files(
            pattern = paste0("^", protein, "_R", replicate, "\\.dcd$"),
            full.names = TRUE
        ) |>
        bio3d::read.dcd(verbose = FALSE)
    list(pdb = pdb, dcd = dcd)
}

# * vrne indekse aminokislin, ki sestavljajo i-to domeno
dom_inds <- function(prot, i) {
  library(magrittr)
  read.csv(paths$domains) %>%
    .[.[]$protein == prot, ] %>%
    .[.[]$domain == i, ] %>%
    {
      d <- .[]
      lapply(
        1:nrow(d),
        \(i) d[i, "start"]:d[i, "end"]
      ) |>
        unlist()
    }
}

plotseg <- function(d, ...) {
  # Define data and group assignments for each point
  n <- length(d)
  x <- 1:n

  pal <- rainbow(n)
  pal_mod <- sapply(pal, adjustcolor, alpha.f = 0.2)

  # Create empty plot
  plot(x, d, type = "n", ...)

  # Loop through each segment (from point i to i+1)
  for (i in 1:(n - 1)) {
    # Color segment based on the group of the starting point
    segments(x[i], d[i], x[i + 1], d[i + 1], col = pal_mod[i], lwd = 2)
  }

  points(x, d, col = pal, pch = 19)
  abline(h = mean(d), col = "black", lwd = 1.5, lty = 2)
}

paths <- list(
    traj        = here::here("atlas_db", "trajectories"),
    pdb_old     = here::here("atlas_db", "PDB"),
    pdb         = here::here("atlas_db", "PDB_chained"),
    sword       = here::here("outputs", "sword_results.csv"),
    sword_clean = here::here("outputs", "sword_results_clean.csv"),
    domains     = here::here("outputs", "two_domains.csv"),
    dist        = here::here("outputs", "COM"),
    angles      = here::here("outputs", "PAI"),
    rmsf        = here::here("outputs", "RMSF")
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

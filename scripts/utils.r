root <- Sys.getenv("ROOT")

paths <- list(
    traj    = file.path(root, "atlas_db/trajectories"),
    pdb     = file.path(root, "atlas_db/PDB"),
    domains = file.path(root, "two_domains.csv"),
    dist    = file.path(root, "atlas_db/COM"),
    angles  = file.path(root, "atlas_db/PAI")
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

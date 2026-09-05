paths <- list(
    domains     = here::here("two_domains.csv"),
    traj        = here::here("atlas_db", "trajectories"),
    pdb         = here::here("atlas_db", "PDB"),
    sword       = here::here("outputs", "sword_results.csv"),
    sword_clean = here::here("outputs", "sword_results_clean.csv"),
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

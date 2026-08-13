#!/bin/Rscript
# NOTE: not used..
# NOTE: not used..
# NOTE: not used..

# seznam proteinov
proteins <- read.csv("two_domains.csv")$protein

# razdalje so v xvg datotekah
xvgs <- list.files(path = "atlas_db/TRAJ", pattern = "xvg", full.names = TRUE)

target <- file.path("atlas_db", "COM")

# csv z imenom {protein}_com_dist_gmx.csv
#
# stolpci R1-R3
# vsaka vrstica vsebuje razdaljo med masnima centroma domen v trenutnem frame-u
#
# frame R1 R2 R3
# 1     x  x  x
# 2     x  x  x
# 3     x  x  x
# ...
run <- function(protein) {
    cat(protein, "...\n")

    # 3 replikati, 3 datoteke na protein
    xvgfiles <- grep(protein, xvgs, value = TRUE)
    assertthat::are_equal(length(xvgfiles), 3)

    r1 <- run_replicate(xvgfiles[1])
    r2 <- run_replicate(xvgfiles[2])
    r3 <- run_replicate(xvgfiles[3])

    min_frames <- min(length(r1), length(r2), length(r3))
    ns <- 1:min_frames

    df <- data.frame(
        "frame" = ns,
        "R1" = r1[ns],
        "R2" = r2[ns],
        "R3" = r3[ns]
    )

    out <- paste0(protein, "_com_dist_gmx.csv")
    csv <- file.path(target, out)
    write.csv(df, csv, quote = FALSE, row.names = FALSE)
}

run_replicate <- function(xvgfile) {
    d <- read.delim(xvgfile)

    # razbij nize na ločene vrednosti
    # vzami samo razdalje
    d[, 1] |>
        trimws() |>
        strsplit("\\s+") |>
        lapply(\(x) {
            as.numeric(x[2])
        }) |>
        unlist()
}

for (p in proteins) run(p)

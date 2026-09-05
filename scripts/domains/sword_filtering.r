#!/bin/Rscript
# WARN: poženi da vidiš če dela pravilno

# Format podatkov ---------------------------------------------------------------------------------------
#
# `sword_results.csv` vsebuje podatke iz JSON datotek, ki jih ustvari SWORD2.
# Vrstice se združujejo po proteinih in proteini po svojih particijah.
#
# protein:   PDB koda ter veriga, ki je bila uporabljena za SWORD2
# aindex:    ambiguity index proteina
# partition: indeks particije. Optimalna ima 0, alternativne 1 ali več
# quality:   ocena particije
# domain:    indeks domene. Prva domena 1, druga domena 2, ...
# AUL:       AUL vrednost domene
# start:     prva aminokislina domene
# end:       zadnja aminokislina domene
#
# Primer za 1a62_A:
#
# protein aindex partition quality domain AUL start end
# 1a62_A  1      0         0       1      81  1     130 <---| opt.
# 1a62_A  1      1         1       1      70  1     47 <----| alt. 1
# 1a62_A  1      1         1       2      0   48    94      |
# 1a62_A  1      1         1       3      46  95    130     |
# 1a62_A  1      2         3       1      72  1     47 <----| alt. 2
# 1a62_A  1      2         3       2      8   48    130     |
# 1a62_A  1      3         1       1      76  1     130 <---| alt. 3
# 1a62_A  1      3         1       2      0   48    94      |
#
# Filtriranje glede na sledeče pogoje
#
# 1. več kot ena domena v particiji
# 2. razmerje med katerokoli domeno v particiji naj ne bo večje od 1:2 (1:3?)
# 3. ambiguity index med 1 in 3
# 4. AUL vrednosti domen vsaj 75 (najmanj 50)
#
# Uporabi samo optimalne particije.
# -------------------------------------------------------------------------------------------------------
library(dplyr)

source(here::here("scripts", "utils.r"))

csv <- read.csv(paths$sword, header = TRUE) |>
    as_tibble() |>
    group_by(protein)

# Število domen -----------------------------------------------------------------------------------------
# Prvi korak je odstraniti vse, ki imajo samo eno domeno.
md_opt <- csv |>
    filter(partition == 0) |>
    filter(max(domain) > 1)

# Razmerja med velikostmi domen--------------------------------------------------------------------------
# Določene domene morajo biti približno enako velike. Uporabil sem razmerji 1 proti 2 (meja 0,5) in 1 proti 3 (meja 0,3).
# `check_ratios` vrne seznam proteinov pri katerih so razmerja med vsako domeno večja od meje.
prot_group <- \(d, p) d[d$protein == p, ]

# razmerje določa cutoff parameter, med 0 in 1
check_ratios <- function(df, cutoff = 0.3) {
    keep <- c()
    for (prot in unique(df$protein)) {
        d2 <- prot_group(df, prot) |> as.data.frame() # tu nočem tibble

        # velikosti domen
        dsizes <- d2$end - d2$start + 1 # ker je 1-based

        # kombinacije vseh domen: (1,2), (1,3), (2,3) ...
        # transponirano, da je dimenzije N×2
        m <- combn(d2$domain, m = 2) |> t()

        # razmerje manjša/večja domena, da je rezultat v intervalu [0,1]
        ratios <- apply(m, 1, \(x) min(dsizes[x]) / max(dsizes[x]))

        # ohrani protein, če so razmerja primerna
        if (all(ratios > cutoff)) {
            keep <- c(keep, prot)
        }
    }
    keep
}

md_keep <- check_ratios(md_opt, cutoff = 0.5)

# A-index -----------------------------------------------------------------------------------------------
# Manjši kot je A-index, boljša je dekompozicija, saj je manj _ambiguous_. Izločil sem vse, ki imajo index manj kot 4.
md_keep <- filter(md_keep, aindex < 4)

# AUL vrednosti -----------------------------------------------------------------------------------------
# Vse domene v particiji morajo imeti dobre AUL vrednosti, na primer >75.
# Koda vzame najslabšo AUL vrednost med domenami v particiji. Če je najslabša vrednost večja
# od meje, potem so vse ostale večje ali enake tej vrednosti. Protein v tem primeru ostane.

# vsebuje eno domeno na protein
bad_doms <- data.frame()

for (prot in unique(md_keep$protein)) {
    df <- prot_group(md_keep, prot)

    # najdi najslabšo domeno
    bad_dom <- which(df$AUL == min(df$AUL))[1]

    # dodaj na seznam
    bad_doms <- rbind(bad_doms, df[bad_dom, ])
}

stopifnot(all(group_keys(bad_doms) == group_keys(md_keep)))  # zajame vse proteine
stopifnot(all(bad_doms$protein == unique(bad_doms$protein))) # ni duplikatov

# proteini, ki ostanejo
sel_prot <- bad_doms[which(bad_doms$AUL > 75), "protein"]
final    <- filter(md_keep, protein %in% sel_prot)

# -------------------------------------------------------------------------------------------------------
write.csv(
    file      = paths$sword_clean,
    x         = final,
    quote     = FALSE,
    row.names = FALSE
)

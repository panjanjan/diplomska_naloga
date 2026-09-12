#!/bin/Rscript
# Filtriranje glede na sledeče pogoje
#
# 1. več kot ena domena v particiji
# 2. razmerje med katerokoli domeno v particiji naj ne bo večje od 1:2
# 3. ambiguity index med 1 in 3
# 4. AUL vrednosti domen večje od 75
#
# Uporabi samo optimalne particije.
# -------------------------------------------------------------------------------------------------------
library(dplyr, warn.conflicts = FALSE)

source(here::here("scripts", "utils.r"))

# -------------------------------------------------------------------------------------------------------
# za izpisovanje sporočil
pad_msg <- \(msg) paste(c(msg, rep(" ", 41 - nchar(msg))), collapse = "")
n_prot  <- \(d) nrow(group_keys(d))

# povleče ven skupino vrstic, ki pripadajo proteinu p v dataframu d
prot_group <- \(d, p) d[d$protein == p, ]

# Data --------------------------------------------------------------------------------------------------
cat(pad_msg("action"), "|", "num prot.\n")
cat(rep("-", 53), "\n", sep = "")
cat(pad_msg("reading sword results"), "| ")

csv <- read.csv(paths$sword, header = TRUE) |>
    as_tibble() |>
    group_by(protein)

cat(n_prot(csv), "\n")

# Število domen -----------------------------------------------------------------------------------------
# Prvi korak je odstraniti vse, ki imajo samo eno domeno.
cat(pad_msg("removing one-domain proteins"), "| ")

md_opt <- csv |>
    filter(partition == 0) |>
    filter(max(domain) > 1)

cat(n_prot(md_opt), "\n")

# Razmerja med velikostmi domen--------------------------------------------------------------------------
# Določene domene morajo biti približno enako velike. Uporabi razmerje 1 proti 2 (meja 0,5).
cat(pad_msg("removing proteins with big domains ratios"), "| ")

# vrne seznam proteinov pri katerih so razmerja med vsako domeno večja od meje.
check_ratios <- function(df, cutoff = 0.5) {
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
md_keep <- filter(md_opt, protein %in% md_keep)

cat(n_prot(md_keep), "\n")

# A-index -----------------------------------------------------------------------------------------------
# Manjši kot je A-index, boljša je dekompozicija, saj je manj _ambiguous_. Izloči vse, ki imajo index izven (1,3).
cat(pad_msg("removing proteins with big A-index"), "| ")

md_keep <- filter(md_keep, aindex > 0 && aindex < 4)

cat(n_prot(md_keep), "\n")

# AUL vrednosti -----------------------------------------------------------------------------------------
# Vse domene v particiji morajo imeti dobre AUL vrednosti, na primer >75.
# Koda vzame najslabšo AUL vrednost med domenami v particiji. Če je najslabša vrednost večja
# od meje, potem so vse ostale večje ali enake tej vrednosti. Protein v tem primeru ostane.
cat(pad_msg("removing proteins with low AUL values"), "| ")

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
sel_prot <- bad_doms[which(bad_doms$AUL > 75), "protein"] |> unlist()
final    <- filter(md_keep, protein %in% sel_prot)

cat(n_prot(final), "\n")

# -------------------------------------------------------------------------------------------------------
cat(rep("-", 53), "\n", sep = "")
cat(pad_msg("saving clean list of proteins"), "|", n_prot(final), "\n")

write.csv(
    file      = paths$sword_clean,
    x         = final,
    quote     = FALSE,
    row.names = FALSE
)

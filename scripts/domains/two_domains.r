#!/bin/Rscript
# Iz prečiščenih SWORD2 rezultatov ohrani proteine z natanko dvema
# domenama.
#
# Format ostane enak kot v sword_results_clean.csv: ena vrstica na
# proteinsko enoto (PU), saj domene niso nujno zvezne. Glej
# build_decompositions_csv.py.
# ------------------------------------------------------------------
library(dplyr, warn.conflicts = FALSE)

source(here::here("scripts", "utils.r"))

# ------------------------------------------------------------------
cat("reading\n")
clean <- read.csv(paths$sword_clean)

cat("trimming\n")
two <- clean |>
  group_by(protein) |>
  filter(n_distinct(domain) == 2) |>
  ungroup()

cat("writing", n_distinct(two$protein), "proteins\n")
write.csv(two, paths$domains, quote = FALSE, row.names = FALSE)

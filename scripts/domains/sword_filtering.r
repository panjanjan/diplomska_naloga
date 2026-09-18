#!/bin/Rscript
# Filtriranje SWORD2 rezultatov glede na sledeče pogoje
#
# 1. optimalna particija (partition == 0)
# 2. več kot ena domena v particiji (two_domains.r nato izbere natanko 2)
# 3. razmerje velikosti manjše in večje domene >= 0,5 (največ 1:2)
# 4. ambiguity index med 1 in 3
# 5. AUL vrednosti vseh domen > 75
#
# Domena je lahko sestavljena iz več proteinskih enot (PU), ki niso
# nujno zvezne (glej build_decompositions_csv.py). Vsaka PU ima svojo
# vrstico, zato je velikost domene vsota dolžin vseh njenih PU-jev.
# Format izhoda je enak formatu vhoda.
# ------------------------------------------------------------------
library(dplyr, warn.conflicts = FALSE)

source(here::here("scripts", "utils.r"))

# ------------------------------------------------------------------
# izpiše sporočilo in število proteinov, d vrne nespremenjen
report <- function(d, msg) {
  cat(formatC(msg, width = -41), "|", n_distinct(d$protein), "\n")
  invisible(d)
}

# ------------------------------------------------------------------
final <- read.csv(paths$sword) |>
  as_tibble() |>
  report("reading sword results") |>
  filter(partition == 0) |>

  # -- [Število domen: odstrani vse, ki imajo samo eno domeno.] --------------
  group_by(protein) |>
  filter(n_distinct(domain) > 1) |>
  report("removing one-domain proteins") |>

  # -- [Razmerja velikosti: domene morajo biti približno enako velike.] ------
  # Če je razmerje najmanjše in največje domene >= meje, potem so vsa razmerja
  # med pari domen >= meje. Npr. domene velikosti 1, 2, 3 aminokisline. Min=1,
  # max=3 => razmerje min-max=1:3. Ostala razmerja, 1:2 in 2:3 so večja od 1:3.
  group_by(protein, partition, domain) |>
  mutate(size = sum(end - start + 1)) |> # 1-based, zato + 1
  group_by(protein) |>
  filter(min(size) / max(size) >= 0.5) |>
  select(-size) |>
  report("removing proteins with big domain ratios") |>

  # -- [A-index: manjši je boljši (manj ambiguous).] ------------------------
  filter(between(aindex, 1, 3)) |>
  report("removing proteins with bad A-index") |>

  group_by(protein, partition) |>
  filter(quality >= 3) |>
  report("removing proteins with low quality") |>

  # -- [AUL: če je najslabša domena nad mejo, so nad mejo vse.] -------------
  group_by(protein) |>
  filter(min(AUL) > 75) |>
  report("removing proteins with low AUL values") |>
  ungroup()

# ------------------------------------------------------------------
cat(strrep("-", 53), "\n", sep = "")
report(final, "saving clean list of proteins")

write.csv(
  file      = paths$sword_clean,
  x         = final,
  quote     = FALSE,
  row.names = FALSE
)

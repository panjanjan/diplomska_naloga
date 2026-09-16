#!/bin/Rscript
# S statističnima testoma t in ks primerja ali se RMSF vrednosti med domenama
# razlikujejo.
#
# Trajektorije so bile poravnane na prvi frame pred depozicijo v ATLAS bazo (I
# think).
# -----------------------------------------------------------------------------
library(bio3d)
library(stringr)
library(parallel)
library(magrittr)

source(here::here("scripts", "utils.r"))

cfg <- list(
  # meja za p-vrednosti
  cutoff = 0.05,

  # število jeder za paralelizacijo
  n_cores = min(detectCores() - 1, 10),

  # TODO: doc
  results_t  = here::here("outputs", "rmsf_ttest_results.csv"),
  results_ks = here::here("outputs", "rmsf_kstest_results.csv")
)

# RMSFji in podatki o domenah
data <- list(
  rmsf    = list.files(paths$rmsf, pattern = ".csv", full.names = TRUE),
  domains = read.csv(paths$domains)
)

n_all <- nrow(data$domains)

# -----------------------------------------------------------------------------
# za vsak protein izvede statistični test nad vsemi selekcijami in replikati.
# Vrne 9-vrstični dataframe, na primer
#
#   protein replicate selection stat pval pass
#   1dd3_A  1         all       0    0    TRUE
#   1dd3_A  1         bb        0    0    TRUE
#   1dd3_A  1         ca        0    0    FALSE
#   1dd3_A  2         all       0    0    TRUE
#   1dd3_A  2         bb        0    0    TRUE
#   1dd3_A  2         ca        0    0    FALSE
#   1dd3_A  3         all       0    0    TRUE
#   1dd3_A  3         bb        0    0    TRUE
#   1dd3_A  3         ca        0    0    FALSE
#
run <- function(i, stat_test) {
  protein <- data$domains$protein[i]

  rmsffiles <- grep(protein, data$rmsf, value = TRUE)
  assertthat::are_equal(length(rmsffiles), 3)

  domain_bounds <- data$domains[i, -1] |> unlist()

  # ---[ združi in vrni rezultate ]---
  #: 3xlist(r1 = 3xdata.frame, r2 = 3xdata.frame, r3 = 3xdata.frame)
  res_l <- list(
    r1 = run_selection(rmsffiles[1], domain_bounds, stat_test),
    r2 = run_selection(rmsffiles[2], domain_bounds, stat_test),
    r3 = run_selection(rmsffiles[3], domain_bounds, stat_test)
  )

  # povleči selekcije iz imen namesto da hardcodeaš vektor, da
  # je vedno pravi vrstni red
  selections <- gsub(".*rmsf_(.*).csv", "\\1", basename(rmsffiles), perl = TRUE)

  d <- data.frame(
    protein = rep(protein, 9),
    replicate = rep(1:3, each = 3),
    selection = rep(selections, 3)
  )

  d <- do.call(rbind, res_l) %>% cbind(d, .)
  d
}

# vrne 3 vrstični dataframe. Vsaka vrstica hrani rezultate statističnega testa
# izvedenega nad RSMFji določenega replikata.
run_selection <- function(rmsffile, domain_bounds, stat_test) {
  csv <- read.csv(rmsffile)

  # razdeli RMSF vrednosti na domene
  #: 2xlist(dom1 = nxvec, dom2 = mxvec)
  part_rmsf <- \(rmsfvals) {
    list(
      dom1 = rmsfvals[domain_bounds[1]:domain_bounds[2]],
      dom2 = rmsfvals[domain_bounds[3]:domain_bounds[4]]
    )
  }

  #: 3xdata.frame(stat, pval, pass)
  rbind(
    run_replicate(part_rmsf(csv$R1), stat_test),
    run_replicate(part_rmsf(csv$R2), stat_test),
    run_replicate(part_rmsf(csv$R3), stat_test)
  )
}

# za t-test logaritmira vrednosti. vrne eno-vrstični dataframe
run_replicate <- function(rmsfpart, stat_test) {
  test_res <- switch(stat_test,
    "t" = {
      t.test(
        log(rmsfpart$dom1),
        log(rmsfpart$dom2),
        alternative = "two.sided"
      )
    },
    "ks" = {
      ks.test(
        rmsfpart$dom1,
        rmsfpart$dom2,
        alternative = "two.sided"
      )
    }
  )

  #: 1xdata.frame(stat, pval, pass)
  data.frame(
    stat = test_res$statistic,
    pval = test_res$p.value,
    pass = (test_res$p.value < cfg$cutoff)
  )
}

# -----------------------------------------------------------------------------
main <- \(stat_test) {
  cat("runing ", stat_test, "-test analysis\n", sep="")

  # vrne list, ki ga zlepi v data.frame/matriko
  results <- mclapply(1:n_all, run, mc.cores = cfg$n_cores, stat_test = stat_test)
  results <- do.call(rbind, results)

  out_name <- paste0("results_", stat_test)
  write.csv(results, cfg[[out_name]], quote = FALSE, row.names = FALSE)
}

cat("using", cfg$n_cores, "cores\n")
cat("using", cfg$cutoff, "as cutoff for p-values\n")

main("t")
main("ks")

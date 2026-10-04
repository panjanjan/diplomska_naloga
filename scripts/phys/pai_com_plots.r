source(here::here("scripts", "utils.r"))

library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)

# -----------------------------------------------------------------------------

protein <- "1dd3_A"

files <- list.files(
  path = paths$dist,
  pattern = paste0(protein, "_.*.csv"),
  full.names = TRUE
)

data <- list(
  com = list(
    noh = read.csv(grep("noh", files, value = TRUE)),
    bb = read.csv(grep("backbone", files, value = TRUE)),
    ca = read.csv(grep("calpha", files, value = TRUE))
  ),
  pai = read.csv(here::here("outputs", "PAI", "1dd3_A_angles.csv"))
)

# scatter of com-razdalje vs prvi principal axis, every frame coloured by
# how crowded its bin is
# bins are coarse on purpose: 1000 frames over 30 bins per axis would leave
# most bins with a single frame and the colours would all look the same
num_bins <- 15

pai <- data$pai$R1_p1
com <- data$com$ca$R1

# [p]lot [d]ata
pd <- list(
  R1 = list(pai = data$pai$R1_p1, com = data$com$ca$R1),
  R2 = list(pai = data$pai$R2_p1, com = data$com$ca$R2),
  R3 = list(pai = data$pai$R3_p1, com = data$com$ca$R3)
)

# -----------------------------------------------------------------------------
# [p]lot [s]catter
plot_scat <- \(pai, com) {
  # WARN: zaokroži na 2 decimalki
  pai <- round(pai, 2)

  # NOTE: mogoče bi razdalje lahko celo na 1..?
  com <- round(com, 2)

  # določi breaks za bins
  bounds <- list(
    pai = seq(min(pai), max(pai), length.out = num_bins + 1),
    com = seq(min(com), max(com), length.out = num_bins + 1)
  )

  # razdeli vrednosti po bins
  bins <- list(
    pai = cut(pai, breaks = bounds$pai, include.lowest = TRUE),
    com = cut(com, breaks = bounds$com, include.lowest = TRUE)
  )

  # izračunaj gostoto (frekvence za vsak bin)
  m <- bins |> table()

  # izriši scatter plot obarvan glede na gostoto (n)
  d <- data.frame(
    pai = pai,
    com = com,
    n = m[cbind(as.integer(bins$pai), as.integer(bins$com))]
  )

  # najgostejši bins zadnji, da so vidni na vrhu
  d <- arrange(d, n)

  d |>
    ggplot(aes(x = com, y = pai, fill = n)) +
    geom_point(shape = 21, size = 6, stroke = 0.5, alpha = 0.8) +
    scale_fill_viridis(option = "viridis", name = "n") +
    labs(x = "com (Å)", y = "pai (rad)") +
    theme(
      legend.position = "right",
      legend.key.height = grid::unit(2, "lines")
    ) +
    theme_bw()
}

# -----------------------------------------------------------------------------
# one panel per replicate
wrap_plots(
  lapply(names(pd), \(r) plot_scat(pd[[r]]$pai, pd[[r]]$com)),
  nrow = 3,
)

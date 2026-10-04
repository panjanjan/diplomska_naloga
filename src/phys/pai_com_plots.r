source(here::here("src", "utils.r"))

library(dplyr)
library(parallel)
library(ggplot2)
library(patchwork)
library(viridis)

# -----------------------------------------------------------------------------

protein <- "1dd3_A"
target <- here::here("outputs", "COM_PAI")
if (!dir.exists(target)) dir.create(target)
n_cores <- min(detectCores() - 1, 8)

# WARN: samo CA
data <- list(
  dist = list.files(paths$dist, "_calpha_dist.csv", full.names = TRUE),
  angles = list.files(paths$angles, ".csv", full.names = TRUE),
  domains = read.csv(paths$domains)
)

# scatter of com-pai vs prvi principal axis, vsak frame pobarvan glede na to
# kako gost je.
# NOTE: nek scaling razdalj in kotov glede na nekaj (npr. velikost proteina)
# bi pomagal...
num_bins <- 15

# -----------------------------------------------------------------------------
# [p]lot [s]catter
plot_scat <- \(pai, com, xlab, ylab, tit, subtit) {
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

  p <- d |>
    ggplot(aes(x = com, y = pai, fill = n)) +
    geom_point(shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
    scale_fill_viridis(option = "viridis", name = "n") +
    # labs(x = "com (Å)", y = "pai (rad)") +
    labs(x = xlab, y = ylab, title = tit, subtitle = subtit) +
    theme(
      legend.position = "right",
      legend.key.height = grid::unit(2, "lines")
    ) +
    theme_bw()

  return(p)
}

# -----------------------------------------------------------------------------
invisible(mclapply(unique(data$domains$protein), \(protein) {
  # WARN: samo prve osi
  pai <- read.csv(grep(protein, data$angles, value = TRUE))
  com <- read.csv(grep(protein, data$dist, value = TRUE))

  # [p]lot [d]ata
  pd <- list(
    R1 = list(pai = pai$R1_p1, com = com$R1),
    R2 = list(pai = pai$R2_p1, com = com$R2),
    R3 = list(pai = pai$R3_p1, com = com$R3)
  )

  # one panel per replicate
  wrap_plots(
    plot_scat(pd$R1$pai, pd$R1$com, "", "", paste(protein, "C-alpha"), "R1"),
    plot_scat(pd$R2$pai, pd$R2$com, "", "pai (rad)", "", "R2"),
    plot_scat(pd$R3$pai, pd$R3$com, "com (Å)", "", "", "R3"),
    nrow = 3
  )

  out <- file.path(target, paste0(protein, "_COM_PAI_ca.png"))
  ggsave(out)
}, mc.cores = n_cores))

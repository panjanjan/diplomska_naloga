#!/bin/Rscript

source(here::here("scripts", "utils.r"))

# -----------------------------------------------------------------------------
# avtokorelacija razdalj med masnima centroma skozi trajektorijo
# vse replike in vse datoteke v enem grafu

files <- list.files(paths$dist, "1k5n_A_.*.csv", full.names = TRUE)
# files <- list.files(paths$dist, "7zh9_A_.*.csv", full.names = TRUE)
reps <- c("R1", "R2", "R3")

# barva po repliki. datoteke (noh/backbone/calpha) rišemo čez isto barvo,
# ker dajo praktično identično avtokorelacijo
cols <- c(R1 = "#D55E00", R2 = "#0072B2", R3 = "#009E73")

# ploščina od y = 0 do vrednosti acf, kot geom_area v ggplot2, prosojna, da
# se prekrivajoči polinomi razlikujejo
fill_alpha <- 0.05

# polygon(x, y) bi se zaključil s poševnico med prvo in zadnjo točko, zato
# y = 0 dodamo sami
area <- function(x, y, col) {
  n <- length(y)
  polygon(c(x, rev(x)), c(y, rep(0, n)), col = col, border = NA)
}

# acf z offsetom ni drugačen, zato razdaljo razmeraj ne potrebujemo
acfs <- list()
for (file in files) {
  d <- read.csv(file)
  for (rep_col in reps) {
    acfs[[rep_col]][[basename(file)]] <-
      acf(d[[rep_col]], lag = length(d[[rep_col]]), plot = FALSE)$acf
  }
}

lags <- seq_along(acfs[[reps[1]]][[1]]) - 1
ylim <- range(unlist(acfs), 0) * 1.05

par(
  mar = c(4.1, 4.1, 3.0, 1.0),
  mgp = c(2.6, 0.7, 0),
  tcl = -0.3,
  cex.axis = 0.8,
  cex.lab = 0.8,
  cex.main = 0.9
)

plot(
  lags,
  acfs[[reps[1]]][[1]],
  type = "n",
  ylim = ylim,
  xlab = "lag (frame)",
  ylab = "acf",
  main = "avtokorelacija razdalje med masnima centroma"
)
grid(col = "grey90")
abline(h = 0, col = "grey30", lwd = 1)

# najprej vse ploščine, nato vse črte, da črte ne pridejo pod ploščine
for (rep_col in reps) {
  for (name in names(acfs[[rep_col]])) {
    area(lags, acfs[[rep_col]][[name]],
         adjustcolor(cols[[rep_col]], alpha.f = fill_alpha))
  }
}
for (rep_col in reps) {
  for (name in names(acfs[[rep_col]])) {
    lines(lags, acfs[[rep_col]][[name]], col = cols[[rep_col]], lwd = 1)
  }
}

legend(
  "topright",
  legend = names(cols),
  col = cols,
  lwd = 2,
  pch = 15,
  pt.cex = 1.6,
  bty = "n",
  cex = 0.8
)

# -----------------------------------------------------------------------------
# files <- list.files(paths$dist, "7zh9_A_.*.csv", full.names = TRUE)
#
# par(mfrow = c(3, 2))
# for (i in 1:3) {
#   d <- read.csv(files[1])[, paste0("R", i)]
#   d2 <- d - d[1]
#   plot(d2,
#        type = "l",
#        lwd = 2,
#        panel.first = grid(),
#        main = basename(files[i]),
#        xlab = "frame",
#        ylab = "razdalja (Å)") |>
#     abline(h = mean(d2), col = "red", lwd = 2)
#   acf(d2, lag = length(d2), lwd = 2, main = "avtokorelacija")
# }

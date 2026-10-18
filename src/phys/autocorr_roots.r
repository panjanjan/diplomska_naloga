# zaznavanje oscilacij v ACF vrednostih
library(magrittr)
source(here::here("src", "utils.r"))

# --- ATTEMPT #1: razdalje med ničlami ---
# če se ničle pojavljajo enakomerno narazen => oscilacije?

# lep primer je:
protein <- "4fch_B"

# `comparing_selections.r` pove, da sta bb in ca podobna, noh ima večje razlike:
#             mean        std
# bb_ca  0.1288206 0.01381483
# bb_noh 1.0284324 0.06822933
# ca_noh 1.1553076 0.06665335

sel <- "backbone"

d <- list.files(
  path = paths$dist,
  pattern = paste0(protein, ".*", sel),
  full.names = TRUE
) |>
  read.csv()

str(d)

abs_dist <- \(x, y) abs(x - y)
err <- \(s) c(mean = mean(s), std = sd(s))

# a je sploh kaj razlike med replikati
combs <- combn(names(d)[-1], 2)
res <- combs |>
  apply(2, \(col) {
    s1 <- d[, col[1]]
    s2 <- d[, col[2]]
    abs_dist(s1, s2) |> err()
  })
res <- rbind(res, res[1, ] - res[2, ])
res <- rbind(res, res[1, ] + res[2, ])
rownames(res) <- c("mean", "std", "min", "max")
colnames(res) <- apply(combs, 2, paste, collapse = "-")

# pa ja..... ampak zelo zelo male razlike
#          R1-R2     R1-R3     R2-R3
# mean 0.7118964 0.7713013 0.4627994
# std  0.5203206 0.5201582 0.3628522
# min  0.1915758 0.2511432 0.0999472
# max  1.2322170 1.2914595 0.8256516

# okej, acf in nato ničle
acfs <- d[, -1] %>%
  apply(2, acf, plot = FALSE, lag.max = nrow(d) * 0.5) %>%
  lapply(., \(x) x$acf) %>%
  as.data.frame()

# najdem indekse, kjer se spremeni predznak ACF
# vključim prvi indeks
zeros <- function(x) c(1, which(diff(sign(x)) != 0))
z <- apply(acfs, 2, zeros)

plot_zeros <- function(r, z, gap = 0) {
  s <- acfs[, r]
  n <- length(s)

  # tam kjer so nični indeksi naj bo rdeče, sicer naj se ne obarva
  c <- rep("transparent", n)
  c[z[[r]]] <- "red"

  plot(s, pch = 19, main = r, ylab = "ACF", panel.first = grid(), frame = FALSE)
  lines(s)
  points(1:n, rep(0, n), col = c, pch = 19)
  abline(h = 0)

  if (gap != 0) {
    abline(v = z[[r]] - gap, lty = 3)
    abline(v = z[[r]] + gap, lty = 3)
  }
}

par(mfrow = c(1, 3))
plot_zeros("R1", z)
plot_zeros("R2", z)
plot_zeros("R3", z)

# najdem predstavnika gruč indeksov; novo ničlo
# indeksi oddaljeni za manj kot gap se gručijo
# gruče so prisotne zaradi šuma v vhodnih podatkih (razdalje)
min_gap <- 10

group_zeros <- function(x, gap = 0) {
  if (gap == 0) return(NA)
  grp <- cumsum(c(TRUE, diff(x) > gap))
  tapply(x, grp, median) |>
    round()
}

z_grp <- lapply(z, group_zeros, gap = min_gap)

par(mfrow = c(1, 3))
plot_zeros("R1", z_grp, gap = min_gap)
plot_zeros("R2", z_grp, gap = min_gap)
plot_zeros("R3", z_grp, gap = min_gap)

# rabim vsaj 3 ničle na replikat, tako se mi zdi prav
# 4 indekse vključeno s prvim
z_ok <- z_grp[which(unlist(lapply(z_grp, length)) >= 4)]
repl_ok <- names(z_ok)

# ali so razdalje med indeksi podobne? dopuščam neko napako

# --- opisne statistike ---
# CV (coefficient of variation): globalna, opisuje celoten vzorec
# en osamelec ne popači rezultata
coeff_var <- function(x) sd(x) / mean(x)

# Gap Tolerance: lokalna, vsaka obravnavana razdalja mora ležati znotraj
# tolerance od povprečja vzorca.
gap_tol <- function(x, tol = 0.1) abs(x - median(x)) <= tol * median(x)

z_stats <- lapply(z_ok, \(x) {
  dx <- diff(x)
  data.frame(
    CV = coeff_var(dx),
    GT = gap_tol(dx)
  )
})

z_stats

# ali so približno enako oddaljeni indeksi...?
lapply(z_stats, \(x) all(x[, "GT"]))

# TODO: vsi proteini
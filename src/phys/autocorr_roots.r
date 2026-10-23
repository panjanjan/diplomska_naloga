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
repls <- grep("R", names(d), value = TRUE)

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
get_acf <- function(df) {
  df %>%
    apply(2, acf, plot = FALSE, lag.max = nrow(df) * 0.5) %>%
    lapply(., \(x) x$acf) %>%
    as.data.frame()
}
acfs <- get_acf(d[, -1])

# najdem indekse, kjer se spremeni predznak ACF
# vključim prvi indeks
zeros <- function(x) c(1, which(diff(sign(x)) != 0))
z <- apply(acfs, 2, zeros)

plot_zeros <- function(acfs, r, z, gap = 0, ...) {
  s <- acfs[, r]
  n <- length(s)

  # tam kjer so nični indeksi naj bo rdeče, sicer naj se ne obarva
  c <- rep("transparent", n)
  c[z[[r]]] <- "red"

  plot(
    s,
    pch = 19,
    main = r,
    panel.first = grid(),
    frame = FALSE,
    ...
  )
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
at_least_four <- function(l) {
  l %>%
    lapply(\(x) length(x) >= 4) %>%
    unlist() %>%
    which() %>%
    l[.]
}
z_ok <- at_least_four(z_grp)

# ali so razdalje med indeksi podobne? dopuščam neko napako

# --- opisne statistike ---
# CV (coefficient of variation): globalna, opisuje vse razdalje. En osamelec
# ne popači rezultata.
# * Visoke vrednosti => sd >> mean => visoka variacija razdalj
# * Nizke vrednosti ==> sd << mean => nizka variacija - to iščem
coeff_var <- function(x) sd(x) / mean(x)

# Gap Tolerance: lokalna, vsaka obravnavana razdalja mora ležati znotraj
# tolerance od povprečja razdalj
gap_tol <- function(x, tol = 0.3) abs(x - median(x)) <= tol * median(x)

gap_stats <- function(l, tol) {
  lapply(l, \(x) {
    dx <- diff(x)
    data.frame(
      CV = coeff_var(dx),
      GT = gap_tol(dx, tol)
    )
  })
}

z_stats <- gap_stats(z_ok)

# ali so približno enako oddaljeni indeksi...?

# ---- vsi proteini -----
library(parallel)
library(viridis)

sel <- "backbone"
min_gap <- 10
tolerance <- 0.3
colors <- viridis(n = 3)
d <- list.files(
  path = paths$dist,
  pattern = sel,
  full.names = TRUE
)
proteins <- read.csv(paths$domains)[, "protein"] |>
  unique() |>
  sort()
out <- here::here("outputs")
pic_dir <- file.path(out, "ACF_roots")
if (!dir.exists(pic_dir)) dir.create(pic_dir)
csv_out <- file.path(out, "autocorr_roots_test.csv")

# protein R1 R2 R3
# xxxx    T  F  T
run <- function(protein) {
  dprot <- grep(protein, d, value = TRUE) |> read.csv()
  repls <- grep("R", names(dprot), value = TRUE)
  acfs <- get_acf(dprot[, repls])
  z <- apply(acfs, 2, zeros)
  z_grp <- lapply(z, group_zeros, gap = min_gap)

  png(
    file.path(pic_dir, paste0(protein, "_acf_roots.png")),
    width = 1500,
    height = 600,
    res = 150
  )
  par(mfrow = c(1, 3))
  plot_zeros(acfs, "R1", z_grp, min_gap, xlab = "", ylab = "ACF", col = colors[1])
  plot_zeros(acfs, "R2", z_grp, min_gap, ylab = "", col = colors[2])
  plot_zeros(acfs, "R3", z_grp, min_gap, xlab = "", ylab = "", col = colors[3])
  dev.off()

  z_ok <- at_least_four(z_grp)
  z_stats <- gap_stats(z_ok, tolerance)

  # recimo, da je 0.2 dobra meja
  cv_cutoff <- 0.2

  # vrne bool vektor, TRUE če passa z vsemi GT in dobrim CV
  res <- lapply(z_stats, \(x) all(x$GT) || all(x$CV <= cv_cutoff)) |> unlist()
  final <- rep(FALSE, 3)
  names(final) <- repls
  for (i in seq_along(res)) {
    final[names(res[i])] <- res[i]
  }

  data.frame(
    protein = protein,
    R1 = final["R1"],
    R2 = final["R2"],
    R3 = final["R3"]
  )
}

lapply(proteins, run) %>%
  do.call(rbind, .) %>%
  write.csv(., csv_out, quote = FALSE, row.names = FALSE)

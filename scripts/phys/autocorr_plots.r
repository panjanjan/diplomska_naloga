source(here::here("scripts", "utils.r"))

# prototipiram na proteinu za katerega vem da se nekaj dogaja
protein <- "1k5n_A"
files <- list.files(
  path = paths$dist,
  pattern = paste0(protein, "_.*.csv"),
  full.names = TRUE
)

data <- list(
  noh = read.csv(grep("noh", files, value = TRUE)),
  bb = read.csv(grep("backbone", files, value = TRUE)),
  ca = read.csv(grep("calpha", files, value = TRUE))
)

n <- nrow(data$noh)

aux_acf <- \(x) acf(x - mean(x), lag = n, pl = FALSE)[["acf"]]

acf_l <- list(
  lags = seq(0, n - 1),
  noh = data.frame(
    R1 = aux_acf(data$noh$R1),
    R2 = aux_acf(data$noh$R2),
    R3 = aux_acf(data$noh$R3)
  ),
  bb = data.frame(
    R1 = aux_acf(data$bb$R1),
    R2 = aux_acf(data$bb$R2),
    R3 = aux_acf(data$bb$R3)
  ),
  ca = data.frame(
    R1 = aux_acf(data$ca$R1),
    R2 = aux_acf(data$ca$R2),
    R3 = aux_acf(data$ca$R3)
  )
)

str(acf_l)

rep_colors <- c(R1 = "#ff8e32", R2 = "#cd68bb", R3 = "#51c3cc")
rep_nums <- names(rep_colors)
sel <- "ca"
line_w <- 1
line_alpha <- 0.2

dark_plot()
par(mfrow = c(2, 1))
sapply(seq_along(rep_colors), \(i) {
  dist_vals <- data[[sel]][[rep_nums[i]]]
  dist_vals <- dist_vals - mean(dist_vals)
  color <- rep_colors[i]
  mod_color <- adjustcolor(color, alpha.f = line_alpha)
  main <- paste0(
    "razdalje masnih centrov, protein: ",
    protein,
    ", selekcija: ",
    toupper(sel)
  )
  if (i == 1) {
    plot(
      1:n,
      dist_vals,
      type = "n",
      xlab = "frame",
      ylab = "razdalja (Å)",
      main = main,
      panel.first = grid()
    )
  }
  polygon(
    x = c(1:n, n:1),
    y = c(dist_vals, rep(0, length(dist_vals))),
    col = mod_color,
    border = color,
    lwd = line_w
  )
  # lines(dist_vals, lwd = line_w, col = color)
  # points(dist_vals, lwd = 0.05, col = mod_color, pch = 19)
}) +
sapply(seq_along(rep_colors), \(i) {
  acf_vals <- acf_l[[sel]][[rep_nums[i]]]
  lags <- acf_l$lags
  color <- rep_colors[i]
  mod_color <- adjustcolor(color, alpha.f = line_alpha)
  main <- "avtokorelacije razdalj"
  if (i == 1) {
    plot(
      lags,
      acf_vals,
      type = "n",
      xlab = "zamik",
      ylab = "ACF",
      main = main,
      panel.first = grid()
    )
  }
  polygon(
    x = c(lags, rev(lags)),
    y = c(acf_vals, rep(0, length(acf_vals))),
    col = mod_color,
    border = color,
    lwd = line_w
  )
})
abline(h = 0, lwd = line_w) +
legend(
  "topright",
  legend = rep_nums,
  col = rep_colors,
  lwd = line_w,
  pch = 15,
  pt.cex = 1.6,
  bty = "n",
  cex = 0.8
)

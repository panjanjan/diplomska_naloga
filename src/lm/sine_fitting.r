source(here::here("src", "utils.r"))

d <- here::here("outputs/COM/1dd3_A_backbone_dist.csv") |>
  read.csv()
d <- d[, -1]
n <- nrow(d)
x <- 1:n

run <- function(r) {
  y <- d[, r]

  # ------ regresija ----------------------------------------------------------
  # scan every frequency: for a FIXED f the model
  #
  #     y = a*sin(2*pi*f*x) + b*cos(2*pi*f*x) + c
  #
  # is linear in (a, b, c), so least squares gives the exact best sine at
  # that f. sqrt(a^2 + b^2) is the amplitude of that frequency

  # sampling: 1 frame = 0.1 ns  ->  Nyquist = 0.5 cycles/frame = 5 GHz
  fs <- seq(1 / n, 0.5, length.out = n / 2) # cycles/frame

  amp_scan <- t(vapply(fs, function(f) {
    si <- sin(2 * pi * f * x)
    co <- cos(2 * pi * f * x)
    fit <- lm(y ~ si + co)
    cf <- coef(fit)
    c(
      amp = sqrt(cf[["si"]]^2 + cf[["co"]]^2),
      r2 = summary(fit)$r.squared,
      phi = atan2(cf[["co"]], cf[["si"]])
    )
  }, numeric(3)))

  freq_ghz <- fs * 10 # cycles/frame -> cycles/ns = GHz
  period_ns <- 1 / freq_ghz

  # local maxima of the amplitude scan = frequencies that are "present"
  amp <- amp_scan[, "amp"]
  inner <- 2:(length(amp) - 1)
  peaks <- inner[amp[inner] > amp[inner - 1] & amp[inner] > amp[inner + 1]]
  peaks <- peaks[order(amp[peaks], decreasing = TRUE)]
  top <- head(peaks, 5)

  cat("\ntop 5 frequencies:\n\n")
  peaks_df <- data.frame(
    freq_ghz = round(freq_ghz[top], 3),
    period_ns = round(period_ns[top], 2),
    amp = round(amp[top], 3),
    r2 = round(amp_scan[top, "r2"], 3)
  )
  print(peaks_df)

  # refine the dominant frequency with nls, starting from the scan result
  best <- peaks[1]
  fit <- nls(
    y ~ a * sin(2 * pi * f * x + phi) + c,
    data = data.frame(x = x, y = y),
    start = list(
      a = unname(amp_scan[best, "amp"]),
      f = fs[best],
      phi = unname(amp_scan[best, "phi"]),
      c = mean(y)
    )
  )

  fit
}

rs <- list(
  R1 = run("R1"),
  R2 = run("R2"),
  R3 = run("R3")
)

png(
  filename = here::here("outputs", "1dd3_A_sine_fitting.png"),
  width = 1700,
  height = 1900,
  res = 180
)
par(mfrow = c(3, 1))
for (i in 1:3) {
  plot(
    d[, i] - mean(d[, i]),
    type = "l",
    panel.first = grid(),
    lwd = 2,
    col = "lightgray",
    main = paste0("1dd3_A R", i),
    xlab = "frame",
    ylab = "COM distance (A)"
  )
  lines(
    predict(rs[[paste0("R", i)]]) - mean(predict(rs[[paste0("R", i)]])),
    lwd = 2
  )
}
dev.off()
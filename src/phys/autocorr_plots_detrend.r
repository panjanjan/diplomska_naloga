source(here::here("src", "utils.r"))

files <- list.files(path = paths$dist, pattern = "_.*.csv", full.names = TRUE)
target <- here::here("outputs", "ACF_detrend")
if (!dir.exists(target)) dir.create(target, recursive = TRUE)

# -----------------------------------------------------------------------------

plot_sel <- function(dist_sel, acf_sel, acf_lags, sel, protein, frames, n) {
  cfg <- list(
    line_w = 1.5,
    area_alpha = 0.15,
    font_scale = 1.5,
    margin = 0,
    rep_colors = c(R1 = "#ff8e32", R2 = "#cd68bb", R3 = "#51c3cc")
  )
  cfg$mod_colors <- sapply(
    X = cfg$rep_colors,
    FUN = adjustcolor,
    alpha.f = cfg$area_alpha
  )
  cfg$rep_nums <- names(cfg$rep_colors)

  # razdalje
  main <- paste0(
    "Razdalje masnih centrov, protein: ",
    protein,
    ", selekcija: ",
    toupper(sel)
  )

  # da nariše vse replikate skupaj, najde min in max med vsemi replikati za
  # y-limite
  d_ylims <- c(
    min(dist_sel) - cfg$margin,
    max(dist_sel) + cfg$margin
  )

  par(mfrow = c(2, 1))
  plot(
    frames,
    rep(0, length(frames)),
    ylim = d_ylims,
    type = "n",
    xlab = "frame",
    ylab = "razdalja (Å)",
    # main = main,
    main = "",
    panel.first = grid(),
    cex = cfg$font_scale,
    frame.plot = FALSE
  )
  for (i in 1:3) {
    polygon(
      x = c(1:n, n:1),
      y = c(dist_sel[, i], rep(0, length(dist_sel[, i]))),
      col = cfg$mod_colors[i],
      border = cfg$rep_colors[i],
      lwd = cfg$line_w,
      cex = cfg$font_scale
    )
  }

  # avtokorelacije
  a_ylims <- c(
    min(acf_sel) - cfg$margin,
    max(acf_sel) + cfg$margin
  )
  plot(
    acf_lags,
    rep(0, length(acf_lags)),
    ylim = a_ylims,
    type = "n",
    xlab = "zamik",
    ylab = "ACF",
    # main = "Avtokorelacije razdalj",
    main = "",
    panel.first = grid(),
    cex = cfg$font_scale,
    frame.plot = FALSE
  )
  for (i in 1:3) {
    polygon(
      x = c(acf_lags, rev(acf_lags)),
      y = c(acf_sel[, i], rep(0, length(acf_sel[, i]))),
      col = cfg$mod_colors[i],
      border = cfg$rep_colors[i],
      lwd = cfg$line_w,
      cex = cfg$font_scale
    )
  }

  legend(
    "topright",
    legend = cfg$rep_nums,
    col = cfg$rep_colors,
    lwd = cfg$line_w,
    pch = 15,
    pt.cex = 1.6,
    bty = "n",
    cex = cfg$font_scale
  )
}

run <- function(protein) {
  protein_dist <- grep(files, pattern = protein, value = TRUE)
  stopifnot(length(protein_dist) == 3)

  dist_l <- list(
    noh = read.csv(grep("noh", protein_dist, value = TRUE)),
    bb = read.csv(grep("backbone", protein_dist, value = TRUE)),
    ca = read.csv(grep("calpha", protein_dist, value = TRUE))
  )

  # make sure, da imajo vsi enako število frame-ov
  stopifnot(length(unique(sapply(dist_l, nrow))) == 1)
  n <- nrow(dist_l[[1]])
  frames <- dist_l[[1]][["frame"]]

  # NOTE: detrend
  dist_l <- list(
    noh = data.frame(
      R1 = dist_l$noh$R1 - mean(dist_l$noh$R1),
      R2 = dist_l$noh$R2 - mean(dist_l$noh$R2),
      R3 = dist_l$noh$R3 - mean(dist_l$noh$R3)
    ),
    bb = data.frame(
      R1 = dist_l$bb$R1 - mean(dist_l$bb$R1),
      R2 = dist_l$bb$R2 - mean(dist_l$bb$R2),
      R3 = dist_l$bb$R3 - mean(dist_l$bb$R3)
    ),
    ca = data.frame(
      R1 = dist_l$ca$R1 - mean(dist_l$ca$R1),
      R2 = dist_l$ca$R2 - mean(dist_l$ca$R2),
      R3 = dist_l$ca$R3 - mean(dist_l$ca$R3)
    )
  )

  aux_acf <- \(x) acf(x, lag = n, pl = FALSE)[["acf"]]

  acf_l <- list(
    lags = seq(0, n - 1),
    noh = data.frame(
      R1 = aux_acf(dist_l$noh$R1),
      R2 = aux_acf(dist_l$noh$R2),
      R3 = aux_acf(dist_l$noh$R3)
    ),
    bb = data.frame(
      R1 = aux_acf(dist_l$bb$R1),
      R2 = aux_acf(dist_l$bb$R2),
      R3 = aux_acf(dist_l$bb$R3)
    ),
    ca = data.frame(
      R1 = aux_acf(dist_l$ca$R1),
      R2 = aux_acf(dist_l$ca$R2),
      R3 = aux_acf(dist_l$ca$R3)
    )
  )

  # -- [plotting] --------------------------------------------------------------
  for (sel in names(dist_l)) {
    fname <- paste0(protein, "_", sel, "_acf.png")
    out <- file.path(target, fname)
    png(out, width = 700, height = 600)
    plot_sel(dist_l[[sel]], acf_l[[sel]], acf_l$lags, sel, protein, frames, n)
    dev.off()
  }
}

# -----------------------------------------------------------------------------

for (protein in read.csv(paths$domain)[["protein"]]) {
  cat(protein, "\r")
  run(protein)
}
cat("\n")

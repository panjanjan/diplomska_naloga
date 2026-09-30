source(here::here("scripts", "utils.r"))

# prototipiram na proteinu za katerega vem da se nekaj dogaja
protein <- "1k5n_A"
sel <- "ca"

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

aux_acf <- \(x) acf(x, lag = n, pl = FALSE)[["acf"]]

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

# -- [plotting] ---------------------------------------------------------------
rep_colors <- c(R1 = "#ff8e32", R2 = "#cd68bb", R3 = "#51c3cc")
rep_nums <- names(rep_colors)
line_w <- 2
line_alpha <- 0.6
area_alpha <- 0.15
font_scale <- 1.5
pic_dims <- c(1496, 1115)

plot_sel <- function(sel) {
  # -- [1. razdalje] ------------------------------------------------------------
  frames <- data[[sel]][["frame"]]
  d <- data.frame(
    R1 = data[[sel]]$R1 - data[[sel]]$R1[1],
    R2 = data[[sel]]$R2 - data[[sel]]$R2[1],
    R3 = data[[sel]]$R3 - data[[sel]]$R3[1]
  )

  main <- paste0(
    "Razdalje masnih centrov, protein: ",
    protein,
    ", selekcija: ",
    toupper(sel)
  )

  # da nariše vse replikate skupaj, najde min in max med vsemi replikati za
  # y-limite
  margin <- 0
  ylims = c(min(d) - margin, max(d) + margin)

  par(mfrow = c(2, 1))
  plot(
    frames,
    rep(0, length(frames)),
    ylim = ylims,
    type = "n",
    xlab = "frame",
    ylab = "razdalja (Å)",
    main = main,
    panel.first = grid(),
    cex = font_scale,
    frame.plot = FALSE
  )
  for (i in 1:3) {
    mod_color <- adjustcolor(rep_colors[i], alpha.f = area_alpha)
    polygon(
      x = c(1:n, n:1),
      y = c(d[, i], rep(0, length(d[, i]))),
      col = mod_color,
      border = rep_colors[i],
      lwd = line_w,
      cex = font_scale
    )
    # lines(d[, i], col = mod_color, lwd = line_w)
    # points(d[, i], col = mod_color, pch = 19, lwd = 1)
  }
  abline(h = 0)

  # -- [2. avtokorelacije] ------------------------------------------------------
  sapply(seq_along(rep_colors), \(i) {
    acf_vals <- acf_l[[sel]][[rep_nums[i]]]
    lags <- acf_l$lags
    color <- rep_colors[i]
    mod_color <- adjustcolor(color, alpha.f = area_alpha)
    main <- "Avtokorelacije razdalj"
    if (i == 1) {
      plot(
        lags,
        acf_vals,
        type = "n",
        xlab = "zamik",
        ylab = "ACF",
        main = main,
        panel.first = grid(),
        frame.plot = FALSE,
        cex = font_scale
      )
    }
    polygon(
      x = c(lags, rev(lags)),
      y = c(acf_vals, rep(0, length(acf_vals))),
      col = mod_color,
      border = color,
      lwd = line_w,
      cex = font_scale
    )
  })
  abline(h = 0)
  legend(
    "topright",
    legend = rep_nums,
    col = rep_colors,
    lwd = line_w,
    pch = 15,
    pt.cex = 1.6,
    bty = "n",
    cex = font_scale
  )
}

png("yo2.png", width = pic_dims[1], height = pic_dims[2])
plot_sel("ca")
dev.off()

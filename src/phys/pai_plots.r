library(ggplot2)
source(here::here("src", "utils.r"))

files <- list.files(paths$angles, ".csv", full.names = TRUE)
target <- here::here("outputs", "PAI")
if (!dir.exists(target)) dir.create(target, recursive = TRUE)

for (filename in files) {
  cat(filename, "\n")
  d <- read.csv(filename)
  p <- sub("_angles.csv", "", basename(filename))
  n <- nrow(d)
  data <- list(
    R1 = d[, "R1_p1"],
    R2 = d[, "R2_p1"],
    R3 = d[, "R3_p1"]
  )
  # samo prva os for now
  acfs <- list(
    R1 = acf(data$R1, lag.max = n * 0.5, plot = FALSE),
    R2 = acf(data$R2, lag.max = n * 0.5, plot = FALSE),
    R3 = acf(data$R3, lag.max = n * 0.5, plot = FALSE)
  )
  # png(file.path(target, paste0(p, "_angle_acf.png")), width = 800)
  par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))
  plot(data$R1, type = "l", frame = FALSE, panel.first = grid(), main = "PAI 1, R1", ylab = "angle", xlab = "frame")
  plot(data$R2, type = "l", frame = FALSE, panel.first = grid(), main = "PAI 1, R2", ylab = "angle", xlab = "frame")
  plot(data$R3, type = "l", frame = FALSE, panel.first = grid(), main = "PAI 1, R3", ylab = "angle", xlab = "frame")
  plot(acfs$R1, type = "l", frame = FALSE, panel.first = grid(), main = "ACF, R1", ylab = "ACF", xlab = "lag")
  abline(h = 0)
  plot(acfs$R2, type = "l", frame = FALSE, panel.first = grid(), main = "ACF, R2", ylab = "ACF", xlab = "lag")
  abline(h = 0)
  plot(acfs$R3, type = "l", frame = FALSE, panel.first = grid(), main = "ACF, R3", ylab = "ACF", xlab = "lag")
  abline(h = 0)
  mtext(p, outer = TRUE, side = 3, line = 0)
  # dev.off()
}

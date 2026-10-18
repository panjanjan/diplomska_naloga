library(dplyr)
library(viridis)
library(parallel)
source(here::here("src", "utils.r"))

# imena proteinov
proteins <- read.csv(paths$domains)$protein %>%
  unique() %>%
  {
    names(.) <- .
    .
  }

res <- mclapply(proteins, \(protein) {
  # dataframes razdalj
  files <- list.files(paths$dist, protein, full.names = TRUE)
  com_dist <- lapply(files, read.csv)
  names(com_dist) <- basename(files)

  # izračuna evklidske razdalje za vsak frame med x in y vzorcem; x je
  # en nabor atomov, y drug. vrne povprečje in standardno deviacijo razlik
  #
  # NOTE: a ni samo absolutna vrednost dovolj, ker so skalarji?
  euc_dist <- function(x, y) {
    (x - y)**2 |>
      rowSums() |>
      sqrt()
  }
  err <- function(s) c(mean = mean(s), std = sd(s))

  combs <- com_dist |>
    names() |>
    combn(2)

  combs %>%
    apply(2, \(col) {
      x <- com_dist[[col[1]]]
      y <- com_dist[[col[2]]]
      r <- euc_dist(x, y) %>% err()
    }) %>%
    {
      colnames(.) <- combs %>%
        apply(2, \(col) {
          gsub(".*_.*_(.*)_dist.csv", "\\1", col, perl = TRUE) %>%
            paste(collapse = "-")
        })
      .
    }
}, mc.cores = 10)

# ---- [CSVs] ------------------------------------------------------------------
# 3 files
# targets:
#  outputs/sel_diff_backbone-calpha.csv
#  outputs/sel_diff_backbone-noh.csv
#  outputs/sel_diff_calpha-noh.csv
# format:
#  protein,mean,std
#  protein1,xxx,yyy
#  protein2,xxx,yyy
#  protein3,xxx,yyy

target <- here::here("outputs")
sels <- colnames(res[[1]])

sapply(sels, \(sel) {
  out <- paste0("sel_diff_", sel, ".csv")
  csv <- file.path(target, out)

  # povleči vrednosti ene kombinacije
  res <- lapply(res, \(m) m[, sel]) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  # premakni proteine iz rownames v svoj stolpec, da ima CSV stolpec "protein"
  res %>%
    mutate(protein = rownames(res)) %>%
    {
      rownames(.) <- NULL
      .
    } %>%
    select(protein, mean, std) %>%
    write.csv(., csv, quote = FALSE, row.names = FALSE)
})

# --- [plotting] ---------------------------------------------------------------
colors <- viridis(n = length(res))

par(mfrow = c(3, 1), oma = c(0, 0, 2, 0))

sapply(colnames(res[[1]]), \(col) {
  lapply(res, \(m) m["mean", col]) %>%
    unlist() %>%
    plot(
      type = "h",
      lwd = 10,
      xaxt = "n",
      xlab = "",
      ylab = "value",
      main = col,
      col = colors,
      frame.plot = FALSE,
      panel.first = grid()
    ) +
    axis(
      1,
      at = seq_along(res),
      labels = names(res),
      las = 2,
      lwd = 0,
      cex.axis = 1
    )
})
mtext(
  "Povprečja Evklidskih razdalj med selekcijami",
  outer = TRUE,
  side = 3,
  line = 0
)

sapply(colnames(res[[1]]), \(col) {
  lapply(res, \(m) m["std", col]) %>%
    unlist() %>%
    plot(
      type = "h",
      lwd = 10,
      xaxt = "n",
      xlab = "",
      ylab = "value",
      main = col,
      col = colors,
      frame.plot = FALSE,
      panel.first = grid()
    ) +
    axis(
      1res,
      at = seq_along(res),
      labels = names(res),
      las = 2,
      lwd = 0,
      cex.axis = 1
    )
})
mtext(
  "Standardna deviacija Evklidskih razdalj med selekcijami",
  outer = TRUE,
  side = 3,
  line = 0
)

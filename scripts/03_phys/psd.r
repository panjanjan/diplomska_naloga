#!/bin/Rscript
# ...
# ------------------------------------------------------------------------------
library(parallel)

setwd(Sys.getenv("ROOT"))
source("./scripts/utils.r")

data <- load_data("dist")

target <- file.path("atlas_db", "PSD")
if (!dir.exists(target)) dir.create(target)

# ------------------------------------------------------------------------------
# * določi moči frekvenc v vhodnem signalu x, ki so razdalje med masnima centroma
#   skozi trajektorijo
# * vrne eno-vrstični data.frame s frekvencami in močmi
psd <- function(x) {
    n <- length(x)

    # nastavimo vse razdalje relativno glede na prvo
    # NOTE: mogoče bi to naredil že v shranjenih razdaljah..?
    x <- x - x[1]

    # frekvenca vzorčenja: kolikokrat na časovno enoto vzorčimo (1 ps)
    dt_sec <- 1 * 1e-12
    fs_hz <- 1 / dt_sec # v Hz

    # frekvenčni korak: razmik med točkami v frekvenčni domeni
    #
    # npr. za 100 meritev s frekvenco vzorčenja 1s je najnižja frekvenca, ki jo
    # lahko določimo 0.01 Hz = 1/100 s^-1 <=> 1 cikel na 100 s (samo 1 cikel skozi vse meritve)
    # najvišja frekvenca je enaka velikosti vzorca, n, in je 1 cikel na 1 s <=> 100/100 s^-1 = 1 Hz
    #
    # za primer 1001 meritev s frekvenco vzorčenja 1 ps (1e-12 s = 1/10000...0 s) je
    # najnižja določljiva frekvenca 1/1001 ps^-1 = 0.000999001 * 1e+12 = 0.999001 GHz,
    # najvišja pa 1001/1001 ps^-1 = 1 * 1e+12 Hz = 1 THz (načeloma...)
    df <- fs_hz / n

    # nyquistova frekvenca (N/2) določa do katere (analizne) frekvence se *aliasing*
    # ne dogaja. potrebujemo samo frekvence do te (included), da zajamemo vse.
    half_n <- 1 + floor(n / 2)
    freq <- seq(0, by = df, length.out = half_n)

    # ker merimo 1.001 ns, lahko pretvorimo frekvence iz Hz v ns (nHz...)
    # 999000999 / s = 999000999 / 1e+9 ns = 0.999001 ns^-1 = 0.999001 nHz
    freq_ns <- freq / 1e+9

    # predpostavka DFT je stacionarnost podatkov. odstranjevanje linearnega trenda
    fit <- lm(x ~ seq_along(x))
    x_detrend <- resid(fit)

    # načeloma naj bi imeli manjši proteini nihanja z višjimi frekvencami in večji z
    # manjšimi, zato jih scale-am na [0,1]
    scale_it <- \(v) (v - min(v)) / (max(v) - min(v))

    x_scaled <- scale_it(x_detrend)

    # TODO: spectral leakage?

    # power
    fft_out <- fft(x_scaled)
    power <- (Mod(fft_out)^2)[1:half_n] |> scale_it()

    data.frame(frequency = freq_ns, power = power)
}

# * csv z imenom {protein}_psd.csv
# * stolpci R1-R3
# * i-ta vrstica vsebuje moč i-te frekvence v signalu
#
# frequency power_R1 power_R2 power_R3
# 0         x        x        x
# 0.999001  x        x        x
# 1.998002  x        x        x
# ...
run <- function(i) {
    d <- read.csv(data[i])
    protein <- sub("_dist.csv", "", basename(data[i]))

    # moči za vsako frekvenco
    r1 <- psd(d$R1)
    r2 <- psd(d$R2)
    r3 <- psd(d$R3)

    powers <- data.frame(
        frequency = r1$frequency,
        power_R1  = r1$power,
        power_R2  = r2$power,
        power_R3  = r3$power
    )

    out <- paste0(protein, "_psd.csv")
    csv <- file.path(target, out)
    write.csv(powers, csv, quote = FALSE, row.names = FALSE)
}

# ------------------------------------------------------------------------------
n_cores <- min(detectCores() - 1, 10)
cat("using", n_cores, "cores\n")
invisible(mclapply(1:length(data), run, mc.cores = n_cores))

# NOTE: zastarelo
# NOTE: zastarelo
# NOTE: zastarelo
# glej rmsf_test.r

library(dplyr)
library(ggplot2)

# iskanje rmsf datotek izbranih proteinov
files <- list.files("../atlas_db/RMSF", pattern = "tsv", full.names = TRUE)
# files <- list.files("atlas_db/RMSF", pattern = "tsv", full.names = TRUE)

# proteini z dvema domenama
domains <- read.csv("two_domains.csv")

# preverjanje distribucij rmsf vrednosti
rmsf <- read.delim("atlas_db/RMSF/1dd3_A_RMSF.tsv")

par(mfrow = c(2, 2), mar = c(3, 3, 1, 1))
rmsf$RMSF_R1 |> hist(main = "")
rmsf$RMSF_R1 |> boxplot()
log(rmsf$RMSF_R1) |> hist(main = "")
log(rmsf$RMSF_R1) |> boxplot()

# očitno ne sledijo normalni porazdelitvi
# zelo so skewed, zato uporabim logaritemsko transformacijo

# velikosti proteinov (število vrstic v RMSF datotekah = št. aminokislin)
ns <- sapply(rmsf_filtered_files, \(x) nrow(as_tibble(read.delim(x))))

prot_group <- function(d, p) {
    d[d$protein == p, ]
}

### 1. KS-test
# ne zahteva log-transformacije
run <- function(protein) {
    # preberi rmsf
    rmsf <- files[grep(protein, files)] |> read.delim()
    rmsf <- rmsf[, -1]

    # dobi podatke od proteina
    d <- prot_group(domains, protein)

    pvals <- sapply(names(rmsf), \(col) {
        # najdi začetek in konec domen
        start <- d$start
        end <- d$end

        # matchaj rmsf
        d1 <- rmsf[start[1]:end[1], col]
        d2 <- rmsf[start[2]:end[2], col]

        # testiraj
        ks.test(d1, d2, alternative = "two.sided")$p.value
    })

    pvals
}

#+ warning=FALSE
test_results <- sapply(filtered_proteins, run) |> t()
test_results |> head()

# kateri passajo z vsemi replikati
keep3 <- which(apply(test_results, 1, \(row) all(row < 0.05)))

# kateri passajo z dvema ...
keep2 <- which(apply(test_results, 1, \(row) sum(row < 0.05) == 2))

# kateri passajo z enim ...
keep1 <- which(apply(test_results, 1, \(row) any(row < 0.05)))

c("3" = length(keep3), "2" = length(keep2), "1" = length(keep1))


### 2. log-transform
run <- function(protein) {
    # preberi rmsf
    rmsf <- files[grep(protein, files)] |> read.delim()
    rmsf <- rmsf[, -1]

    # log transform
    rmsf <- log(rmsf)

    # dobi podatke od proteina
    d <- prot_group(domains, protein)

    pvals <- sapply(names(rmsf), \(col) {
        # najdi začetek in konec domen
        start <- d$start
        end <- d$end

        # matchaj rmsf
        d1 <- rmsf[start[1]:end[1], col]
        d2 <- rmsf[start[2]:end[2], col]

        # testiraj
        t.test(d1, d2, alternative = "two.sided")$p.value
    })

    pvals
}

#+ warning=FALSE
test_results <- sapply(filtered_proteins, run) |> t()
test_results |> head()

# kateri passajo z vsemi replikati
keep3 <- which(apply(test_results, 1, \(row) all(row < 0.05)))

# kateri passajo z dvema ...
keep2 <- which(apply(test_results, 1, \(row) sum(row < 0.05) == 2))

# kateri passajo z enim ...
keep1 <- which(apply(test_results, 1, \(row) any(row < 0.05)))

c("3" = length(keep3), "2" = length(keep2), "1" = length(keep1))

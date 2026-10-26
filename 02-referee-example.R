## Discussion after Proposition 1: the graph Z -> X <- H -> Y <- W, with H unobserved.
## The district factor needs Z in the second factor: P(x | z) P(y | x, z, w).
## Dropping Z, i.e. using P(x | z) P(y | x, w), gives the wrong value.
source("common.R")

pH <- c("00" = 1/2, "01" = 1/6, "10" = 1/6, "11" = 1/6)   # distribution of H = (U, V)
pZW <- function(a) ifelse(a == 1, 1/3, 2/3)
tab <- CJ(z = 0:1, w = 0:1, u = 0:1, v = 0:1)
tab[, p := pZW(z) * pZW(w) * pH[paste0(u, v)]]
tab[, `:=`(x = (z + u) %% 2, y = (w + v) %% 2)]
P <- tab[, .(p = sum(p)), by = .(z, w, x, y)]
m <- function(...) { cond <- list(...); P[Reduce(`&`, Map(function(n, v) P[[n]] == v, names(cond), cond)), sum(p)] }

for (z in 0:1) {
  w <- 0
  wrong <- sum(sapply(0:1, \(x) m(x = x, z = z) / m(z = z) * m(x = x, y = 1, w = w) / m(x = x, w = w)))
  right <- sum(sapply(0:1, \(x) m(x = x, z = z) / m(z = z) * m(x = x, y = 1, z = z, w = w) / m(x = x, z = z, w = w)))
  cat(sprintf("z = %d, w = 0:  P(Y = 1) with Z dropped = %s,  correct = %s\n",
              z, MASS::fractions(wrong), MASS::fractions(right)))
}

section("Conditioning sets used by the software")
rg <- dhvm(caugi(Z %-->% X, X %-->% Y, W %-->% Y, H %-->% X + Y), latent.pattern = "^H")
print(dhvm_district_params(rg)$plist)

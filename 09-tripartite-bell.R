## Section 4.3: the tripartite Bell graph (Figure 7). This takes up to about half an hour.
source("common.R")

bell <- dhvm(caugi(X %-->% V1, Y %-->% V2, Z %-->% V3, U %-->% V1 + V2 + V3))
t <- system.time(cn <- meraconstraints(bell))
print(t)
H <- cn$polytope$Hrep[[1]]
cat("Total constraints:", nrow(H), "  inequalities:", sum(H[, 1] == 0),
    "  equalities:", sum(H[, 1] == 1), "\n")
fl <- cn$nontrivial[[1]]
cat("Nontrivial inequalities:", sum(fl & H[, 1] == 0), "  nontrivial equalities:", sum(fl & H[, 1] == 1), "\n")
k <- nonnegativity_entries(cn, 1)
cat("Nonnegativity constraints:", sum(!is.na(k)), "  distinct entries:", length(unique(k[!is.na(k)])), "\n")

section("Independencies")
print(cn$independences)

section("Equality constraints: the two-party margins do not depend on the third setting")
lev <- entry_levels(cn, 1)
E   <- H[H[, 1] == 1, -(1:2), drop = FALSE]
## For each pair of parties (a, b), with s the setting of the third: the margin of (a, b),
## summed over the other two settings' values held fixed, is the same at s = 0 and s = 1
cols <- list(c("V1", "V3", "Y"), c("V1", "V2", "Z"), c("V2", "V3", "X"))
rows <- do.call(rbind, lapply(cols, \(cc) {
  a <- cc[1]; b <- cc[2]; s <- cc[3]; others <- setdiff(c("X", "Y", "Z"), s)
  key <- do.call(paste, lev[c(a, b, others)])
  do.call(rbind, lapply(unique(key), \(u) {
    as.numeric(key == u & lev[[s]] == 1) - as.numeric(key == u & lev[[s]] == 0)
  }))
}))
N <- t(vapply(unique(do.call(paste, lev[c("X", "Y", "Z")])), \(u) as.numeric(do.call(paste, lev[c("X", "Y", "Z")]) == u), numeric(nrow(lev))))
cat("rank(H_e):", qr(E)$rank, "  rank(normalization, margins):", qr(rbind(N, rows))$rank,
    "  rank of all together:", qr(rbind(N, rows, E))$rank, "\n")

section("The inequality displayed in Section 4.3")
coefficient <- function(v, s) as.numeric(do.call(paste0, lev[c("V1", "V2", "V3")]) == v & do.call(paste0, lev[c("X", "Y", "Z")]) == s)
a <- coefficient("001", "010") + coefficient("010", "101") - coefficient("110", "001") - coefficient("010", "100") -
  coefficient("000", "111") - coefficient("001", "000") - coefficient("011", "000")
i <- which(H[, 1] == 0 & H[, 2] == 0 & apply(H[, -(1:2)], 1, \(r) all(-r == a)))
cat("Row of H_i:", i, "  nontrivial:", if (length(i)) fl[i] else NA, "\n")

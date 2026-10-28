## Section 4.4, Table 1: size, constraint counts and computation times.
## Usage: Rscript 10-benchmarks.R [names of cases to skip]
## e.g.   Rscript 10-benchmarks.R "Tripartite Bell"
## The redundancy check (rcdd::redundant) is slow: about 4 hours for Figure 7 and
## 7.5 hours for the ternary IV model. Set check_redundancy <- FALSE to skip it.
## The screening step times prob_simplex_screen(), the closed-form screening of
## Proposition 4 in meraconstraints (>= 0.3.0).
source("common.R")
check_redundancy <- FALSE

graphs <- list(
  "IV, |Z|=|X|=|Y|=2"      = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)),
  "IV, |Z|=3, |X|=|Y|=2"   = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3)),
  "IV, |Z|=4, |X|=|Y|=2"   = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 4)),
  "IV, |Z|=|X|=|Y|=3"      = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = 3, X = 3, Y = 3)),
  "Fig 5 right"            = dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)),
  "X -> A <- U -> B <- Y"  = dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)),
  "Fig 7 merged"           = dhvm_merge_latents(dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                                                          U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5))),
  "Tripartite Bell"        = dhvm(caugi(X %-->% A, Y %-->% B, Z %-->% C, U %-->% A + B + C))
)

## Optionally skip slow cases, e.g. Rscript benchmarks.R "Tripartite Bell"
skip <- commandArgs(trailingOnly = TRUE)
graphs <- graphs[setdiff(names(graphs), skip)]

res <- list()
for (nm in names(graphs)) {
  cat("Running:", nm, "\n"); flush.console()
  g   <- graphs[[nm]]
  rel <- dhvm_linear_relations(g)
  t_vh <- system.time(tope <- dhvm_relations_to_polytope(rel))["elapsed"]
  t_ir <- if (check_redundancy) system.time(Hirr <- lapply(tope$Hrep, \(h) rcdd::redundant(rcdd::d2q(h))$output |> rcdd::q2d()))["elapsed"] else { Hirr <- tope$Hrep; NA }
  t_sc <- system.time(flag <- prob_simplex_screen(rel, Hirr))["elapsed"]
  for (j in seq_along(Hirr)) {
    V <- tope$Vrep[[j]]
    res[[length(res) + 1]] <- data.frame(
      graph = nm, district = j,
      rows_B = ncol(V) - 2, cols_B = nrow(V),
      n_ineq = sum(Hirr[[j]][, 1] == 0), n_eq = sum(Hirr[[j]][, 1] == 1),
      n_redundant_removed = if (check_redundancy) nrow(tope$Hrep[[j]]) - nrow(Hirr[[j]]) else NA,
      n_flagged = sum(flag[[j]]),
      sec_VtoH = t_vh, sec_irredundant = t_ir, sec_screen = t_sc)
  }
}
res <- do.call(rbind, res)
print(res)
write.csv(res, "table1-benchmarks.csv", row.names = FALSE)

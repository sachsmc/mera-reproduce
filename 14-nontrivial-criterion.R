## Section 3.3, Proposition 6: which constraints are nontrivial, its relation to the nested
## Markov model, and why the equality constraints must be used.
source("common.R")

graphs <- list(
  "IV"                    = dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y)),
  "Front door"            = dhvm(caugi(A %-->% M %-->% Y, U1 %-->% A + Y)),
  "Fig 5 left"            = dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V1 %-->% V4, U2 %-->% V2 + V4)),
  "Fig 5 right"           = dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4)),
  "X -> A <- U -> B <- Y" = dhvm(caugi(X %-->% A, Y %-->% B, U %-->% A + B)),
  "Verma graph"           = dhvm(caugi(A %-->% B %-->% C %-->% D, U %-->% B + D)),
  "Fig 7 merged"          = dhvm_merge_latents(dhvm(caugi(V1 %-->% V2 + V6, V2 %-->% V3, V3 %-->% V4, V4 %-->% V5, V6 %-->% V3,
                                                          U1 %-->% V2 + V4, U2 %-->% V1 + V3 + V6, U3 %-->% V4 + V5)))
)
d2 <- meraconstraints:::district_levels; ir <- meraconstraints:::indicator_rows
cs <- meraconstraints:::closed_member_sets; mr <- meraconstraints:::margin_rows
rk <- function(m) if (is.null(m) || !nrow(m)) 0L else qr(m)$rank

## ---- counts, and the margin equations of the sets closed in the district ------------
## N and O are the equations of Proposition 6. M adds the sets that are closed in the
## district (they contain every member that is a parent of a member) but not closed.
section("Counts by district, and the equations that span the equality constraints")
rows <- list()
for (nm in names(graphs)) {
  g   <- graphs[[nm]]; cn <- meraconstraints(g); rel <- cn$functional_relations
  obs <- subset(g$nodes_data, latent == FALSE)$name
  anc <- lapply(setNames(obs, obs), \(v) intersect(ancestors(g$caugi_graph, nodes = v), obs))
  for (j in seq_along(cn$polytope$Hrep)) {
    H   <- cn$polytope$Hrep[[j]]; f <- cn$nontrivial[[j]]
    lev <- d2(rel, j); mem <- rel$district_parms$poterms[[j]]$out; ext <- rel$district_parms$poterms[[j]]$do
    pa  <- lapply(setNames(mem, mem), \(m) intersect(parents(g$caugi_graph, nodes = m), obs))
    ext_pa <- lapply(pa, intersect, ext)
    closed <- cs(mem, anc)
    in_district <- Filter(\(S) all(intersect(unlist(pa[S]), mem) %in% S), cs(mem, lapply(anc, \(a) character())))
    extra <- Filter(\(S) !any(vapply(closed, setequal, logical(1), S)), in_district)
    N <- ir(lev, ext)
    O <- do.call(rbind, lapply(closed, mr, lev = lev, ext = ext, ext_pa = ext_pa))
    M <- do.call(rbind, lapply(extra, mr, lev = lev, ext = ext, ext_pa = ext_pa))
    E <- -H[H[, 1] == 1, -(1:2), drop = FALSE]
    ## Proposition 6(ii): trivial inequalities are exactly the nonnegativity constraints
    k <- nonnegativity_entries(cn, j)
    stopifnot(identical(is.na(k), unname(f[H[, 1] == 0])), length(unique(k[!is.na(k)])) == sum(!is.na(k)))
    rows[[length(rows) + 1]] <- data.frame(graph = nm, district = names(cn$polytope$Hrep)[j],
      ineq = sum(H[, 1] == 0), nontrivial_ineq = sum(f & H[, 1] == 0), eq = nrow(E), nontrivial_eq = sum(f & H[, 1] == 1),
      rank_NO = rk(rbind(N, O)), rank_He = rk(rbind(N, O, E)), rank_NOM = rk(rbind(N, O, M)),
      closed_in_district_only = paste(vapply(extra, paste, "", collapse = ","), collapse = "; "))
  }
}
print(do.call(rbind, rows), row.names = FALSE)

## ---- the row used as an example in Section 3.3 ---------------------------------------
section("X -> A <- U -> B <- Y: the row p(00 | 11) <= p(00 | 10) + p(01 | 10)")
g <- graphs[["X -> A <- U -> B <- Y"]]; cn <- meraconstraints(g)
H <- cn$polytope$Hrep[[1]]; lev <- entry_levels(cn, 1)
p <- function(a, b, x, y) with(lev, as.numeric(A == a & B == b & X == x & Y == y))
h <- p(0, 0, 1, 1) - p(0, 0, 1, 0) - p(0, 1, 1, 0)
i <- which(H[, 1] == 0 & H[, 2] == 0 & apply(H[, -(1:2)], 1, \(r) all(-r == h)))
cat("Row of the H-representation:", i, "  nontrivial:", cn$nontrivial[[1]][i], "\n")
V <- rcdd::q2d(cn$polytope$Vrep[[1]])[, -(1:2)]
on <- V[as.vector(V %*% h) == 0, , drop = FALSE]
cat("Vertices on its facet:", nrow(on), "  of which p(01 | 11) = 0:", sum(on %*% p(0, 1, 1, 1) == 0), "\n")
cat("Maximum over the product of simplices (bound 0):",
    sum(tapply(h, paste(lev$X, lev$Y), max)), "\n")
E <- H[H[, 1] == 1, -1, drop = FALSE]
cat("Equal to p(01 | 11) >= 0 after adding equality constraints:",
    qr(rbind(E, c(0, h + p(0, 1, 1, 1))))$rank == qr(E)$rank, "\n")

## ---- invariance: rewriting rows changes the old screen, not Proposition 6 -------------
section("Adding random combinations of the equality rows to the inequality rows")
set.seed(1)
rel  <- cn$functional_relations; tope <- cn$polytope
eq   <- H[, 1] == 1
H2   <- H
H2[!eq, -1] <- H[!eq, -1] + matrix(sample(-2:2, sum(!eq) * sum(eq), TRUE), sum(!eq)) %*% H[eq, -1]
tope2 <- tope; tope2$Hrep[[1]] <- H2; tope2$HrepQ[[1]] <- rcdd::d2q(H2)
cat("Inequalities flagged by the simplex screen, as returned:", sum(prob_simplex_screen(rel, tope$Hrep)[[1]][!eq]),
    "  rewritten:", sum(prob_simplex_screen(rel, tope2$Hrep)[[1]][!eq]), "\n")
cat("Nontrivial by Proposition 6, as returned:", sum(cn$nontrivial[[1]][!eq]),
    "  rewritten:", sum(nontrivial_screen(g, rel, tope2)[[1]][!eq]), "\n")

## ---- linear programming check of Proposition 6(ii) ------------------------------------
section("Maximizing each inequality over the simplices and the equality constraints")
for (nm in setdiff(names(graphs), "Fig 7 merged")) {
  cn <- meraconstraints(graphs[[nm]])
  agree <- sapply(seq_along(cn$polytope$Hrep), \(j) {
    H <- cn$polytope$Hrep[[j]]; eq <- H[, 1] == 1
    Q <- rbind(H[eq, , drop = FALSE], cbind(0, 0, diag(ncol(H) - 2)))
    lp <- vapply(which(!eq), \(r) {
      s <- rcdd::lpcdd(rcdd::d2q(Q), rcdd::d2q(-H[r, -(1:2)]), minimize = FALSE)
      rcdd::q2d(s$optimal.value) > H[r, 2] + 1e-9
    }, logical(1))
    identical(unname(cn$nontrivial[[j]][!eq]), lp)
  })
  cat(sprintf("%-24s agrees in every district: %s\n", nm, all(agree)))
}

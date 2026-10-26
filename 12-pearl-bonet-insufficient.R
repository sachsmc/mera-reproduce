## Kedagni and Mourifie (2020) state that Pearl's instrumental inequality and Bonet's (2001)
## inequality are not sufficient to detect all violations of the instrumental variable
## assumptions (their numerical illustration is in their Supplementary Material). Here
## "Bonet's inequality" is a single inequality, their (10) with y = d = 0, in one labelling
## of the outcome and treatment values. This script constructs such an example for binary
## treatment D (= X), binary outcome Y and a three-level instrument Z, and checks it exactly
## against the complete constraints. It then shows that Pearl's inequality together with
## all relabellings of Bonet's inequality (all four versions of (10)) is complete, as Bonet
## showed using the same polytope method.
##
## Notation: q_yd(z) = P(Y = y, D = d | Z = z) = P*(X = d, Y = y | Z = z).
source("common.R")

nz <- 3
g  <- dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Z = nz))
cn <- meraconstraints(g)
lev <- cn$functional_relations$problevels[[1]][, c("X", "Y", "Z")]
n   <- nrow(lev); Zs <- 0:(nz - 1)
idx <- function(d, y, z) which(lev$X == d & lev$Y == y & lev$Z == z)
q   <- function(y, d, z, s = 1) { a <- numeric(n); a[idx(d, y, z)] <- s; a }
grid <- function(k) as.list(as.data.frame(t(expand.grid(rep(list(Zs), k)))))

## Each family returns rows c(b, a) meaning a . q <= b; every sup or inf over z becomes
## one inequality for each choice of instrument values.
pearl <- function() do.call(rbind, lapply(0:1, \(d) do.call(rbind, lapply(grid(2), \(z)
  c(1, q(0, d, z[1]) + q(1, d, z[2]))))))
## (10): q_yd(z) <= inf_w [q_yd(w) + q_00-or-.. ] ... in the form
##   d = 1: q_y1(z) <= [q_00(w1) + q_y1(w1)] + [q_10(w2) + q_y1(w2)]
##   d = 0: q_y0(z) <= [q_y0(w1) + q_01(w1)] + [q_y0(w2) + q_11(w2)]
ineq10 <- function(y, d) do.call(rbind, lapply(grid(3), \(z) {
  w1 <- z[2]; w2 <- z[3]
  rhs <- if (d == 1) q(0, 0, w1) + q(y, 1, w1) + q(1, 0, w2) + q(y, 1, w2)
         else        q(y, 0, w1) + q(0, 1, w1) + q(y, 0, w2) + q(1, 1, w2)
  c(0, q(y, d, z[1]) - rhs)
}))
bonet <- function() ineq10(0, 0)                    # q_00(z) <= q_0(w1) + q_00+11(w2)
## (11): 1 <= sum_{i,j} [q_i0(w_ij) + q_j1(w_ij)]
ineq11 <- function() do.call(rbind, lapply(grid(4), \(w) {
  a <- numeric(n); k <- 0
  for (i in 0:1) for (j in 0:1) { k <- k + 1; a <- a - q(i, 0, w[k]) - q(j, 1, w[k]) }
  c(-1, a)
}))
ineq10_all <- function() rbind(ineq10(0, 0), ineq10(1, 0), ineq10(0, 1), ineq10(1, 1))

polytope_vertices <- function(R) {
  simplex <- t(sapply(Zs, \(z) as.numeric(lev$Z == z)))
  H <- rbind(cbind(0, R[, 1], -R[, -1]), cbind(0, 0, diag(n)), cbind(1, -1, simplex))
  V <- rcdd::scdd(rcdd::d2q(unique(H)), representation = "H")$output
  rcdd::q2d(V[V[, 1] == "0" & V[, 2] == "1", -(1:2), drop = FALSE])
}
holds <- function(R, p) all(R[, -1, drop = FALSE] %*% p <= R[, 1] + 1e-12)
H_our <- cn$polytope$Hrep[[1]]
in_model <- function(p) hrep_test(H_our, p)

V_pb <- polytope_vertices(rbind(pearl(), bonet()))
outside <- which(!apply(V_pb, 1, in_model))
cat("Vertices of the Pearl + Bonet polytope:", nrow(V_pb), "\n")
cat("  of which incompatible with the instrumental variable model:", length(outside), "\n")

## Turn one such vertex into a strictly positive distribution by mixing it with the
## uniform distribution (q_yd(z) = 1/4), which lies inside every polytope here.
v <- V_pb[outside[1], ]
u <- rep(1/4, n)
p <- 0.9 * v + 0.1 * u
cat("\nExample: q = 0.9 * vertex + 0.1 * uniform\n")
tab <- data.table(lev)[, q := as.character(MASS::fractions(p))]
print(dcast(tab, Z ~ paste0("q_", Y, X), value.var = "q"))

cat("\nPearl's inequality (9) holds:            ", holds(pearl(), p), "\n")
cat("Bonet's inequality holds:                ", holds(bonet(), p), "\n")
for (yd in list(c(0, 0), c(1, 0), c(0, 1), c(1, 1)))
  cat(sprintf("Inequality (10) with y = %d, d = %d holds:  ", yd[1], yd[2]), holds(ineq10(yd[1], yd[2]), p), "\n")
cat("Inequality (11) holds:                   ", holds(ineq11(), p), "\n")
cat("Compatible with the IV model (complete): ", in_model(p), "\n")

viol <- !hrep_test(H_our, p, which = TRUE)
section("Complete constraints violated by the example, in terms of P*")
cat(hrep_to_text(H_our[viol, , drop = FALSE], cn$functional_relations$probstarform[[1]]), sep = "\n")

section("Which sets of inequalities are complete? (|Z| = 3)")
V_our <- unique(rcdd::q2d(cn$polytope$Vrep[[1]])[, -(1:2)])
key <- function(M) apply(round(M, 10), 1, paste, collapse = ",")
compare_set <- function(R, label) {
  V <- polytope_vertices(R)
  cat(sprintf("%-40s vertices %3d, outside the model %3d, identical to the complete set: %s\n",
              label, nrow(V), sum(!apply(V, 1, in_model)), setequal(key(V), key(V_our))))
}
compare_set(rbind(pearl(), bonet()),                "Pearl + Bonet (y = d = 0 only)")
compare_set(rbind(pearl(), ineq10_all()),           "Pearl + all relabellings of Bonet (10)")
compare_set(rbind(pearl(), ineq10_all(), ineq11()), "Pearl + (10) + (11)")

## Comparison with Kedagni and Mourifie (2020), "Generalized instrumental inequalities".
##
## Setting: binary treatment D (= X in the software), discrete outcome Y and instrument Z
## with finitely many values, and the exclusion restriction. Kedagni and Mourifie (KM)
## derive inequalities (12)-(14), one family for each partition of the outcome's values,
## that must hold under joint statistical independence (Y0, Y1) _||_ Z. For binary Y they
## prove these are also sufficient (their Theorem 1).
##
## Our method gives the complete set of constraints for the instrumental variable graph,
## which corresponds to the stronger random assignment assumption {Y0, Y1, (D_z)} _||_ Z.
## Both sets of constraints are linear in q_yd(z) = P(Y = y, D = d | Z = z), which is
## our P*(X = d, Y = y | Z = z). So both define polytopes in the same space, and we can
## compare them exactly by computing their vertices (in exact rational arithmetic).
source("common.R")

## All set partitions of a vector
set_partitions <- function(x) {
  if (length(x) == 1) return(list(list(x)))
  out <- list()
  for (p in set_partitions(x[-1])) {
    for (i in seq_along(p)) { q <- p; q[[i]] <- c(x[1], q[[i]]); out[[length(out) + 1]] <- q }
    out[[length(out) + 1]] <- c(list(x[1]), p)
  }
  out
}

km_compare <- function(ny, nz) {
  g  <- dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y), nvals = c(Y = ny, Z = nz))
  cn <- meraconstraints(g)
  lev <- cn$functional_relations$problevels[[1]][, c("X", "Y", "Z")]
  n  <- nrow(lev)
  idx <- function(d, y, z) which(lev$X == d & lev$Y == y & lev$Z == z)
  Y <- 0:(ny - 1); Z <- 0:(nz - 1)

  rows <- list(); rhs <- c()
  add <- function(a, b) { rows[[length(rows) + 1]] <<- a; rhs <<- c(rhs, b) }
  qA <- function(a, A, d, z, s = 1) { for (y in A) a[idx(d, y, z)] <- a[idx(d, y, z)] + s; a }

  for (d in 0:1) {                       # (12): sum_y sup_z q_yd(z) <= 1
    for (zs in as.list(as.data.frame(t(expand.grid(rep(list(Z), ny)))))) {
      a <- numeric(n); for (k in seq_along(Y)) a <- qA(a, Y[k], d, zs[k]); add(a, 1)
    }
  }
  for (P in set_partitions(Y)) {
    np <- length(P)
    for (d in 0:1) for (A in P) {       # (13): sum_{y in A} sup_z q_yd(z) <= sum_B inf_z [q_Ad(z) + q_B,1-d(z)]
      for (zs in as.list(as.data.frame(t(expand.grid(rep(list(Z), length(A) + np)))))) {
        a <- numeric(n)
        for (k in seq_along(A)) a <- qA(a, A[k], d, zs[k])
        for (b in seq_len(np)) {
          w <- zs[length(A) + b]
          a <- qA(a, A, d, w, -1); a <- qA(a, P[[b]], 1 - d, w, -1)
        }
        add(a, 0)
      }
    }
    pairs <- expand.grid(i = seq_len(np), j = seq_len(np))   # (14): 1 <= sum_{A,B} inf_z [q_A0(z) + q_B1(z)]
    for (zs in as.list(as.data.frame(t(expand.grid(rep(list(Z), nrow(pairs))))))) {
      a <- numeric(n)
      for (k in seq_len(nrow(pairs))) {
        a <- qA(a, P[[pairs$i[k]]], 0, zs[k], -1); a <- qA(a, P[[pairs$j[k]]], 1, zs[k], -1)
      }
      add(a, -1)
    }
  }
  A_km <- unique(cbind(rhs, do.call(rbind, rows)))
  message(sprintf("[%s] KM inequalities built: %d distinct", format(Sys.time(), "%H:%M:%S"), nrow(A_km)))

  ## KM polytope = KM inequalities + probability simplices, in rcdd H format [l, b, -A]
  simplex_eq <- t(sapply(Z, \(z) as.numeric(lev$Z == z)))
  H_km <- rbind(cbind(0, A_km[, 1], -A_km[, -1]),            # a q <= b
                cbind(0, 0, diag(n)),                          # q >= 0
                cbind(1, -1, simplex_eq))                      # sum_{d,y} q(z) = 1
  message(sprintf("[%s] computing vertices of KM polytope", format(Sys.time(), "%H:%M:%S")))
  V_km <- rcdd::scdd(rcdd::d2q(H_km), representation = "H")$output
  message(sprintf("[%s] done", format(Sys.time(), "%H:%M:%S")))
  V_km <- V_km[V_km[, 1] == "0" & V_km[, 2] == "1", -(1:2), drop = FALSE]   # vertices
  V_our <- rcdd::q2d(cn$polytope$Vrep[[1]])
  V_our <- unique(V_our[, -(1:2), drop = FALSE])
  key <- function(M) apply(M, 1, paste, collapse = ",")
  k_km  <- key(rcdd::q2d(V_km)); k_our <- key(V_our)

  ## Our facets that some KM vertex violates (constraints missing from KM's set)
  H <- cn$polytope$Hrep[[1]]
  Vk <- rcdd::q2d(V_km)
  lhs <- Vk %*% t(-H[, -(1:2)])
  viol <- which(apply(sweep(lhs, 2, H[, 2], "-") > 1e-9, 2, any))

  list(ny = ny, nz = nz, n_km_ineq = nrow(A_km), n_our_ineq = sum(H[, 1] == 0),
       v_km = length(k_km), v_our = length(k_our),
       km_valid = all(sweep(V_our %*% t(A_km[, -1]), 2, A_km[, 1], "-") <= 1e-9),
       equal = setequal(k_km, k_our), extra_vertices = Vk[!(k_km %in% k_our), , drop = FALSE],
       missing = if (length(viol)) hrep_to_text(H[viol, , drop = FALSE], cn$functional_relations$probstarform[[1]]) else character(),
       lev = lev)
}

## Binary outcome with binary and ternary instruments (a few seconds). Larger cases,
## e.g. c(2, 4) or a ternary outcome c(3, 2), take much longer.
for (case in list(c(2, 2), c(2, 3))) {
  section(sprintf("|Y| = %d, |D| = 2, |Z| = %d", case[1], case[2]))
  t <- system.time(r <- km_compare(case[1], case[2]))["elapsed"]
  cat(sprintf("KM inequalities (distinct): %d   our inequality constraints: %d\n", r$n_km_ineq, r$n_our_ineq))
  cat(sprintf("Vertices of KM polytope: %d   vertices of our polytope: %d\n", r$v_km, r$v_our))
  cat("All KM inequalities hold on our polytope (KM valid for the graph):", r$km_valid, "\n")
  cat("Polytopes identical:", r$equal, sprintf("  (%.1f s)\n", t))
  if (!r$equal) {
    cat("Our constraints violated by some vertex of the KM polytope:", length(r$missing), "\n")
    cat(head(r$missing, 5), sep = "\n")
    cat("One KM vertex outside our polytope, P*(X, Y | Z) in the order X, Y, Z:\n")
    print(cbind(r$lev, q = MASS::fractions(r$extra_vertices[1, ])))
  }
}

## Section 4.1: two nested Markov equivalent graphs with different inequality constraints
## (Figure 5). Both graphs have the same conditional independence and the same, trivial,
## equality constraints. The 8 nontrivial inequality constraints of the right graph are the
## instrumental inequalities for each value of V3. Distributions sampled from the left graph
## never violate the left graph's constraints, but do violate these 8.
source("common.R")
set.seed(2024)

left  <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V1 %-->% V4, U2 %-->% V2 + V4))
right <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4))
cn_left  <- meraconstraints(left)
cn_right <- meraconstraints(right)
print(cn_left)
print(cn_right)

j  <- "V2_V4"
Hl <- cn_left$polytope$Hrep[[j]]
Hr <- cn_right$polytope$Hrep[[j]]
lev <- entry_levels(cn_right, j)
stopifnot(identical(do.call(paste, lev), do.call(paste, entry_levels(cn_left, j)[, names(lev)])))

section("Equality constraints")
El <- Hl[Hl[, 1] == 1, -1]; Er <- Hr[Hr[, 1] == 1, -1]
cat("Rank left:", qr(El)$rank, "  right:", qr(Er)$rank, "  together:", qr(rbind(El, Er))$rank, "\n")
cat("Nontrivial, left:", sum(cn_left$nontrivial[[j]] & Hl[, 1] == 1),
    "  right:", sum(cn_right$nontrivial[[j]] & Hr[, 1] == 1), "\n")

section("Inequality constraints")
cat("Left:", sum(Hl[, 1] == 0), "of which nontrivial", sum(cn_left$nontrivial[[j]] & Hl[, 1] == 0), "\n")
cat("Right:", sum(Hr[, 1] == 0), "of which nontrivial", sum(cn_right$nontrivial[[j]] & Hr[, 1] == 0), "\n")

## The form in the paper:
## P*(v2, v4 | V1 = 0, v3) + P*(v2, 1 - v4 | V1 = 1, v3) <= 1, for v2, v3, v4 in {0, 1}
fam <- list(); lab <- character()
for (v3 in 0:1) for (v2 in 0:1) for (v4 in 0:1) {
  fam[[length(fam) + 1]] <- with(lev, as.numeric(V2 == v2 & V3 == v3 & ((V4 == v4 & V1 == 0) | (V4 == 1 - v4 & V1 == 1))))
  lab <- c(lab, sprintf("v2 = %d, v3 = %d, v4 = %d", v2, v3, v4))
}
fam <- do.call(rbind, fam)

## Each returned nontrivial row must equal one of these rows after adding equality constraints
section("Each nontrivial row of the right graph, and its form in the paper")
nt   <- which(cn_right$nontrivial[[j]] & Hr[, 1] == 0)
rk   <- qr(Er)$rank
pbar <- colMeans(rcdd::q2d(cn_right$polytope$Vrep[[j]])[, -(1:2)])
same <- function(row, f) {
  fr <- c(1, -f)
  qr(rbind(Er, row, fr))$rank == rk + 1 && sum(row * c(1, pbar)) * sum(fr * c(1, pbar)) > 0
}
hit <- sapply(nt, \(r) which(apply(fam, 1, \(f) same(Hr[r, -1], f))))
stopifnot(is.numeric(hit), setequal(hit, seq_len(nrow(fam))))
txt <- hrep_to_text(Hr, cn_right$functional_relations$probstarform[[j]])[nt]
cat(sprintf("%s  <-  %s", lab[hit], txt), sep = "\n")

Vr <- rcdd::q2d(cn_right$polytope$Vrep[[j]])[, -(1:2)]
Vl <- rcdd::q2d(cn_left$polytope$Vrep[[j]])[, -(1:2)]
cat("\nLargest left side of these inequalities over the vertices, right graph:", max(fam %*% t(Vr)),
    "  left graph:", max(fam %*% t(Vl)), "\n")

section("Example: V2 = 0 with probability 0.9, V4 = V1 with probability 0.9")
p <- with(lev, ifelse(V2 == 0, .9, .1) * ifelse(V4 == V1, .9, .1))
cat("All constraints of the left graph hold:", hrep_test(Hl, p), "\n")
cat("Rows of the right graph violated:", sum(!hrep_test(Hr, p, which = TRUE)), "\n")
print(setNames(as.vector(fam %*% p), lab))

## P*(V2, V4 | V1, V3) = P(V2 | V1) P(V4 | V1, V2, V3), computed from a sampled distribution
pstar_from <- function(probs, cn) {
  probs <- copy(probs)
  probs[, prob.V1V2V3 := prob.V3_V2 * prob.V2_V1 * prob.V1]
  probs[, prob.V4_rest := prob.joint / prob.V1V2V3]
  probs[, prob.test := prob.V4_rest * prob.V2_V1]
  as_entries(cn, j, probs, "prob.test")
}

nsim <- 500
viol_right <- matrix(FALSE, nsim, nrow(Hr))
viol_left  <- matrix(FALSE, nsim, nrow(Hl))
for (i in 1:nsim) {
  probs <- data.table(dhvm_sample_distribution(left, parmsamp = \(k) rnorm(k, sd = 5)))
  viol_right[i, ] <- !hrep_test(Hr, pstar_from(probs, cn_right), which = TRUE)
  viol_left[i, ]  <- !hrep_test(Hl, pstar_from(probs, cn_left), which = TRUE)
}
section("Sampling from the left graph")
cat("Constraints of the left graph ever violated:", sum(colSums(viol_left) > 0), "\n")
cat("Constraints of the right graph violated:", sum(colSums(viol_right) > 0), "\n")
cat("These are exactly the nontrivial rows:", setequal(which(colSums(viol_right) > 0), nt), "\n")

## Section 4.2: two nested Markov equivalent graphs with different inequality
## constraints (Figure 5). Distributions sampled from the left graph never violate the
## left graph's constraints, but do violate 8 inequality constraints of the right graph.
source("common.R")
set.seed(2024)

left  <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V1 %-->% V4, U2 %-->% V2 + V4))
right <- dhvm(caugi(V1 %-->% V2 %-->% V3 %-->% V4, V2 %-->% V4, U2 %-->% V2 + V4))
cn_left  <- meraconstraints(left)
cn_right <- meraconstraints(right)
print(cn_left)
print(cn_right)

## P*(V2, V4 | V1, V3) = P(V2 | V1) P(V4 | V1, V2, V3), computed from a sampled distribution
pstar_from <- function(probs, cn) {
  probs <- copy(probs)
  probs[, prob.V1V2V3 := prob.V3_V2 * prob.V2_V1 * prob.V1]
  probs[, prob.V4_rest := prob.joint / prob.V1V2V3]
  probs[, prob.test := prob.V4_rest * prob.V2_V1]
  merge(cn$functional_relations$problevels[["V2_V4"]], probs, sort = FALSE)$prob.test
}

nsim <- 500
viol_right <- matrix(FALSE, nsim, nrow(cn_right$polytope$Hrep[["V2_V4"]]))
viol_left  <- matrix(FALSE, nsim, nrow(cn_left$polytope$Hrep[["V2_V4"]]))
for (i in 1:nsim) {
  probs <- data.table(dhvm_sample_distribution(left, parmsamp = \(k) rnorm(k, sd = 5)))
  p <- pstar_from(probs, cn_right)
  viol_right[i, ] <- !hrep_test(cn_right$polytope$Hrep[["V2_V4"]], p, which = TRUE)
  viol_left[i, ]  <- !hrep_test(cn_left$polytope$Hrep[["V2_V4"]],  p, which = TRUE)
}
cat("\nConstraints of the left graph ever violated by samples from the left graph:",
    sum(colSums(viol_left) > 0), "\n")
cat("Constraints of the right graph violated by samples from the left graph:",
    sum(colSums(viol_right) > 0), "\n")
section("Those constraints of the right graph, in terms of interventional probabilities")
cat(hrep_to_text(cn_right$polytope$Hrep[["V2_V4"]],
                 cn_right$functional_relations$probstarform[["V2_V4"]])[colSums(viol_right) > 0], sep = "\n")

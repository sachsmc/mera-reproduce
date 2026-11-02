## Section 3 and Appendix B: the instrumental variable model (Figure 4).
source("common.R")

iv <- dhvm(caugi(Z %-->% X, X %-->% Y, U %-->% X + Y))
cn <- meraconstraints(iv)
print(cn)

section("Matrix B (rows: P*(x, y | z); columns: joint values of the response function variables)")
rel <- cn$functional_relations$relation[[1]]
r   <- unique(unlist(rel))
B   <- t(sapply(rel, \(rr) as.integer(r %in% rr)))
print(unname(B))
cat("Dimensions:", dim(B), "  distinct columns:", ncol(unique(B, MARGIN = 2)), "\n")

section("H-representation (rcdd format: columns are [equality flag, b, -A])")
print(cn$polytope$Hrep[[1]])

section("All constraints, in terms of observed probabilities")
cat(constraint_text(cn, 1, nontrivial_only = FALSE), sep = "\n")

section("Nontrivial rows of H_i (Proposition 6)")
H  <- cn$polytope$Hrep[[1]]
Hi <- which(H[, 1] == 0)
cat("Rows of H_i:", which(cn$nontrivial[[1]][Hi]), "\n")
cat(constraint_text(cn, 1), sep = "\n")

section("The other rows of H_i are nonnegativity constraints")
k <- nonnegativity_entries(cn, 1)
keep <- !is.na(k)
print(data.frame(row = which(keep), entry_levels(cn, 1)[k[keep], c("X", "Y", "Z")], row.names = NULL))

section("The code in Section 6")
print(hrep_test(H, c(0.101, 0.338, 0.013, 0.548, 0.317, 0.51, 0.004, 0.169), which = TRUE))
